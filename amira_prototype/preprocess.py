import os
import pysam
from tqdm import tqdm
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from construct_unitig import parse_fastq

def parse_sam_file(sam_file_path,
                genesOfInterest):
    """
    Reads a SAM file and returns a dictionary of read ID to query locations and sequences.
    sam_file_path: str, path to SAM file to parse
    """
    sam_file = pysam.AlignmentFile(sam_file_path, "r")
    # Initialize an empty dictionary to hold the parsed results
    read_data = {}
    # Iterate over each read in the SAM file
    for read in sam_file:
        if read.is_mapped:
            read_id = read.query_name
            #if read.reference_name.replace(".aln.fas", "").replace(".fasta", "") in genesOfInterest:
            if read.cigartuples[0][0] == 5:
                start_bp = read.cigartuples[0][1]
            else:
                start_bp = 0
            end_bp = start_bp + read.query_alignment_length
            sequence = read.query_alignment_sequence
            strand = '-' if read.is_reverse else '+'
            # Create a dictionary of the read's data
            read_dict = {"sequence": sequence, "strand": strand, "identifier": read.reference_name.replace(".aln.fas", "").replace(".fasta", "")}
            # Add the read's data to the dictionary of read data
            if read_id not in read_data:
                read_data[read_id] = {}
            read_data[read_id][(start_bp, end_bp)] = read_dict
    sam_file.close()
    to_delete = []
    for readId in read_data:
        positions = read_data[readId]
        if not any(positions[p]["identifier"] in genesOfInterest for p in positions):
            to_delete.append(readId)
    for r in to_delete:
        del read_data[r]
    return read_data

def supplement_read_data(read_amr_genes, fastq_content, interval_length):
    supplemented_read_data = {}
    for readId in tqdm(read_amr_genes):
        positions = list(read_amr_genes[readId].keys())
        sorted_positions = sorted(positions, key=lambda x: x[0])
        readLength = len(fastq_content[readId]["sequence"])
        new_positions = []
        #if not sorted_positions[0][0] == 0:
        #    new_positions += split_interval(sorted_positions[0][0] - 1 - interval_length, sorted_positions[0][0] - 1, interval_length)
        #     if not (sorted_positions[0][0] - 1 > interval_length and sorted_positions[0][0] - 1 < interval_length * 2):
        #         new_positions += split_interval(0, sorted_positions[0][0] - 1, interval_length)
        #     else:
        #         new_positions += split_interval(sorted_positions[0][0] - 1 - interval_length, sorted_positions[0][0] - 1, interval_length)
        #if not sorted_positions[-1][1] == readLength - 1:
        #    new_positions += split_interval(sorted_positions[-1][1] + 1, sorted_positions[-1][1] + 1 + interval_length, interval_length)
        #     loci_length = readLength - 1 - sorted_positions[-1][1] + 1
        #     if not (loci_length > interval_length and loci_length < interval_length * 2):
        #         new_positions += split_interval(sorted_positions[-1][1] + 1, readLength - 1, interval_length)
        #     else:
        #         new_positions += split_interval(sorted_positions[-1][1] + 1, sorted_positions[-1][1] + 1 + interval_length, interval_length)
        for p in range(len(sorted_positions)-1):
            if not sorted_positions[p][1] >= sorted_positions[p+1][0]:
                new_positions += split_interval(sorted_positions[p][1] + 1, sorted_positions[p+1][0] - 1, interval_length)
        supplemented_read_data[readId] = {}
        for pos in new_positions:
            if pos:
                supplemented_read_data[readId][pos] = {"sequence": fastq_content[readId]["sequence"][pos[0]: pos[1]]}
            #else:
             #   supplemented_read_data[read_id][pos] = read_amr[pos]
    return supplemented_read_data

import numpy as np

def split_interval(start, end, sub_interval_length):
    # If the length of the interval is less than or equal to the sub-interval length, just return the original interval
    interval_length = end - start + 1
    if interval_length == sub_interval_length:
        return [(start, end)]
    if interval_length < sub_interval_length:
        return [None]
    n = round((end - start + 1) / sub_interval_length)
    all_bases = np.array([i for i in range(start, end + 2)])
    chunks = np.array_split(all_bases, n)
    sub_intervals = [(c[0], c[-1] - 1) for c in chunks]
    return sub_intervals

def mmseqs2_easy_cluster(supplemented_read_data,
                        threads,
                        outputDir):
    # Write sequences to input file
    seq_file_path = os.path.join(outputDir, "supplemented_read_data.fasta")
    with open(seq_file_path, "w") as seq_file:
        for read_id, sub_intervals in supplemented_read_data.items():
            for pos, data in sub_intervals.items():
                seq_file.write(f">{read_id}_{pos[0]}_{pos[1]}\n{data['sequence']}\n")
    # Run mmseqs2 easy-cluster
    cluster_file_path = os.path.join(outputDir, "supplemented_read_data")
    subprocess.run(["/home/daniel/Documents/GitHub/amira_prototype/panRG_building_tools/MMseqs2/build/bin/mmseqs", "easy-cluster", seq_file_path, cluster_file_path, "tmp", "--threads", str(threads), "-c", "0.7", "--cov-mode", "0", "--kmer-per-seq", "100", "--cluster-reassign", "1"],
                check=True)

def map_sequences_to_representatives(input_fasta, cluster_tsv, minimap2_output, threads, outputDir):
    # Parse the cluster file
    sequence_to_representative = {}
    representatives = set()
    with open(cluster_tsv, "r") as cluster_file:
        for line in cluster_file:
            cluster_representative, sequence_id = line.strip().split("\t")
            sequence_to_representative[sequence_id] = cluster_representative
            representatives.add(cluster_representative)

    # Write a new FASTA file with only the cluster representative sequences
    representative_fasta = os.path.join(outputDir, "representative_sequences.fasta")
    with open(representative_fasta, "w") as rep_fasta_file:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in representatives:
                record.seq = Seq("N" * 50) + record.seq + Seq("N" * 50)
                SeqIO.write(record, rep_fasta_file, "fasta")

    # Run minimap2 to map the input sequences to the cluster representatives
    minimap2_command = [
        "minimap2",
        "-a",  # Output in SAM format
        "-t",
        str(threads),
        representative_fasta,  # Reference file with cluster representatives
        input_fasta,  # Query file with input sequences
        "-o", minimap2_output,  # Output SAM file
    ]
    subprocess.run(minimap2_command, check=True)

def determine_cluster_strands(minimap2_output, cluster_file_path):
    # import clustering information
    clustered = {}
    cluster_mapping = {}
    clusterId = 0
    with open(cluster_file_path, "r") as cluster_file:
        for line in cluster_file:
            parts = line.strip().split("\t")
            cluster_name = parts[0]
            if not cluster_name in clustered:
                clustered[cluster_name] = set()
                cluster_mapping[cluster_name] = clusterId
                clusterId += 1
            clustered[cluster_name].add(parts[1])
    # import the minimap2 mapping
    sam_content = pysam.AlignmentFile(minimap2_output, "r")
    # store the cluster and strand information
    clusterDict = {}
    for query in tqdm(sam_content.fetch()):
        if query.is_mapped:
            if query.query_name in clustered[query.reference_name]:
                if query.is_reverse:
                    strand = "-"
                else:
                    strand = "+"
                readId, pos_start, pos_end = query.query_name.split("_")
                if not readId in clusterDict:
                    clusterDict[readId] = {}
                clusterDict[readId][(int(pos_start), int(pos_end))] = strand + str(cluster_mapping[query.reference_name])
    return clusterDict

def qc_readfile(readfile,
                filtered_readfile):
    filtlong_command = "panRG_building_tools/Filtlong/bin/filtlong --min_length 1000 --keep_percent 95 " + readfile + " | gzip > " + filtered_readfile
    subprocess.run(filtlong_command,
                shell=True,
                check=True)

def process_sam(samfile,
                readfile,
                outputDir,
                threads,
                genesOfInterest):
    read_amr_genes = parse_sam_file(samfile,
                                    genesOfInterest)
    filtered_readfile = os.path.join(os.path.dirname(readfile), os.path.basename(readfile).replace(".fastq", ".filtered.fastq"))
    qc_readfile(readfile,
                filtered_readfile)
    fastq_content = parse_fastq(filtered_readfile)
    to_delete = []
    for r in read_amr_genes:
        if not r in fastq_content:
            to_delete.append(r)
    for r in to_delete:
        del read_amr_genes[r]
    supplemented_read_data = supplement_read_data(read_amr_genes,
                                                fastq_content,
                                                1000)
    mmseqs2_easy_cluster(supplemented_read_data,
                        threads,
                        outputDir)
    map_sequences_to_representatives(os.path.join(outputDir, "supplemented_read_data.fasta"),
                            os.path.join(outputDir, "supplemented_read_data_cluster.tsv"),
                            os.path.join(outputDir, "minimap2_output.sam"),
                            threads,
                            outputDir)
    clusterDict = determine_cluster_strands(os.path.join(outputDir, "minimap2_output.sam"),
                                os.path.join(outputDir, "supplemented_read_data_cluster.tsv"))
    readGenes = {}
    for readId in tqdm(read_amr_genes):
        amr_loci = read_amr_genes[readId]
        all_loci = []
        for p in amr_loci:
            all_loci.append((p[0], amr_loci[p]["strand"] + amr_loci[p]["identifier"]))
        if readId in clusterDict:
            cluster_loci = clusterDict[readId]
            for p in cluster_loci:
                all_loci.append((p[0], cluster_loci[p]))
        loci_list = [l[1] for l in sorted(all_loci, key=lambda x: x[0])]
        readGenes[readId] = loci_list
    return readGenes