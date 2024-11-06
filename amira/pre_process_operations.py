import json
import os
import subprocess
import sys

import pysam
from tqdm import tqdm


def clean_gene(g):
    chars_to_remove = set(["|", "(", ")", "-", "*", "+", "#", ":", "=", "/", ",", "'"])
    cleaned_gene = "".join(char for char in g if char not in chars_to_remove)
    return cleaned_gene


def process_pandora_json(
    pandoraJSON: str, genesOfInterest: list[str], gene_positions: str
) -> tuple[dict[str, list[str]], list[str]]:
    with open(pandoraJSON) as i:
        annotatedReads = json.loads(i.read())
    with open(gene_positions) as i:
        gene_position_dict = json.loads(i.read())
    to_delete = []
    subsettedGenesOfInterest = set()
    for read in tqdm(annotatedReads):
        containsAMRgene = False
        for g in range(len(annotatedReads[read])):
            if annotatedReads[read][g][1:] in genesOfInterest:
                containsAMRgene = True
                subsettedGenesOfInterest.add(annotatedReads[read][g][1:])
        if not containsAMRgene:
            to_delete.append(read)
    # for read in to_delete:
    #     del annotatedReads[read]
    genesOfInterest = list(subsettedGenesOfInterest)

    return annotatedReads, genesOfInterest, gene_position_dict


def get_read_start(cigar: list[tuple[int, int]]) -> int:
    """return an int of the 0 based position where the read region starts mapping to the gene"""
    # check if there are any hard clipped bases at the start of the mapping
    if cigar[0][0] == 5:
        regionStart = cigar[0][1] - 1
    else:
        regionStart = 0
    return regionStart


def get_read_end(cigar: list[tuple[int, int]], regionStart: int) -> tuple[int, int]:
    """return an int of the 0 based position where the read region stops mapping to the gene"""
    regionLength = 0
    for tuple in cigar:
        if not tuple[0] == 5:
            regionLength += tuple[1]
    regionEnd = regionStart + regionLength - 1
    return regionEnd, regionLength


def determine_gene_strand(read: pysam.libcalignedsegment.AlignedSegment) -> tuple[str, str]:
    strandlessGene = (
        read.reference_name.replace("~~~", ";")
        .replace(".aln.fas", "")
        .replace(".fasta", "")
        .replace(".fa", "")
    )
    if not read.is_forward:
        gene_name = "-" + strandlessGene
    else:
        gene_name = "+" + strandlessGene
    return gene_name, strandlessGene


def remove_poorly_mapped_genes(
    pandora_consensus,
    zero_coverage_threshold,
    genesOfInterest,
    read_path,
    cores,
    output_dir,
    consensus_file,
    minimap2_path,
):
    sys.stderr.write("\nAmira: removing genes with low coverage\n")
    # map the reads to the pandora consensus
    map_command = f"{minimap2_path} -x map-ont -a --MD --secondary=no -t {cores} "
    map_command += f"-o {os.path.join(output_dir, 'mapped_to_consensus.sam')} "
    map_command += f"{consensus_file} {read_path} && "
    map_command += (
        f"samtools sort -@ {cores} {os.path.join(output_dir, 'mapped_to_consensus.sam')} > "
    )
    map_command += f"{os.path.join(output_dir, 'mapped_to_consensus.bam')} && "
    map_command += f"samtools index {os.path.join(output_dir, 'mapped_to_consensus.bam')}"
    subprocess.run(map_command, shell=True, check=True)
    # Load the BAM file
    bam_path = os.path.join(output_dir, "mapped_to_consensus.bam")
    bam_file = pysam.AlignmentFile(bam_path, "rb")
    # Initialize coverage dictionary for each gene
    gene_coverage = {gene: [0] * bam_file.get_reference_length(gene) for gene in pandora_consensus}
    # Iterate through each read in the BAM file
    minimap2_annotatedReads = {}
    for read in tqdm(bam_file.fetch()):
        if read.is_mapped:
            if read.query_name not in minimap2_annotatedReads:
                minimap2_annotatedReads[read.query_name] = set()
            # Calculate the alignment length of the read
            alignment_length = read.query_alignment_end - read.query_alignment_start
            reference_length = bam_file.get_reference_length(read.reference_name)
            # Check if alignment is at least 80% of the reference length
            if alignment_length >= 0.9 * reference_length:
                # Add the reference name to the set for the read
                minimap2_annotatedReads[read.query_name].add(read.reference_name)

            if read.reference_name in gene_coverage:
                # Get the start and end positions of the read mapping to the reference
                start = read.reference_start
                end = read.reference_end
                # Mark the covered positions in the gene coverage list
                for pos in range(start, end):
                    gene_coverage[read.reference_name][pos] = 1
    # Close the BAM file
    bam_file.close()
    # Calculate coverage percentage for each gene
    count = 0
    for gene in gene_coverage:
        if gene in genesOfInterest:
            continue
        if (len(gene_coverage[gene]) - sum(gene_coverage[gene])) / len(
            gene_coverage[gene]
        ) > zero_coverage_threshold:
            del pandora_consensus[gene]
            count += 1
    # clean up the files
    os.remove(os.path.join(output_dir, "mapped_to_consensus.sam"))
    os.remove(os.path.join(output_dir, "mapped_to_consensus.bam"))
    os.remove(os.path.join(output_dir, "mapped_to_consensus.bam.bai"))


def convert_pandora_output(
    pandoraSam: str,
    pandora_consensus: dict[str, list[str]],
    genesOfInterest: set[str],
    geneMinCoverage: int,
    gene_length_lower_threshold: float,
    gene_length_upper_threshold: float,
    read_path: str,
    cores: int,
    output_dir: str,
    minimap2_path: str,
) -> tuple[dict[str, list[str]], list[str]]:
    # load the pseudo SAM
    pandora_sam_content = pysam.AlignmentFile(pandoraSam, "rb")
    annotatedReads: dict[str, list[str]] = {}
    gene_position_dict: dict[str, list[tuple[int, int]]] = {}
    geneCounts: dict[str, int] = {}
    # remove genes that have a high proportion of unmapped bases in the pandora consensus
    remove_poorly_mapped_genes(
        pandora_consensus,
        0.2,
        genesOfInterest,
        read_path,
        cores,
        output_dir,
        os.path.join(os.path.dirname(pandoraSam), "pandora.consensus.fq.gz"),
        minimap2_path,
    )
    # iterate through the read regions
    read_tracking = {}
    distances = {}
    for read in pandora_sam_content.fetch():
        # convert the cigarsting to a Cigar object
        cigar = read.cigartuples
        # check if the read has mapped to any regions
        if read.is_mapped:
            if read.query_name not in read_tracking:
                read_tracking[read.query_name] = {"end": 0, "index": 0}
            # get the start base that the region maps to on the read
            regionStart = get_read_start(cigar)
            # get the end base that the region maps to on the read
            regionEnd, regionLength = get_read_end(cigar, regionStart)
            # append the strand of the match to the name of the gene
            gene_name, strandlessGene = determine_gene_strand(read)
            # skip this read if there are no reads with >= 80% coverage
            # if read.query_name not in minimap2_annotatedReads:
            #    continue
            # exclude genes that do not have a pandora consensus
            if strandlessGene in genesOfInterest or (
                strandlessGene in pandora_consensus
                and gene_length_lower_threshold * len(pandora_consensus[strandlessGene]["sequence"])
                <= regionLength
                <= gene_length_upper_threshold * len(pandora_consensus[strandlessGene]["sequence"])
            ):
                read_name = read.query_name  # + "_" + str(read_tracking[read.query_name]["index"]
                if read_name not in annotatedReads:
                    annotatedReads[read_name] = []
                    gene_position_dict[read_name] = []
                # count how many times we see each gene
                if strandlessGene not in geneCounts:
                    geneCounts[strandlessGene] = 0
                geneCounts[strandlessGene] += 1
                # store the per read gene names, gene starts and gene ends
                gene_position_dict[read_name].append((regionStart, regionEnd))
                # store the per read gene names
                annotatedReads[read_name].append(gene_name)
                if read_name not in distances:
                    distances[read_name] = []
                distances[read_name].append(regionStart - read_tracking[read.query_name]["end"])
                read_tracking[read.query_name]["end"] = regionEnd
                # if strandlessGene not in fp_genes and read.query_name not in fp_reads:
                # proportion_gene_length.append(
                #     regionLength / len(pandora_consensus[strandlessGene]["sequence"])
                # )
    subsettedGenesOfInterest = set()
    for r in tqdm(annotatedReads):
        annotatedReads[r] = [
            gene for gene in annotatedReads[r] if geneCounts[gene[1:]] > geneMinCoverage - 1
        ]
        for g in range(len(annotatedReads[r])):
            if annotatedReads[r][g][1:] in genesOfInterest:
                subsettedGenesOfInterest.add(annotatedReads[r][g][1:])
    assert not len(annotatedReads) == 0
    return annotatedReads, subsettedGenesOfInterest, gene_position_dict


def process_reference_alleles(path_to_interesting_genes):
    # import the list of genes of interest
    with open(path_to_interesting_genes, "r") as i:
        reference_content = i.read().split(">")[1:]
    genesOfInterest = set()
    reference_alleles = {}
    for allele in reference_content:
        newline_split = allele.split("\n")
        assert (
            newline_split[0].count(";") == 1
        ), "Reference FASTA headers can only contain 1 semicolon"
        gene_name, allele_name = newline_split[0].split(";")
        genesOfInterest.add(gene_name)
        sequence = "".join(newline_split[1:])
        if gene_name not in reference_alleles:
            reference_alleles[gene_name] = {}
        reference_alleles[gene_name][allele_name] = sequence
    return reference_alleles, genesOfInterest
