import glob
import json
import os
import statistics
import subprocess
import sys

import pysam
from tqdm import tqdm


def run_pandora_map(pandora_path, panRG_path, readfile, outdir, cores, seed):
    command = f"{pandora_path} map -t {cores} --min-gene-coverage-proportion 0.5 --max-covg 10000 "
    command += (
        f"-o {os.path.join(outdir, 'pandora_output')} {panRG_path} {readfile} --rng-seed {seed} "
    )
    command += "--min-abs-gene-coverage 1"
    # check that the panRG file exists
    if not os.path.exists(panRG_path):
        sys.stderr.write("\nAmira: panRG file does not exist.\n")
        sys.exit(1)
    if ".panidx.zip" not in panRG_path:
        sys.stderr.write("\nAmira: panRG file does not end in .panidx.zip.\n")
        sys.exit(1)
    # run pandora
    subprocess.run(command, shell=True, check=True)
    pandoraSam = glob.glob(os.path.join(outdir, "pandora_output", "*.filtered.sam"))[0]
    pandoraConsensus = os.path.join(outdir, "pandora_output", "pandora.consensus.fq.gz")
    return pandoraSam, pandoraConsensus


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


def load_species_specific_files(species, AMR_gene_reference_FASTA, sequence_names, core_genes):
    from pathlib import Path

    if AMR_gene_reference_FASTA is None or sequence_names is None or core_genes is None:
        species_dir = Path(__file__).resolve().parent / "assets"
        if not os.path.exists(f"{species_dir}/{species}"):
            sys.stderr.write(f"\nAmira: {species} is not a supported species name.\n")
            sys.exit(1)
        else:
            if AMR_gene_reference_FASTA is None:
                AMR_gene_reference_FASTA = os.path.join(
                    f"{species_dir}/{species}", "AMR_alleles_unified.fa"
                )
            if sequence_names is None:
                sequence_names = os.path.join(f"{species_dir}/{species}", "AMR_calls.json")
            if core_genes is None:
                core_genes = os.path.join(f"{species_dir}/{species}", "core_genes.txt")
        return AMR_gene_reference_FASTA, sequence_names, core_genes
    else:
        return AMR_gene_reference_FASTA, sequence_names, core_genes


def remove_poorly_mapped_genes(
    pandora_consensus,
    zero_coverage_threshold,
    genesOfInterest,
    read_path,
    cores,
    output_dir,
    consensus_file,
    minimap2_path,
    samtools_path,
):
    sys.stderr.write("\nAmira: removing genes with low coverage\n")
    # map the reads to the pandora consensus
    map_command = f"{minimap2_path} -x map-ont -a --MD --secondary=no -t {cores} "
    map_command += f"-o {os.path.join(output_dir, 'mapped_to_consensus.sam')} "
    map_command += f"{consensus_file} {read_path} && "
    map_command += (
        f"{samtools_path} sort -@ {cores} {os.path.join(output_dir, 'mapped_to_consensus.sam')} > "
    )
    map_command += f"{os.path.join(output_dir, 'mapped_to_consensus.bam')} && "
    map_command += f"{samtools_path} index {os.path.join(output_dir, 'mapped_to_consensus.bam')}"
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


def convert_pandora_output(
    pandoraSam: str,
    pandora_consensus: dict[str, list[str]],
    genesOfInterest: set[str],
    relativeMinGeneThreshold: int,
    gene_length_lower_threshold: float,
    gene_length_upper_threshold: float,
    read_path: str,
    cores: int,
    output_dir: str,
    minimap2_path: str,
    samtools_path: str,
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
        samtools_path,
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
            # exclude genes that do not have a pandora consensus
            if strandlessGene in genesOfInterest or (
                strandlessGene in pandora_consensus
                and gene_length_lower_threshold * len(pandora_consensus[strandlessGene]["sequence"])
                <= regionLength
                <= gene_length_upper_threshold * len(pandora_consensus[strandlessGene]["sequence"])
            ):
                read_name = read.query_name
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
    # get the min relative gene frequency
    geneMinCoverage = statistics.mean(geneCounts.values()) * relativeMinGeneThreshold
    # filter genes that do not meet the minimum coverage threshold
    subsettedGenesOfInterest = set()
    for r in tqdm(annotatedReads):
        new_calls = []
        new_positions = []
        for i in range(len(annotatedReads[r])):
            gene = annotatedReads[r][i]
            if geneCounts[gene[1:]] >= geneMinCoverage:
                new_calls.append(gene)
                new_positions.append(gene_position_dict[r][i])
                if gene[1:] in genesOfInterest:
                    subsettedGenesOfInterest.add(gene[1:])
        annotatedReads[r] = new_calls
        gene_position_dict[r] = new_positions
    assert not len(annotatedReads) == 0
    return annotatedReads, subsettedGenesOfInterest, gene_position_dict


def process_reference_alleles(path_to_interesting_genes, promoters):
    # import the list of genes of interest
    with open(path_to_interesting_genes, "r") as i:
        reference_content = i.read().split(">")[1:]
    reference_alleles = {}
    genesOfInterest = set()
    promoter_alleles = []
    for allele in reference_content:
        newline_split = allele.split("\n")
        assert (
            newline_split[0].count(";") == 1
        ), "Reference FASTA headers can only contain 1 semicolon"
        gene_name, allele_name = newline_split[0].split(";")
        sequence = "".join(newline_split[1:])
        if "promoter" in gene_name:
            promoter_alleles.append((gene_name.replace("_promoter", ""), allele_name, sequence))
            continue
        genesOfInterest.add(gene_name)
        if gene_name not in reference_alleles:
            reference_alleles[gene_name] = {}
        reference_alleles[gene_name][allele_name] = sequence
    # if promoters are specified then add them
    if promoters is True:
        promoters_to_add = {}
        for gene_name in reference_alleles:
            for promoter_allele in promoter_alleles:
                if promoter_allele[0] in gene_name:
                    promoter_name = gene_name + "_promoter"
                    if promoter_name not in promoters_to_add:
                        promoters_to_add[promoter_name] = {}
                    promoters_to_add[promoter_name][promoter_allele[1]] = promoter_allele[2]
        reference_alleles.update(promoters_to_add)
    return reference_alleles, genesOfInterest


def samtools_get_mean_depth(bam_file, core_genes, samtools_path):
    # Run samtools coverage and capture output
    result = subprocess.run(
        [samtools_path, "coverage", bam_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True,
    )
    # read the output
    mean_depth_per_contig = {}
    for line in result.stdout.strip().split("\n"):
        if line.startswith("#"):
            continue
        rname, startpos, endpos, numreads, covbases, coverage, meandepth, meanbaseq, meanmapqs = (
            line.split("\t")
        )
        if rname in core_genes:
            mean_depth_per_contig[rname] = float(meandepth)
    return mean_depth_per_contig


def get_core_gene_mean_depth(bam_file, core_gene_file, samtools_path):
    # load the core genes
    with open(core_gene_file) as i:
        core_genes = set(i.read().split("\n"))
    # get the mean depth across each core gene
    mean_depth_per_core_gene = samtools_get_mean_depth(bam_file, core_genes, samtools_path)
    os.remove(bam_file)
    os.remove(bam_file + ".bai")
    return statistics.mean(list(mean_depth_per_core_gene.values()))
