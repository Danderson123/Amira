import json

import pysam
from tqdm import tqdm


def clean_gene(g):
    chars_to_remove = set(["|", "(", ")", "-", "*", "+", "#", ":", "=", "/", ",", "'"])
    cleaned_gene = "".join(char for char in g if char not in chars_to_remove)
    return cleaned_gene


def process_pandora_json(
    pandoraJSON: str, positionJSON: str, genesOfInterest: list[str]
) -> tuple[dict[str, list[str]], list[str]]:
    with open(pandoraJSON) as i:
        annotatedReads = json.loads(i.read())
    with open(positionJSON) as i:
        genePositions = json.loads(i.read())
    ###################
    # # just needed for running on simulated data
    # newAnnotatedReads = {}
    # for r in annotatedReads:
    #     newAnnotatedReads[r] = ["_".join(g.split("_")[:-1]) for g in annotatedReads[r]]
    # annotatedReads = newAnnotatedReads
    ###################
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
    #    del annotatedReads[read]
    genesOfInterest = list(subsettedGenesOfInterest)
    return annotatedReads, genesOfInterest, genePositions


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


def remove_poorly_mapped_genes(pandora_consensus, zero_coverage_threshold, genesOfInterest):
    # iterate through the read regions
    zero_coverage_base_proportion = []
    genes_to_delete = []
    for gene in pandora_consensus:
        zero_coverage_proportion = pandora_consensus[gene]["quality"].count("!") / len(
            pandora_consensus[gene]["quality"]
        )
        zero_coverage_base_proportion.append(zero_coverage_proportion)
        if zero_coverage_proportion > zero_coverage_threshold and gene not in genesOfInterest:
            genes_to_delete.append(gene)
    for gene in genes_to_delete:
        del pandora_consensus[gene]


def convert_pandora_output(
    pandoraSam: str,
    pandora_consensus: dict[str, list[str]],
    genesOfInterest: set[str],
    geneMinCoverage: int,
    gene_length_lower_threshold,
    gene_length_upper_threshold,
) -> tuple[dict[str, list[str]], list[str]]:
    # load the pseudo SAM
    pandora_sam_content = pysam.AlignmentFile(pandoraSam, "rb")
    annotatedReads: dict[str, list[str]] = {}
    gene_position_dict: dict[str, list[tuple[int, int]]] = {}
    geneCounts: dict[str, int] = {}
    # remove genes that have a high proportion of unmapped bases in the pandora consensus
    remove_poorly_mapped_genes(pandora_consensus, 0.2, genesOfInterest)
    # iterate through the read regions
    read_tracking = {}
    distances = {}
    proportion_gene_length = []
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
            if strandlessGene in pandora_consensus and (
                gene_length_lower_threshold * len(pandora_consensus[strandlessGene]["sequence"])
                <= regionLength
                <= gene_length_upper_threshold * len(pandora_consensus[strandlessGene]["sequence"])
                or strandlessGene in genesOfInterest
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
                if not read_name in distances:
                    distances[read_name] = []
                distances[read_name].append(regionStart - read_tracking[read.query_name]["end"])
                read_tracking[read.query_name]["end"] = regionEnd
                # if strandlessGene not in fp_genes and read.query_name not in fp_reads:
                proportion_gene_length.append(
                    regionLength / len(pandora_consensus[strandlessGene]["sequence"])
                )
    to_delete = []
    subsettedGenesOfInterest = set()
    for r in tqdm(annotatedReads):
        annotatedReads[r] = [
            gene for gene in annotatedReads[r] if geneCounts[gene[1:]] > geneMinCoverage - 1
        ]
        containsAMRgene = False
        for g in range(len(annotatedReads[r])):
            if annotatedReads[r][g][1:] in genesOfInterest:
                containsAMRgene = True
                subsettedGenesOfInterest.add(annotatedReads[r][g][1:])
            # split_names = annotatedReads[r][g][1:].split(".")
            # if any(subgene in genesOfInterest for subgene in split_names):
            #     containsAMRgene = True
            #     subsettedGenesOfInterest.add(annotatedReads[r][g][1:])
        if not containsAMRgene:
            to_delete.append(r)
    # for t in to_delete:
    #    del annotatedReads[t]
    assert not len(annotatedReads) == 0
    return annotatedReads, subsettedGenesOfInterest, gene_position_dict
