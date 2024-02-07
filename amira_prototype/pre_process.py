import json

import pysam
from tqdm import tqdm


def process_pandora_json(
    pandoraJSON: str, genesOfInterest: list[str]
) -> tuple[dict[str, list[str]], list[str]]:
    with open(pandoraJSON) as i:
        annotatedReads = json.loads(i.read())
    to_delete = []
    subsettedGenesOfInterest = set()
    for read in tqdm(annotatedReads):
        containsAMRgene = False
        for g in range(len(annotatedReads[read])):
            split_names = annotatedReads[read][g][1:].split(".")
            if any(subgene in genesOfInterest for subgene in split_names):
                containsAMRgene = True
                subsettedGenesOfInterest.add(annotatedReads[read][g][1:])
        if not containsAMRgene:
            to_delete.append(read)
    # for read in to_delete:
    #    del annotatedReads[read]
    genesOfInterest = list(subsettedGenesOfInterest)
    return annotatedReads, genesOfInterest


def get_read_start(cigar: list[tuple[int, int]]) -> int:
    """return an int of the 0 based position where the read region starts mapping to the gene"""
    # check if there are any hard clipped bases at the start of the mapping
    if cigar[0][0] == 5:
        regionStart = cigar[0][1]
    else:
        regionStart = 0
    return regionStart


def get_read_end(cigar: list[tuple[int, int]], regionStart: int) -> tuple[int, int]:
    """return an int of the 0 based position where the read region stops mapping to the gene"""
    regionLength = 0
    for cig_tuple in cigar:  # Changed variable name from 'tuple' to 'cig_tuple'
        if cig_tuple[0] != 5:  # Using '!=' for consistency
            regionLength += cig_tuple[1]
    regionEnd = regionStart + regionLength
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


def convert_pandora_output(
    pandoraSam: str,
    pandora_consensus: dict[str, list[str]],
    genesOfInterest: set[str],
    geneMinCoverage: int,
) -> tuple[dict[str, list[str]], list[str]]:
    # load the pseudo SAM
    pandora_sam_content = pysam.AlignmentFile(pandoraSam, "rb")
    annotatedReads: dict[str, list[str]] = {}
    readLengthDict: dict[str, list[tuple[int, int]]] = {}
    geneCounts: dict[str, int] = {}
    # iterate through the read regions
    read_tracking = {}
    distances = []
    for read in pandora_sam_content.fetch():
        # convert the cigarsting to a Cigar object
        cigar = read.cigartuples
        # check if the read has mapped to any regions
        if read.is_mapped:
            if not read.query_name in read_tracking:
                read_tracking[read.query_name] = {"end": 0, "index": 0}
            # get the start base that the region maps to on the read
            regionStart = get_read_start(cigar)
            # get the end base that the region maps to on the read
            regionEnd, regionLength = get_read_end(cigar, regionStart)
            # append the strand of the match to the name of the gene
            gene_name, strandlessGene = determine_gene_strand(read)
            if regionStart - read_tracking[read.query_name]["end"] > 7000:
                read_tracking[read.query_name]["index"] += 1
            distances.append(regionStart - read_tracking[read.query_name]["end"])
            # exclude genes that do not have a pandora consensus
            if strandlessGene in pandora_consensus:
                if (
                    read.query_name + "_" + str(read_tracking[read.query_name]["index"])
                    not in annotatedReads
                ):
                    annotatedReads[
                        read.query_name + "_" + str(read_tracking[read.query_name]["index"])
                    ] = []
                    readLengthDict[
                        read.query_name + "_" + str(read_tracking[read.query_name]["index"])
                    ] = []
                # count how many times we see each gene
                if strandlessGene not in geneCounts:
                    geneCounts[strandlessGene] = 0
                geneCounts[strandlessGene] += 1
                # store the per read gene names, gene starts and gene ends
                readLengthDict[
                    read.query_name + "_" + str(read_tracking[read.query_name]["index"])
                ].append((regionStart, regionEnd))
                # store the per read gene names
                annotatedReads[
                    read.query_name + "_" + str(read_tracking[read.query_name]["index"])
                ].append(gene_name)
                read_tracking[read.query_name]["end"] = regionEnd
    to_delete = []
    subsettedGenesOfInterest = set()
    for r in tqdm(annotatedReads):
        annotatedReads[r] = [
            gene for gene in annotatedReads[r] if geneCounts[gene[1:]] > geneMinCoverage - 1
        ]
        containsAMRgene = False
        for g in range(len(annotatedReads[r])):
            split_names = annotatedReads[r][g][1:].split(".")
            if any(subgene in genesOfInterest for subgene in split_names):
                containsAMRgene = True
                for subgene in split_names:
                    if subgene in genesOfInterest:
                        annotatedReads[r][g] = annotatedReads[r][g][0] + subgene
                        break
                subsettedGenesOfInterest.add(annotatedReads[r][g][1:])
        if not containsAMRgene:
            to_delete.append(r)
    # for t in to_delete:
    #    del annotatedReads[t]
    assert not len(annotatedReads) == 0
    return annotatedReads, subsettedGenesOfInterest, distances
