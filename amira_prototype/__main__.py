import argparse
import os
import pysam
import sys

from construct_graph import GeneMerGraph
from construct_unitig import UnitigTools
from preprocess import process_sam

def get_options():
    """define args from the command line"""
    parser = argparse.ArgumentParser(description='Build a prototype gene de Bruijn graph.')
    parser.add_argument('--pandoraSam', dest='pandoraSam',
                        help='Pandora map SAM file path', required=True)
    parser.add_argument('--readfile', dest='readfile',
                        help='path of gzipped long read fastq', required=True)
    parser.add_argument('--output', dest='output_dir', type=str, default="gene_de_Bruijn_graph",
                        help='directory for Amira outputs')
    parser.add_argument('-k', dest='geneMer_size', type=int, default=3,
                        help='kmer length for the gene de Bruijn graph')
    parser.add_argument('-n', dest='node_min_coverage', type=int, default=1,
                        help='minimum threshold for gene-mer coverage')
    parser.add_argument('-e', dest='edge_min_coverage', type=int, default=1,
                        help='minimum threshold for edge coverage between gene-mers')
    parser.add_argument('-g', dest='gene_min_coverage', type=int, default=1,
                        help='minimum threshold for gene filtering')
    parser.add_argument('--gene-path', dest='path_to_interesting_genes',
                        help='path to a newline delimited file of genes of interest', required=True)
    parser.add_argument('--flye-path', dest='flye_path',
                        help='path to Flye binary', default=None, required=False)
    parser.add_argument('--raven-path', dest='raven_path',
                        help='path to Raven binary', default=None, required=False)
    parser.add_argument('--use-consensus', dest='use_pandora_consensus',
                        help='polish the pandora consensus of each gene to recover AMR alleles',
                        action='store_true' , default=None, required=False)
    parser.add_argument("--pandora-consensus", dest="pandoraConsensus",
                        help='path to the pandora consensus fastq for this sample', default=None, required=False)
    parser.add_argument('--racon-path', dest='racon_path',
                        help='path to Racon binary', default=None, required=False)
    parser.add_argument('--threads', dest='threads', type=int, default=1,
                        help='number of threads to use')
    parser.add_argument('--debug', dest='debug', action='store_true', default=False,
                        help='Amira debugging')
    args = parser.parse_args()
    return args

def get_read_start(cigar):
    """ return an int of the 0 based position where the read region starts mapping to the gene """
    # check if there are any hard clipped bases at the start of the mapping
    if cigar[0][0] == 5:
        regionStart = cigar[0][1]
    else:
        regionStart = 0
    return regionStart

def get_read_end(cigar,
                regionStart):
    """ return an int of the 0 based position where the read region stops mapping to the gene """
    regionLength = 0
    for tuple in cigar:
        if not tuple[0] == 5:
            regionLength += tuple[1]
    regionEnd = regionStart + regionLength
    return regionEnd, regionLength

def determine_gene_strand(read):
    strandlessGene = read.reference_name.replace("~~~", ";").replace(".aln.fas", "").replace(".fasta", "").replace(".fa", "")
    if not read.is_forward:
        gene_name = "-" + strandlessGene
    else:
        gene_name = "+" + strandlessGene
    return gene_name, strandlessGene

def convert_pandora_output(pandoraSam,
                        genesOfInterest,
                        geneMinCoverage):
    # load the pseudo SAM
    pandora_sam_content = pysam.AlignmentFile(pandoraSam, "rb")
    annotatedReads = {}
    readLengthDict = {}
    geneCounts = {}
    # iterate through the read regions
    for read in pandora_sam_content.fetch():
        # convert the cigarsting to a Cigar object
        cigar = read.cigartuples
        # check if the read has mapped to any regions
        if read.is_mapped:
            # get the start base that the region maps to on the read
            regionStart = get_read_start(cigar)
            # get the end base that the region maps to on the read
            regionEnd, regionLength = get_read_end(cigar,
                                                regionStart)
            # append the strand of the match to the name of the gene
            gene_name, strandlessGene = determine_gene_strand(read)
            if not read.query_name in annotatedReads:
                annotatedReads[read.query_name] = []
                readLengthDict[read.query_name] = []
            # count how many times we see each gene
            if not strandlessGene in geneCounts:
                geneCounts[strandlessGene] = 0
            geneCounts[strandlessGene] += 1
            # store the per read gene names, gene starts and gene ends
            readLengthDict[read.query_name].append((regionStart, regionEnd))
            # store the per read gene names
            annotatedReads[read.query_name].append(gene_name)
    to_delete = []
    for r in annotatedReads:
        annotatedReads[r] = [gene for gene in annotatedReads[r] if geneCounts[gene[1:]] > geneMinCoverage - 1]
        if not any(g[1:] in genesOfInterest for g in annotatedReads[r]):
            to_delete.append(r)
    for t in to_delete:
        del annotatedReads[t]
    assert not len(annotatedReads) == 0
    return annotatedReads, readLengthDict

def plot_gene_counts(annotatedReads,
                    outputDir):
    import matplotlib.pyplot as plt
    import numpy as np
    # make the output dir if it doesn't exist
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)
    # keep track of the number of genes per read
    geneCounts = []
    geneCountDict = {}
    for r in annotatedReads:
        numberOfGenes = len(annotatedReads[r])
        geneCounts.append(numberOfGenes)
        geneCountDict[r] = numberOfGenes
    # plot a bar chart of the number of genes per read
    q25, q75 = np.percentile(geneCounts, [25, 75])
    bin_width = 2 * (q75 - q25) * len(geneCounts) ** (-1/3)
    numBins = round((max(geneCounts) - min(geneCounts)) / bin_width)
    fig, ax = plt.subplots()
    ax.hist(geneCounts, bins=numBins, density = False)
    plt.margins(x=0)
    ax.set_ylabel('Absolute frequency')
    ax.set_xlabel('Number of genes')
    fig.savefig(os.path.join(outputDir, "read_gene_count_distribution.png"))
    plt.close(fig)

def plot_read_lengths(readFile,
                    annotatedReads,
                    outputDir):
    from construct_unitig import parse_fastq
    import matplotlib.pyplot as plt
    import numpy as np
    readDict = parse_fastq(readFile)
    readLengths = []
    for r in readDict:
        if r in annotatedReads:
            sequence = readDict[r]["sequence"]
            readLengths.append(len(sequence))
    # plot a histogram of read lengths in base pairs
    q25, q75 = np.percentile(readLengths, [25, 75])
    bin_width = 2 * (q75 - q25) * len(readLengths) ** (-1/3)
    numBins = round((max(readLengths) - min(readLengths)) / bin_width)
    fig, ax = plt.subplots()
    ax.hist(readLengths, bins=numBins, density = False)
    plt.margins(x=0)
    ax.set_ylabel('Absolute frequency')
    ax.set_xlabel('Length of read (bp)')
    fig.savefig(os.path.join(outputDir, "read_length_distribution.png"))
    plt.close(fig)

def main():
    # get command line options
    args = get_options()
    # make the output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # import the list of genes of interest
    with open(args.path_to_interesting_genes, "r") as i:
        genesOfInterest = i.read().splitlines()
    # convert the Pandora SAM file to a dictionary
    sys.stderr.write("\nAmira: loading Pandora SAM\n")
    annotatedReads, readDict = convert_pandora_output(args.pandoraSam,
                                                    set(genesOfInterest),
                                                    args.gene_min_coverage)
    # build the graph
    sys.stderr.write("\nAmira: building pre-correction gene-mer graph\n")
    graphToCorrect = GeneMerGraph(annotatedReads,
                        args.geneMer_size)
    if args.debug:
        sys.stderr.write("\nAmira: writing pre-correction gene-mer graph\n")
        graphToCorrect.generate_gml(os.path.join(args.output_dir, "pre_correction_gene_mer_graph"),
                                    args.geneMer_size,
                                    1,
                                    1)
    sys.stderr.write("\nCorrecting annotation errors by removing appendages and popping bubbles\n")
    annotatedReads = graphToCorrect.correct_errors(5)
    sys.stderr.write("\nAmira: building corrected gene-mer graph\n")
    graph = GeneMerGraph(annotatedReads,
                        args.geneMer_size)
    sys.stderr.write("\nAmira: filtering gene-mer graph\n")
    graph.filter_graph(args.node_min_coverage,
                    args.edge_min_coverage)
    if args.debug:
        # color nodes in the graph
        for node in graph.all_nodes():
            node.color_node(genesOfInterest)
        # plot_gene_counts(annotatedReads,
        #                 args.output_dir)
        # plot_read_lengths(args.readfile,
        #                 annotatedReads,
        #                 args.output_dir)
    sys.stderr.write("\nAmira: writing gene-mer graph\n")
    graph.generate_gml(os.path.join(args.output_dir, "gene_mer_graph"),
                    args.geneMer_size,
                    args.node_min_coverage,
                    args.edge_min_coverage)
    # initialise the UnitigBuilder class
    sys.stderr.write("\nAmira: getting viable paths\n")
    unitigTools = UnitigTools(graph,
                            genesOfInterest)
    # generate a visualisation of the unitigs
    sys.stderr.write("\nAmira: separating paralog reads\n")
    readFiles, unitig_mapping = unitigTools.separate_reads(os.path.join(args.output_dir, "unitigs"),
                                                        args.readfile)
    # run flye on the subsetted reads
    if args.flye_path:
        sys.stderr.write("\nAmira: assembling paralog reads with Flye\n")
        unitigTools.multithread_flye(readFiles,
                                    args.flye_path,
                                    args.threads)
    # run raven on subsetted reads and pandora consensus
    if args.raven_path:
        sys.stderr.write("\nAmira: assembling paralog reads with Raven\n")
        unitigTools.multithread_raven(readFiles,
                                    args.raven_path,
                                    args.threads)
    # use each gene's pandora consensus as the basis for polishing
    if args.use_pandora_consensus:
        sys.stderr.write("\nAmira: polishing Pandora AMR gene consensus sequences\n")
        unitigTools.polish_pandora_consensus(readFiles,
                                            args.racon_path,
                                            args.pandoraConsensus,
                                            genesOfInterest,
                                            args.threads)
    # make plots to visualise unitigs
    sys.stderr.write("\nAmira: generating unitig plots\n")
    unitigTools.visualise_unitigs(readDict,
                                unitig_mapping,
                                args.output_dir)
    sys.exit(0)

if __name__ == "__main__":
    main()