import argparse
import os
import pysam
import sys
from tqdm import tqdm

from construct_graph import GeneMerGraph
from construct_unitig import UnitigTools, parse_fastq

from test_functions import TestUnitigTools

def get_options():
    """define args from the command line"""
    parser = argparse.ArgumentParser(description='Build a prototype gene de Bruijn graph.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--pandoraSam', dest='pandoraSam',
                        help='Pandora map SAM file path')
    group.add_argument('--pandoraJSON', dest='pandoraJSON',
                        help='Pandora map JSON file path')
    parser.add_argument('--pandoraConsensus', dest='pandoraConsensus',
                        help='path to Pandora consensus fastq', required=False)
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
    parser.add_argument('-p', dest='bubble_popper_threshold', type=float, default=1.5,
                        help='minimum difference in coverage threshold to collapse paths in the graph (default = 1.5)')
    parser.add_argument('-c', dest='cleaning_iterations', type=int, default=2,
                        help='number of gene-mer graph cleaning operations (default = 2)')
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
    parser.add_argument('--eval', dest='eval', action='store_true', default=False,
                        help='Amira evaluation')
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
                        pandora_consensus,
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
            # exclude genes that do not have a pandora consensus
            if strandlessGene in pandora_consensus:
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
    subsettedGenesOfInterest = set()
    for r in annotatedReads:
        #annotatedReads[r] = [gene for gene in annotatedReads[r] if geneCounts[gene[1:]] > geneMinCoverage - 1]
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
    #for t in to_delete:
     #   del annotatedReads[t]
    assert not len(annotatedReads) == 0
    return annotatedReads, readLengthDict, list(subsettedGenesOfInterest)

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

def plot_node_coverages(onemer_graph_coverages,
                    threemer_graph_coverages,
                    fivemer_graph_coverages,
                    sevenmer_graph_coverages,
                    ninemer_graph_coverages,
                    output_dir):
    import matplotlib.pyplot as plt
    for to_plot in tqdm([("one", onemer_graph_coverages),
                        ("three", threemer_graph_coverages),
                        ("five", fivemer_graph_coverages),
                        ("seven", sevenmer_graph_coverages),
                        ("nine", ninemer_graph_coverages)]):
        plt.hist(to_plot[1], bins=[i for i in range(0, 100, 1)], alpha=0.4, label=to_plot[0])
        plt.margins(x=0)
        plt.xlim(0, 100)
        plt.ylim(0,1000)
        #plt.legend(loc='upper right')
        plt.ylabel('Absolute frequency')
        plt.xlabel('Node coverage')
        plt.title('Node coverage distributions over different gene-mer lengths')
        plt.savefig(os.path.join(output_dir, f"node_coverage_distributions_{to_plot[0]}.png"))
        plt.close()

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
    if args.pandoraSam:
        # load the pandora consensus and convert to a dictionary
        pandora_consensus = parse_fastq(args.pandoraConsensus)
        sys.stderr.write("\nAmira: loading Pandora SAM\n")
        annotatedReads, readDict, genesOfInterest = convert_pandora_output(args.pandoraSam,
                                                                        pandora_consensus,
                                                                        set(genesOfInterest),
                                                                        args.gene_min_coverage)
    if args.pandoraJSON:
        import json
        with open(args.pandoraJSON) as i:
            annotatedReads = json.loads(i.read())
        to_delete = []
        subsettedGenesOfInterest = set()
        for read in tqdm(annotatedReads):
            containsAMRgene = False
            for g in range(len(annotatedReads[read])):
                split_names = annotatedReads[read][g][1:].split(".")
                if any(subgene in genesOfInterest for subgene in split_names):
                    containsAMRgene = True
                    #for subgene in split_names:
                    #    if subgene in genesOfInterest:
                    #        annotatedReads[read][g] = annotatedReads[read][g][0] + subgene
                    #        subsettedGenesOfInterest.add(annotatedReads[read][g][1:])
                    #        break
                    subsettedGenesOfInterest.add(annotatedReads[read][g][1:])
            if not containsAMRgene:
                to_delete.append(read)
        for read in to_delete:
            del annotatedReads[read]
        genesOfInterest = subsettedGenesOfInterest
    print(genesOfInterest)
    if args.debug:
        raw_graph = GeneMerGraph(annotatedReads,
                                args.geneMer_size,
                                args.eval)
        # color nodes in the graph
        for node in raw_graph.all_nodes():
            node.color_node(genesOfInterest)
        import json
        with open(os.path.join(args.output_dir, "genesAnnotatedOnReads.json"), "w") as o:
            o.write(json.dumps(annotatedReads))
        sys.stderr.write("\nAmira: writing pre-correction gene-mer graph\n")
        raw_graph.generate_gml(os.path.join(args.output_dir, "pre_correction_gene_mer_graph"),
                                    args.geneMer_size,
                                    1,
                                    1)
        # sys.stderr.write("\nAmira: plotting node coverage distributions at different values of k\n")
        # onemer_graph_coverages = GeneMerGraph(annotatedReads, 1).get_all_node_coverages()
        # threemer_graph_coverages = GeneMerGraph(annotatedReads, 3).get_all_node_coverages()
        # fivemer_graph_coverages = GeneMerGraph(annotatedReads, 5).get_all_node_coverages()
        # sevenmer_graph_coverages = GeneMerGraph(annotatedReads, 7).get_all_node_coverages()
        # ninemer_graph_coverages = GeneMerGraph(annotatedReads, 9).get_all_node_coverages()
        # plot_node_coverages(onemer_graph_coverages,
        #                     threemer_graph_coverages,
        #                     fivemer_graph_coverages,
        #                     sevenmer_graph_coverages,
        #                     ninemer_graph_coverages,
        #                     args.output_dir)
        #sys.exit(0)
    # clean the graph iteratively
    i = 1
    while i <= args.cleaning_iterations:
        sys.stderr.write(f"\nAmira: Correcting annotation errors by trimming hairs and popping bubbles {i}/{args.cleaning_iterations}\n")
        graphToCorrect = GeneMerGraph(annotatedReads,
                                    args.geneMer_size,
                                    args.eval)
        annotatedReads = graphToCorrect.correct_errors(10,
                                                    args.bubble_popper_threshold)
        i += 1
    sys.stderr.write("\nAmira: building corrected gene-mer graph\n")
    graph = GeneMerGraph(annotatedReads,
                        args.geneMer_size,
                        args.eval)
    graph.filter_graph(args.node_min_coverage,
                    args.edge_min_coverage)
    if args.debug:
        # color nodes in the graph
        for node in graph.all_nodes():
            node.color_node(genesOfInterest)
    # #######################
    readNodes = graph.get_readNodes()
    annotatedReads = {}
    for readId in tqdm(readNodes):
        for i in range(len(readNodes[readId]) - 1):
            sourceNode = graph.get_node_by_hash(readNodes[readId][i])
            targetNode = graph.get_node_by_hash(readNodes[readId][i+1])
            if not graph.check_if_nodes_are_adjacent(sourceNode, targetNode):
                graph.add_edge(sourceNode.get_geneMer(), targetNode.get_geneMer())
        annotatedReads[readId] = graph.follow_path_to_get_annotations(readNodes[readId])
    graph = GeneMerGraph(annotatedReads,
                        args.geneMer_size,
                        args.eval)
    # ####
    sys.stderr.write("\nAmira: writing gene-mer graph\n")
    graph.generate_gml(os.path.join(args.output_dir, "gene_mer_graph"),
                    args.geneMer_size,
                    args.node_min_coverage,
                    args.edge_min_coverage)
    # initialise the UnitigBuilder class
    unitigTools = TestUnitigTools(graph,
                            genesOfInterest,
                            args.readfile,
                            args.output_dir)
    # subset the reads corresponding to each paralog
    readFiles = unitigTools.assign_reads_to_amr_genes()
    sys.exit(0)
    sys.stderr.write("\nAmira: separating paralog reads\n")
    # initialise the UnitigBuilder class
    unitigTools = UnitigTools(graph,
                            genesOfInterest,
                            args.readfile,
                            args.output_dir)
    # subset the reads corresponding to each paralog
    readFiles = unitigTools.separate_paralogs()
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
    # sys.stderr.write("\nAmira: generating unitig plots\n")
    # unitigTools.visualise_unitigs(readDict,
    #                             unitig_mapping,
    #                             args.output_dir)
    sys.exit(0)

if __name__ == "__main__":
    main()