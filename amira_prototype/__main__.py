import argparse
from cigar import Cigar
import os
import pysam
import sys
from tqdm import tqdm

from construct_graph import GeneMerGraph
from construct_unitig import UnitigTools, Unitigs
from construct_read import Read

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
    parser.add_argument('--gene-path', dest='path_to_interesting_genes',
                        help='path to a newline delimited file of genes of interest', required=True)
    parser.add_argument('--flye-path', dest='flye_path',
                        help='path to Flye binary', default=None, required=False)
    parser.add_argument('--raven-path', dest='raven_path',
                        help='path to Raven binary', default=None, required=False)
    parser.add_argument('--use-consensus', dest='use_pandora_consensus',
                        help='polish the pandora consensus of each gene to recover AMR alleles',
                        action='store_true', default=False, required=False)
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

def get_read_start(cigarString):
    """ return an int of the 0 based position where the read region starts mapping to the gene """
    cigarSplit = list(cigarString.items())
    # check if there are any hard clipped bases at the start of the mapping
    if "H" in cigarSplit[0]:
        regionStart = cigarSplit[0][0]
    else:
        regionStart = 0
    return regionStart

def get_read_end(cigarString,
                regionStart):
    """ return an int of the 0 based position where the read region stops mapping to the gene """
    regionLength = len(cigarString)
    regionEnd = regionStart + regionLength
    return regionEnd, regionLength

def determine_gene_strand(read):
    if not read.is_forward:
        gene_name = "-" + read.reference_name.replace("~~~", ";")
    else:
        gene_name = "+" + read.reference_name.replace("~~~", ";")
    return gene_name

def convert_pandora_output(pandoraSam):
    # load the pseudo SAM
    pandora_sam_content = pysam.AlignmentFile(pandoraSam, "rb")
    annotatedReads = {}
    readDict = {}
    # iterate through the read regions
    for read in pandora_sam_content.fetch():
        # convert the cigarsting to a Cigar object
        cigarString = Cigar(read.cigarstring)
        # check if the read has mapped to any regions
        if read.is_mapped:
            # get the start base that the region maps to on the read
            regionStart = get_read_start(cigarString)
            # get the end base that the region maps to on the read
            regionEnd, regionLength = get_read_end(cigarString,
                                                regionStart)
            # append the strand of the match to the name of the gene
            gene_name = determine_gene_strand(read)
            if not read.query_name in annotatedReads:
                annotatedReads[read.query_name] = []
            if not gene_name[1:] in readDict:
                readDict[gene_name[1:]] = []
            # store the per read gene names, gene starts and gene ends
            readDict[gene_name[1:]].append(regionLength)
            # store the per read gene names
            annotatedReads[read.query_name].append(gene_name)
    return annotatedReads, readDict

def plot_gene_counts(annotatedReads,
                    outputDir):
    import matplotlib.pyplot as plt
    # make the output dir if it doesn't exist
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)
    # keep track of the number of genes per read
    geneCounts = {}
    for readId in annotatedReads:
        numberOfGenes = len(annotatedReads[readId])
        if not numberOfGenes in geneCounts:
            geneCounts[numberOfGenes] = 0
        geneCounts[numberOfGenes] += 1
    # plot a bar chart of the number of genes per read
    plt.hist(list(geneCounts.keys()), weights=list(geneCounts.values()), bins=50, density = False)
    plt.ylabel('Frequency')
    plt.xlabel('Number of genes')
    # save the bar chart
    plt.savefig(os.path.join(outputDir, "read_gene_count_distribution.pdf"))

def enumerate_paths(graph,
                    outputDir):
    uniqueReadPaths = {}
    # iterate through nodes in the graph
    nonUniqueReads = {}
    nodePaths = {}
    for node in graph.all_nodes():
        # if this node contains an AMR gene and is branching
        if node.get_color() == 2:
            # get forward nodes
            possible_forward_nodes = [graph.get_edge_by_hash(h) for h in node.get_forward_edge_hashes()]
            # get backward nodes
            possible_backward_nodes = [graph.get_edge_by_hash(h) for h in node.get_backward_edge_hashes()]
            # iterate through forward nodes
            for fw in possible_forward_nodes:
                if not fw.get_targetNode() == node:
                    fw_targetNode = fw.get_targetNode()
                else:
                    fw_targetNode = fw.get_sourceNode()
                if not graph.get_degree(fw_targetNode) > 2:
                    # iterate through backward nodes
                    for bw in possible_backward_nodes:
                        if not bw.get_targetNode() == node:
                            bw_targetNode = bw.get_targetNode()
                        else:
                            bw_targetNode = bw.get_sourceNode()
                        if not graph.get_degree(bw_targetNode) > 2:
                            # hash the path
                            pathHash = hash(tuple(list(sorted([fw_targetNode.__hash__(), bw_targetNode.__hash__()]))))
                            # store the pathHash as being associated with these nodes
                            for h in [fw_targetNode.__hash__(), bw_targetNode.__hash__(), node.__hash__()]:
                                if not h in nodePaths:
                                    nodePaths[h] = set()
                                nodePaths[h].add(str(pathHash))
                            # store the reads in the path
                            fw_reads = [r for r in fw_targetNode.get_reads()]
                            bw_reads = [r for r in bw_targetNode.get_reads()]
                            for read in list(set(fw_reads + bw_reads)):
                                if not read in uniqueReadPaths:
                                    uniqueReadPaths[read] = pathHash
                                else:
                                    if not uniqueReadPaths[read] == pathHash:
                                        if not read in nonUniqueReads:
                                            nonUniqueReads[read] = set()
                                        nonUniqueReads[read].add(pathHash)
    for d in nonUniqueReads:
        del uniqueReadPaths[d]
    # add non unique reads
    uniqueCounts = {}
    for read in uniqueReadPaths:
        if not uniqueReadPaths[read] in uniqueCounts:
            uniqueCounts[uniqueReadPaths[read]] = {"unique": [], "non-unique": []}
        uniqueCounts[uniqueReadPaths[read]]["unique"].append(read)
    for read in nonUniqueReads:
        for path in list(nonUniqueReads[read]):
            if not path in uniqueCounts:
                uniqueCounts[path] = {"unique": [], "non-unique": []}
            uniqueCounts[path]["non-unique"].append(read)
    # write to a tab-delimited file
    if os.path.exists(os.path.join(outputDir, "path_reads.tab")):
        os.remove(os.path.join(outputDir, "path_reads.tab"))
    for pathHash in uniqueCounts:
        uniqueProportion = str(len(uniqueCounts[pathHash]["unique"])) + "/" + str(len(uniqueCounts[pathHash]["unique"]) + len(uniqueCounts[pathHash]["non-unique"]))
        with open(os.path.join(outputDir, "path_reads.tab"), "a+") as o:
            o.write(str(pathHash) + "\t" + uniqueProportion + "\t" + ", ".join(uniqueCounts[pathHash]["unique"]) + "\n")
    return nodePaths

def get_path_unique_reads(unitig,
                    outputDir):
    # get all nodes containing AMR genes
    nodeOfInterest = {}
    seen_paths = []
    graph = unitig.get_graph()
    for node in tqdm(graph.all_nodes()):
        if node.get_color() == 2:
            # get forward nodes
            possible_forward_nodes = [graph.get_edge_by_hash(h) for h in node.get_forward_edge_hashes()]
            # get backward nodes
            possible_backward_nodes = [graph.get_edge_by_hash(h) for h in node.get_backward_edge_hashes()]
            # iterate through forward nodes
            all_fw_paths = []
            for fw in possible_forward_nodes:
                if not fw.get_targetNode() == node:
                    fw_targetNode = fw.get_targetNode()
                else:
                    fw_targetNode = fw.get_sourceNode()
                if not graph.get_degree(fw_targetNode) > 2:
                    fw_path = unitig.get_node_forward_path_from_node(fw_targetNode)
                else:
                    fw_path = None
            all_fw_paths.append(fw_path)
            # iterate through backward nodes
            all_bw_paths = []
            for bw in possible_backward_nodes:
                if not bw.get_targetNode() == node:
                    bw_targetNode = bw.get_targetNode()
                else:
                    bw_targetNode = bw.get_sourceNode()
                if not graph.get_degree(bw_targetNode) > 2:
                    bw_path = unitig.get_node_backward_path_from_node(bw_targetNode)
                else:
                    bw_path = None
                all_bw_paths.append(bw_path)
            # enumerate all possible paths
            for f in all_fw_paths:
                if f:
                    for b in all_bw_paths:
                        if b:
                            # hash the path
                            path = b + [node] + f
                            pathHash = hash(tuple([n.__hash__() for n in path]))
                            nodeOfInterest[pathHash] = path
    for node in tqdm(graph.all_nodes()):
        if node.get_color() == 1:
            fw_path = unitig.get_node_forward_path_from_node(node)
            bw_path = unitig.get_node_backward_path_from_node(node)
            pathNodes = bw_path + [node] + fw_path
            pathHash = hash(tuple([n.__hash__() for n in pathNodes if n]))
            if not any(all(pathNodes in s) for s in seen_paths):
                nodeOfInterest[pathHash] = pathNodes
    # get all reads unique to each path
    readPaths = {}
    pathReads = {}
    readablePaths = {}
    for p in tqdm(nodeOfInterest):
        all_reads = set()
        for node in nodeOfInterest[p]:
            for r in node.get_reads():
                all_reads.add(r)
        for r in list(all_reads):
            if not r in readPaths:
                readPaths[r] = p
            else:
                if not readPaths[r] == p:
                    readPaths[r] = None
                    pathReads[p] = []
        readablePaths[p] = [graph.get_gene_mer_label(n) for n in nodeOfInterest[p]]
    for r in readPaths:
        if readPaths[r]:
            if not readPaths[r] in pathReads:
                pathReads[readPaths[r]] = []
            pathReads[readPaths[r]].append(r)
    if os.path.exists(os.path.join(outputDir, "unique_reads.tab")):
        os.remove(os.path.join(outputDir, "unique_reads.tab"))
    for path in pathReads:
        with open(os.path.join(outputDir, "unique_reads.tab"), "a+") as o:
            o.write(str(path) + "\t" + ", ".join(pathReads[path]) + "\n")
    import json
    with open(os.path.join(outputDir, "readable_paths.json"), "w") as o:
        o.write(json.dumps(readablePaths))

def get_unique_reads_for_branched(unitig,
                                outputDir):
    graph = unitig.get_graph()
    readNode = {}
    for node in tqdm(graph.all_nodes()):
        if node.get_color() == 2:
            for read in node.get_reads():
                if not read in readNode:
                    readNode[read] = graph.get_gene_mer_label(node)
                else:
                    if not readNode[read] == graph.get_gene_mer_label(node):
                        readNode[read] = None
    nodeReads = {}

def main():
    # get command line options
    args = get_options()
    # make the output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # convert the Pandora SAM file to a dictionary
    sys.stderr.write("\nAmira: loading Pandora SAM\n")
    annotatedReads, readDict = convert_pandora_output(args.pandoraSam)
    # import the list of genes of interest
    with open(args.path_to_interesting_genes, "r") as i:
        genesOfInterest = i.read().splitlines()
    # build the graph
    sys.stderr.write("\nAmira: building gene-mer graph\n")
    graph = GeneMerGraph(annotatedReads,
                        args.geneMer_size)
    sys.stderr.write("\nAmira: filtering gene-mer graph\n")
    graph.filter_graph(args.node_min_coverage,
                    args.edge_min_coverage)
    # if debugging then color nodes with AMR genes
    if args.debug:
        # plot a histogram of gene counts per read
        plot_gene_counts(annotatedReads,
                        os.path.join(args.output_dir,
                                    "debug"))
        # add a color attribute to each node
        for node in tqdm(graph.all_nodes()):
            node.color_node(genesOfInterest)
        # get the nodes containing AMR genes
        unitig = Unitigs(graph,
                        genesOfInterest)
        # get reads unique to each node with a degree > 2 and containing an AMR gene
        get_unique_reads_for_branched(graph,
                                    os.path.join(args.output_dir,
                                            "debug"))
       # get_path_unique_reads(unitig,
        #                    os.path.join(args.output_dir,
         #                       "debug"))
        # enumerate paths
        #nodePaths = enumerate_paths(graph,
         #                       os.path.join(args.output_dir,
          #                      "debug"))
        # write a debug graph
        graph_data = graph.generate_debug_gml(os.path.join(args.output_dir,
                                                        "debug",
                                                        "debug_gene_mer_graph"),
                                                args.geneMer_size,
                                                args.node_min_coverage,
                                                args.edge_min_coverage,
                                                nodePaths)
        sys.exit(0)
    # write the graph
    sys.stderr.write("\nAmira: writing gene-mer graph\n")
    graph.generate_gml(os.path.join(args.output_dir, "gene_mer_graph"),
                    args.geneMer_size,
                    args.node_min_coverage,
                    args.edge_min_coverage)
    # initialise the UnitigBuilder class
    unitigTools = UnitigTools(graph,
                            genesOfInterest)
    # generate a visualisation of the unitigs
    sys.stderr.write("\nAmira: separating paralog reads\n")
    readFiles, unitig_mapping = unitigTools.separate_reads(args.output_dir,
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