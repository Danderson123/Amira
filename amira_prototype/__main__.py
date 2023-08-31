import argparse
import json
import os
import statistics
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
    args = parser.parse_args()
    return args

def process_pandora_json(pandoraJSON, genesOfInterest):
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
    for read in to_delete:
        del annotatedReads[read]
    genesOfInterest = subsettedGenesOfInterest
    return annotatedReads, genesOfInterest

def write_debug_files(annotatedReads,
                    geneMer_size,
                    genesOfInterest,
                    output_dir):
    raw_graph = GeneMerGraph(annotatedReads,
                            geneMer_size)
    # color nodes in the graph
    for node in raw_graph.all_nodes():
        node.color_node(genesOfInterest)
    import json
    with open(os.path.join(output_dir, "genesAnnotatedOnReads.json"), "w") as o:
        o.write(json.dumps(annotatedReads))
    sys.stderr.write("\nAmira: writing pre-correction gene-mer graph\n")
    raw_graph.generate_gml(os.path.join(output_dir, "pre_correction_gene_mer_graph"),
                                geneMer_size,
                                1,
                                1)
    return raw_graph

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
    if args.pandoraJSON:
        annotatedReads, genesOfInterest = process_pandora_json(args.pandoraJSON,
                                                            genesOfInterest)
    print(list(sorted(list(genesOfInterest))))
    if args.debug:
        raw_graph = write_debug_files(annotatedReads,
                                    args.geneMer_size,
                                    genesOfInterest,
                                    args.output_dir)
        short_annotated_reads = raw_graph._shortReads
        short_graph = GeneMerGraph(short_annotated_reads,
                                3)
        short_graph.generate_gml(os.path.join(args.output_dir, "short_graph"),
                                3,
                                1,
                                1)
    # remove components that have a coverage of 1
    raw_graph.remove_low_coverage_components(4)
    # get the new read annotations
    readNodes = raw_graph.get_readNodes()
    annotatedReads = {}
    for readId in readNodes:
        annotatedReads[readId] = raw_graph.follow_path_to_get_annotations(readNodes[readId])
    # iteratively replace the lowest coverage path in a component with the highest coverage path
    for i in range(1):
        # build the updated graph without the removed nodes
        graph = GeneMerGraph(annotatedReads,
                        args.geneMer_size)
        # remove components that have a coverage of 1
        graph.remove_low_coverage_components(4)
        for component in graph.components():
            # get the heaviest path through each component
            heaviest_path = graph.new_get_heaviest_path_through_component(component)
            # get the nodes that are not in the heaviest path
            reads_to_correct = set()
            path_set = set(heaviest_path)
            for node in graph.get_nodes_in_component(component):
                node_hash = node.__hash__()
                if not node_hash in path_set:
                    for readId in node.get_reads():
                        reads_to_correct.add(readId)
            # correct the reads to the heaviest path
            graph.correct_read_nodes_to_heaviest_path(heaviest_path, reads_to_correct)
            #heaviest_path, reverse_heaviest_path = graph.get_heaviest_path_through_component(component)
            #lightest_path, heaviest_path = graph.get_lightest_path_through_component(heaviest_path, reverse_heaviest_path)
            #graph.correct_reads_to_heaviest_path(heaviest_path,
            #                                lightest_path)
            #lightest_path_from_start, lightest_coverages_from_start = graph.new_get_lightest_path_through_component(heaviest_path)
            #lightest_path_from_end, lightest_coverages_from_end = graph.new_get_lightest_path_through_component(reverse_heaviest_path)
            #lightest_path_from_end = list(reversed(lightest_path_from_end))
            # decide which lightest path we are correcting based on their coverages
            #lightest_coverages_from_start = statistics.mean(lightest_coverages_from_start)
            #lightest_coverages_from_end = statistics.mean(lightest_coverages_from_end)
            #if lightest_coverages_from_start < lightest_coverages_from_end:
            #    graph.correct_reads_to_heaviest_path(heaviest_path,
            #                                lightest_path_from_start)
            #elif lightest_coverages_from_start > lightest_coverages_from_end:
            #    graph.correct_reads_to_heaviest_path(heaviest_path,
            #                                lightest_path_from_end)
            #else:
            #    pass
        # get the new read annotations
        readNodes = graph.get_readNodes()
        annotatedReads = {}
        for readId in readNodes:
            for i in range(len(readNodes[readId]) - 1):
                sourceNode = graph.get_node_by_hash(readNodes[readId][i])
                targetNode = graph.get_node_by_hash(readNodes[readId][i+1])
                if not graph.check_if_nodes_are_adjacent(sourceNode, targetNode):
                    graph.add_edge(sourceNode.get_geneMer(), targetNode.get_geneMer())
            annotatedReads[readId] = graph.follow_path_to_get_annotations(readNodes[readId])
    # build the corrected graph
    sys.stderr.write("\nAmira: building corrected gene-mer graph\n")
    graph = GeneMerGraph(annotatedReads,
                        args.geneMer_size)
    # filter low coverage things
    graph.filter_graph(args.node_min_coverage,
                    args.edge_min_coverage)
    graph.remove_low_coverage_components(4)
    # write out the gene mer graph as a gml
    sys.stderr.write("\nAmira: writing gene-mer graph\n")
    if args.debug:
        # color nodes in the graph
        for node in graph.all_nodes():
            node.color_node(genesOfInterest)
    graph.generate_gml(os.path.join(args.output_dir, "gene_mer_graph"),
                    args.geneMer_size,
                    args.node_min_coverage,
                    args.edge_min_coverage)
    # initialise the UnitigBuilder class
    unitigTools = TestUnitigTools(graph,
                            genesOfInterest,
                            args.readfile,
                            args.output_dir)
    unitigTools.assign_reads_to_amr_genes()
    sys.exit(0)

if __name__ == "__main__":
    main()