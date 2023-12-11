import argparse
import json
import os
import sys

import pysam
from tqdm import tqdm

from amira_prototype.construct_graph import GeneMerGraph
from amira_prototype.construct_unitig import UnitigTools, parse_fastq
from amira_prototype.test_functions import TestUnitigTools, Unitig


def get_options():
    """define args from the command line"""
    parser = argparse.ArgumentParser(description="Build a prototype gene de Bruijn graph.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pandoraSam", dest="pandoraSam", help="Pandora map SAM file path")
    group.add_argument("--pandoraJSON", dest="pandoraJSON", help="Pandora map JSON file path")
    parser.add_argument(
        "--pandoraConsensus",
        dest="pandoraConsensus",
        help="path to Pandora consensus fastq",
        required=False,
    )
    parser.add_argument(
        "--readfile", dest="readfile", help="path of gzipped long read fastq", required=True
    )
    parser.add_argument(
        "--output",
        dest="output_dir",
        type=str,
        default="gene_de_Bruijn_graph",
        help="directory for Amira outputs",
    )
    parser.add_argument(
        "-k",
        dest="geneMer_size",
        type=int,
        default=3,
        help="kmer length for the gene de Bruijn graph",
    )
    parser.add_argument(
        "-n",
        dest="node_min_coverage",
        type=int,
        default=1,
        help="minimum threshold for gene-mer coverage",
    )
    parser.add_argument(
        "-e",
        dest="edge_min_coverage",
        type=int,
        default=1,
        help="minimum threshold for edge coverage between gene-mers",
    )
    parser.add_argument(
        "-g",
        dest="gene_min_coverage",
        type=int,
        default=1,
        help="minimum threshold for gene filtering",
    )
    parser.add_argument(
        "-p",
        dest="bubble_popper_threshold",
        type=float,
        default=1.5,
        help="minimum difference in coverage threshold to collapse paths in the graph (default = 1.5)",
    )
    parser.add_argument(
        "-c",
        dest="cleaning_iterations",
        type=int,
        default=2,
        help="number of gene-mer graph cleaning operations (default = 2)",
    )
    parser.add_argument(
        "--gene-path",
        dest="path_to_interesting_genes",
        help="path to a newline delimited file of genes of interest",
        required=True,
    )
    parser.add_argument(
        "--flye-path", dest="flye_path", help="path to Flye binary", default=None, required=False
    )
    parser.add_argument(
        "--raven-path", dest="raven_path", help="path to Raven binary", default=None, required=False
    )
    parser.add_argument(
        "--use-consensus",
        dest="use_pandora_consensus",
        help="polish the pandora consensus of each gene to recover AMR alleles",
        action="store_true",
        default=None,
        required=False,
    )
    parser.add_argument(
        "--racon-path", dest="racon_path", help="path to Racon binary", default=None, required=False
    )
    parser.add_argument(
        "--threads", dest="threads", type=int, default=1, help="number of threads to use"
    )
    parser.add_argument(
        "--debug", dest="debug", action="store_true", default=False, help="Amira debugging"
    )
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


def get_read_start(cigar):
    """return an int of the 0 based position where the read region starts mapping to the gene"""
    # check if there are any hard clipped bases at the start of the mapping
    if cigar[0][0] == 5:
        regionStart = cigar[0][1]
    else:
        regionStart = 0
    return regionStart


def get_read_end(cigar, regionStart):
    """return an int of the 0 based position where the read region stops mapping to the gene"""
    regionLength = 0
    for tuple in cigar:
        if not tuple[0] == 5:
            regionLength += tuple[1]
    regionEnd = regionStart + regionLength
    return regionEnd, regionLength


def determine_gene_strand(read):
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


def convert_pandora_output(pandoraSam, pandora_consensus, genesOfInterest, geneMinCoverage):
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
            regionEnd, regionLength = get_read_end(cigar, regionStart)
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
    for t in to_delete:
        del annotatedReads[t]
    assert not len(annotatedReads) == 0
    return annotatedReads, list(subsettedGenesOfInterest)


def write_debug_files(annotatedReads, geneMer_size, genesOfInterest, output_dir):
    raw_graph = GeneMerGraph(annotatedReads, geneMer_size)
    # color nodes in the graph
    for node in raw_graph.all_nodes():
        node.color_node(genesOfInterest)
    import json

    with open(os.path.join(output_dir, "genesAnnotatedOnReads.json"), "w") as o:
        o.write(json.dumps(annotatedReads))
    sys.stderr.write("\nAmira: writing pre-correction gene-mer graph\n")
    raw_graph.generate_gml(
        os.path.join(output_dir, "pre_correction_gene_mer_graph"), geneMer_size, 1, 1
    )
    return raw_graph


def make_unitig_plot(readsAsUnitigs, output_dir):
    readsAsUnitigIds = {}
    for read in readsAsUnitigs:
        readsAsUnitigIds[read] = ["+" + str(u.get_unitig_ID()) for u in readsAsUnitigs[read]]
    unitig_graph = GeneMerGraph(readsAsUnitigIds, 1)
    unitig_graph.generate_gml(os.path.join(output_dir, "unitig_graph.gml"), 1, 1, 1)


def shared_indices(list1, list2):
    indices = [i for i, elem in enumerate(list1) if elem in list2]
    return (indices[0], indices[-1]) if indices else (None, None)


def calculate_unitig_distance_matrix(unitigs_in_component):
    from collections import Counter

    fw_genes_per_unitig = {}
    rv_genes_per_unitig = {}
    for i in range(len(unitigs_in_component)):
        fw_genes_per_unitig[i] = unitigs_in_component[i].get_genes()
        rv_genes_per_unitig[i] = unitigs_in_component[i].get_reverse_genes()
    unitig_correction = {}
    for i in tqdm(fw_genes_per_unitig):
        distance_matrix = [0 for _ in range(len(unitigs_in_component))]
        for j in fw_genes_per_unitig:
            if not i == j:
                sorted_matches = sorted(
                    [fw_genes_per_unitig[j][:], rv_genes_per_unitig[j][:]],
                    key=lambda x: len(set(fw_genes_per_unitig[i]).intersection(set(x))),
                    reverse=True,
                )
                counter1 = Counter(fw_genes_per_unitig[i])
                counter2 = Counter(sorted_matches[0])
                shared_count = len(list((counter1 & counter2).elements()))
                union_count = len(
                    fw_genes_per_unitig[i]
                )  # len(list((counter1 | counter2).elements()))
                distance_matrix[j] = shared_count / union_count
        with open("dist_mat.txt", "a+") as o:
            o.write("\t".join([str(round(c, 1)) for c in distance_matrix]) + "\n")
        max_index, max_value = max(enumerate(distance_matrix), key=lambda x: x[1])
        # only correct unitigs if there is another one that is more than 80% similar
        if (
            max_value > 0.7
            and unitigs_in_component[max_index].get_coverage()
            >= unitigs_in_component[i].get_coverage()
        ):
            unitig_correction[unitigs_in_component[i]] = unitigs_in_component[max_index]
        else:
            unitig_correction[unitigs_in_component[i]] = unitigs_in_component[i]
    return unitig_correction


def get_unitig_boundaries(nodes_in_component, graph):
    """Identify unitig boundaries within a component based on node degrees."""
    unitig_boundaries = set()
    node_mapping = {}
    for node in nodes_in_component:
        node_hash = node.__hash__()
        node_mapping[node_hash] = node
        if not graph.get_degree(node) > 2:
            forward_neighbors = graph.get_forward_neighbors(node)
            backward_neighbors = graph.get_backward_neighbors(node)
            if len(backward_neighbors) == 0 or len(forward_neighbors) == 0:
                unitig_boundaries.add(node_hash)
        else:
            unitig_boundaries.add(node_hash)
    return unitig_boundaries, node_mapping


def generate_unitigs(node_mapping, unitig_boundaries, graph, component, current_unitig_id):
    """Generate unitigs based on boundaries and graph structure."""
    unitigs_in_component = []
    node_to_unitig_mapping = {}
    seen_unitigs = set()
    for node_hash in unitig_boundaries:
        node = node_mapping[node_hash]
        for neighbour_node in graph.get_all_neighbors(node):
            linear_path = graph.get_linear_path_for_node(neighbour_node, True)
            assert linear_path[0] in unitig_boundaries and linear_path[-1] in unitig_boundaries
            list_of_genes = (
                graph.follow_path_to_get_annotations(linear_path)
                if not len(linear_path) == 1
                else graph.get_gene_mer_label(graph.get_node_by_hash(linear_path[0])).split("~~~")
            )
            if not (
                tuple(linear_path) in seen_unitigs
                or tuple(list(reversed(linear_path))) in seen_unitigs
            ):
                unitig = Unitig(
                    [node_mapping[h] for h in linear_path],
                    list_of_genes,
                    component,
                    current_unitig_id,
                )
                for node_hash in linear_path:
                    node_to_unitig_mapping[node_hash] = unitig
                unitigs_in_component.append(unitig)
                current_unitig_id += 1
                seen_unitigs.add(tuple(linear_path))
    return unitigs_in_component, node_to_unitig_mapping, current_unitig_id


def get_reads_to_correct(correct_per_unitig):
    """Get the set of reads that need to be corrected based on unitig correction."""
    reads_to_correct = set()
    for unitig in correct_per_unitig:
        for read in unitig.get_reads():
            reads_to_correct.add(read)
    return reads_to_correct


def update_reads_as_unitigs(graph, reads_to_correct, node_to_unitig_mapping, reads_as_unitigs):
    """Update reads to unitigs mapping based on corrections."""
    for read in reads_to_correct:
        node_hashes = graph.get_readNodes()[read]
        unitigs_on_read, previous_unitig = [], ""
        for node_hash in node_hashes:
            if node_hash in node_to_unitig_mapping:
                unitig = node_to_unitig_mapping[node_hash]
                if not unitig == previous_unitig:
                    unitigs_on_read.append(unitig)
                    previous_unitig = unitig
        if not unitigs_on_read == []:
            reads_as_unitigs[read] = unitigs_on_read
    return reads_as_unitigs


def get_unitigs(graph, output_dir):
    """Main function to extract unitigs from a graph and handle corrections."""
    current_unitig_id = 0
    reads_as_unitigs = {}
    for component in tqdm(graph.components()):
        nodes_in_component = graph.get_nodes_in_component(component)
        unitig_boundaries, node_mapping = get_unitig_boundaries(nodes_in_component, graph)
        # get the unitigs in the component
        unitigs_in_component, node_to_unitig_mapping, current_unitig_id = generate_unitigs(
            node_mapping, unitig_boundaries, graph, component, current_unitig_id
        )
        # make a dictionary to correct obviously false unitigs
        unitig_correction = calculate_unitig_distance_matrix(unitigs_in_component)
        # correct the obviously false unitigs
        correct_per_unitig = {
            u: unitig_correction[u] for u in unitig_correction if not unitig_correction[u] == u
        }
        # update the node to unitig
        # for node_hash in node_to_unitig_mapping:
        #    if node_to_unitig_mapping[node_hash] in correct_per_unitig:
        #        node_to_unitig_mapping[node_hash] = correct_per_unitig[node_to_unitig_mapping[node_hash]]
        reads_to_correct = get_reads_to_correct(correct_per_unitig)
        reads_as_unitigs = update_reads_as_unitigs(
            graph, reads_to_correct, node_to_unitig_mapping, reads_as_unitigs
        )
    make_unitig_plot(reads_as_unitigs, output_dir)


def get_reads_that_span_full_component(graph):
    longest_reads = {}
    for component in tqdm(graph.components()):
        node_mapping = {}
        component_boundaries = set()
        nodes_in_component = graph.get_nodes_in_component(component)
        for node in nodes_in_component:
            node_hash = node.__hash__()
            node_mapping[node_hash] = node
            if not graph.get_degree(node) > 2:
                forward_neighbors = graph.get_forward_neighbors(node)
                backward_neighbors = graph.get_backward_neighbors(node)
                if len(backward_neighbors) == 0 or len(forward_neighbors) == 0:
                    component_boundaries.add(node_hash)
        # get the reads that have at least 2 component boundaries
        boundary_reads = []
        seen_reads = set()
        for nodeHash in component_boundaries:
            for read in node_mapping[nodeHash].get_reads():
                if not read in seen_reads:
                    boundaryNodesOnRead = [
                        n for n in graph.get_readNodes()[read] if n in component_boundaries
                    ]
                    if len(boundaryNodesOnRead) > 0:
                        boundary_reads.append((read, graph.get_readNodes()[read]))
                        seen_reads.add(read)
        boundary_reads = sorted(boundary_reads, key=lambda x: len(x[1]), reverse=True)[:50]
        for r in boundary_reads:
            longest_reads[r[0]] = graph.follow_path_to_get_annotations(r[1])
    long_graph = GeneMerGraph(longest_reads, 5)
    long_graph.generate_gml("amira.output/long_read_graph", 5, 1, 1)


def main():
    # get command line options
    args = get_options()
    # make the output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # import the list of genes of interest
    with open(args.path_to_interesting_genes, "r") as i:
        genesOfInterest = i.read().splitlines()
    cleaned_genesOfInterest = []
    chars_to_remove = ".|()-*+#:=/,'"
    for g in genesOfInterest:
        cleaned_gene = "".join(char for char in g if char not in chars_to_remove)
        cleaned_genesOfInterest.append(cleaned_gene)  # convert the Pandora SAM file to a dictionary
    genesOfInterest = cleaned_genesOfInterest
    if args.pandoraJSON:
        annotatedReads, genesOfInterest = process_pandora_json(args.pandoraJSON, genesOfInterest)
    if args.pandoraSam:
        # load the pandora consensus and convert to a dictionary
        pandora_consensus = parse_fastq(args.pandoraConsensus)
        annotatedReads, genesOfInterest = convert_pandora_output(
            args.pandoraSam, pandora_consensus, genesOfInterest, args.gene_min_coverage
        )
    print(list(sorted(list(genesOfInterest))))
    if args.debug:
        raw_graph = write_debug_files(
            annotatedReads, args.geneMer_size, genesOfInterest, args.output_dir
        )
        short_annotated_reads = raw_graph._shortReads
        short_graph = GeneMerGraph(short_annotated_reads, 3)
        # short_graph.generate_gml(os.path.join(args.output_dir, "short_graph"),
        #                        3,
        #                        1,
        #                        1)
    # remove components that have a coverage of 1
    # raw_graph.remove_low_coverage_components(3)
    # # get the new read annotations
    # readNodes = raw_graph.get_readNodes()
    # annotatedReads = {}
    # for readId in readNodes:
    #     annotatedReads[readId] = raw_graph.follow_path_to_get_annotations(readNodes[readId])
    # # iteratively replace the lowest coverage path in a component with the highest coverage path
    # for i in range(3):
    #     # build the updated graph without the removed nodes
    #     graph = GeneMerGraph(annotatedReads,
    #                     args.geneMer_size)
    #     # remove components that have a coverage of 1
    #     graph.remove_low_coverage_components(3)
    #     for component in graph.components():
    #         # get the heaviest path through each component
    #         heaviest_path = graph.new_get_heaviest_path_through_component(component)
    #         # get the nodes that are not in the heaviest path
    #         reads_to_correct = set()
    #         path_set = set(heaviest_path)
    #         for node in graph.get_nodes_in_component(component):
    #             node_hash = node.__hash__()
    #             if not node_hash in path_set:
    #                 for readId in node.get_reads():
    #                     reads_to_correct.add(readId)
    #         # correct the reads to the heaviest path
    #         graph.correct_read_nodes_to_heaviest_path(heaviest_path, reads_to_correct)
    #     # get the new read annotations
    #     readNodes = graph.get_readNodes()
    #     annotatedReads = {}
    #     for readId in readNodes:
    #         for i in range(len(readNodes[readId]) - 1):
    #             sourceNode = graph.get_node_by_hash(readNodes[readId][i])
    #             targetNode = graph.get_node_by_hash(readNodes[readId][i+1])
    #             if not graph.check_if_nodes_are_adjacent(sourceNode, targetNode):
    #                 graph.add_edge(sourceNode.get_geneMer(), targetNode.get_geneMer())
    #         annotatedReads[readId] = graph.follow_path_to_get_annotations(readNodes[readId])
    # build the corrected graph
    sys.stderr.write("\nAmira: building corrected gene-mer graph\n")
    graph = GeneMerGraph(annotatedReads, args.geneMer_size)
    # # filter low coverage things
    graph.filter_graph(args.node_min_coverage, args.edge_min_coverage)
    # #graph.remove_low_coverage_components(3)
    # # get the new read annotations
    # readNodes = graph.get_readNodes()
    # annotatedReads = {}
    # for readId in readNodes:
    #     for i in range(len(readNodes[readId]) - 1):
    #         sourceNode = graph.get_node_by_hash(readNodes[readId][i])
    #         targetNode = graph.get_node_by_hash(readNodes[readId][i+1])
    #         if not graph.check_if_nodes_are_adjacent(sourceNode, targetNode):
    #             graph.add_edge(sourceNode.get_geneMer(), targetNode.get_geneMer())
    #     annotatedReads[readId] = graph.follow_path_to_get_annotations(readNodes[readId])
    # graph = GeneMerGraph(annotatedReads,
    #                     1)
    # graph.generate_gml(os.path.join(args.output_dir, "gene_mer_graph"),
    #                 1,
    #                 1,
    #                 1)
    # annotatedReads = graph.cut_triangles()
    # graph = GeneMerGraph(annotatedReads,
    #                     1)
    # graph.generate_gml(os.path.join(args.output_dir, "gene_mer_graph"),
    #                 1,
    #                 1,
    #                 1)
    # graph = GeneMerGraph(annotatedReads,
    #                    args.geneMer_size)
    # write out the gene mer graph as a gml
    sys.stderr.write("\nAmira: writing gene-mer graph\n")
    if args.debug:
        # color nodes in the graph
        for node in graph.all_nodes():
            node.color_node(genesOfInterest)
    graph.generate_gml(
        os.path.join(args.output_dir, "gene_mer_graph"),
        args.geneMer_size,
        args.node_min_coverage,
        args.edge_min_coverage,
    )
    # get_reads_that_span_full_component(graph)
    # return a list of unitigs in the graph
    # unitigs = get_unitigs(graph,
    #                    args.output_dir)
    # initialise the UnitigBuilder class
    unitigTools = TestUnitigTools(graph, genesOfInterest, args.readfile, args.output_dir)
    # unitigTools.assign_reads_to_amr_genes()

    sys.exit(0)


if __name__ == "__main__":
    main()
