import argparse
import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks, savgol_filter

from amira_prototype.construct_graph import GeneMerGraph
from amira_prototype.construct_unitig import parse_fastq
from amira_prototype.pre_process import convert_pandora_output, process_pandora_json

# from test_functions import TestUnitigTools, Unitig


def get_options() -> argparse.Namespace:
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
        default=None,
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
        "--gene-path",
        dest="path_to_interesting_genes",
        help="path to a newline delimited file of genes of interest",
        required=True,
    )
    parser.add_argument(
        "--debug", dest="debug", action="store_true", default=False, help="Amira debugging"
    )
    args = parser.parse_args()
    return args


def find_trough(node_coverages, filename):
    # Calculate the frequency of each coverage value
    coverages = {}
    for cov in node_coverages:
        coverages[cov] = coverages.get(cov, 0) + 1
    # Sort the coverage values and their frequencies
    vals = sorted(coverages.items())
    x_values = [v[0] for v in vals]
    y_values = [v[1] for v in vals]
    # Apply log transformation to the counts, adding 1 to avoid log(0)
    log_counts = np.log(np.array(y_values) + 1)
    # Smooth the log-transformed histogram counts using a Savitzky-Golay filter
    window_length, poly_order = 15, 3  # Example values; need to be chosen based on your data
    if len(log_counts) < window_length:  # Ensure we have enough data points for the chosen window
        window_length = (
            len(log_counts) // 2 * 2 + 1
        )  # Make the window length the next odd number less than the data length
    smoothed_log_counts = savgol_filter(log_counts, window_length, poly_order)
    plt.figure(figsize=(10, 6))
    plt.bar(x_values, log_counts, label="Counts", color="white", edgecolor="black")
    plt.plot(
        x_values, smoothed_log_counts, color="red", label="Smoothed counts"
    )  # Exponentiate to undo log transform for plotting
    plt.title("Histogram of node coverages with Smoothed Curve")
    plt.xlabel("Node Coverage")
    plt.ylabel("Log of absolute frequency")
    plt.xlim([0, max(x_values)])
    plt.legend()
    plt.savefig(filename)
    plt.close()

    # Check if the first or last data point is a peak and handle it
    modified_smoothed_log_counts = np.array([0] + list(smoothed_log_counts))
    # Find the indices of the peaks
    peak_indices, _ = find_peaks(modified_smoothed_log_counts)
    # Find the two most prominent peaks
    if len(peak_indices) < 2:
        raise ValueError("Not enough peaks to find a trough.")
    prominent_peaks = peak_indices[np.argsort(modified_smoothed_log_counts[peak_indices])[-2:]]
    prominent_peaks.sort()
    # Find the index of the minimum value (trough) between the two peaks
    trough_index = (
        np.argmin(modified_smoothed_log_counts[prominent_peaks[0] : prominent_peaks[1]])
        + prominent_peaks[0]
    )
    trough_value = x_values[trough_index - 1]
    # Plot the histogram and the trough
    plt.figure(figsize=(10, 6))
    plt.bar(x_values, log_counts, label="Counts", color="white", edgecolor="black")
    plt.plot(
        x_values, smoothed_log_counts, color="red", label="Smoothed counts"
    )  # Exponentiate to undo log transform for plotting
    plt.axvline(x=trough_value, color="r", linestyle="--", label=f"Trough at x={trough_value:.2f}")
    plt.title("Histogram of node coverages with Smoothed Curve")
    plt.xlabel("Node Coverage")
    plt.ylabel("Log of absolute frequency")
    plt.xlim([0, max(x_values)])
    plt.legend()
    plt.savefig(filename)
    plt.close()
    return trough_value


def plot_log_histogram(distances, filename):
    import random

    distances = random.sample(distances, 10000)
    # Calculate histogram data (counts and bin edges) without plotting
    counts, bin_edges = np.histogram(
        distances, bins=len(distances) - 1
    )  # Adjust bin count as needed
    # Calculate the log of counts to plot, adding 1 to avoid log(0)
    log_counts = np.log(counts + 1)
    # Calculate bin centers from bin edges for bar placement
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    # Create the bar plot
    plt.figure(figsize=(10, 6))  # Optional: Adjust figure size
    # Plot bars with bin centers as x values and logged counts as heights
    plt.bar(
        bin_centers,
        log_counts,
        align="center",
        width=np.diff(bin_edges),
        color="white",
        edgecolor="black",
    )
    # Set plot title and labels
    plt.title("Histogram of distances between genes")
    plt.xlabel("Distance")
    plt.ylabel("Log of absolute frequency")
    # Save the plot as a PNG file
    plt.savefig(filename, dpi=300)  # Adjust dpi for higher resolution if needed
    plt.close()


def write_debug_files(
    annotatedReads: dict[str, list[str]],
    geneMer_size: int,
    genesOfInterest: list[str],
    output_dir: str,
) -> GeneMerGraph:
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
    import statistics

    unitigs_in_component = []
    node_to_unitig_mapping = {}
    seen_unitigs = set()
    for node_hash in unitig_boundaries:
        node = node_mapping[node_hash]
        for neighbour_node in graph.get_all_neighbors(node):
            assert node_hash == node.__hash__()
            if not graph.get_degree(neighbour_node) > 2:
                linear_path = graph.get_linear_path_for_node(neighbour_node, True)
            else:
                if (
                    graph.get_edges_between_nodes(node, neighbour_node)[0].get_sourceNodeDirection()
                    == 1
                ):
                    linear_path = [node_hash, neighbour_node.__hash__()]
                if (
                    graph.get_edges_between_nodes(node, neighbour_node)[0].get_sourceNodeDirection()
                    == -1
                ):
                    linear_path = [neighbour_node.__hash__(), node_hash]
            assert linear_path[0] in unitig_boundaries and linear_path[-1] in unitig_boundaries
            list_of_genes = (
                graph.follow_path_to_get_annotations(linear_path, node.get_list_of_reads()[0])
                if not len(linear_path) == 1
                else graph.get_gene_mer_label(graph.get_node_by_hash(linear_path[0])).split("~~~")
            )
            if not (
                tuple(linear_path) in seen_unitigs
                or tuple(list(reversed(linear_path))) in seen_unitigs
            ):
                # if len(linear_path) == 1:
                #     print(component, linear_path, node_hash, neighbour_node.__hash__())
                #     print([n.__hash__() for n in graph.get_all_neighbors(node)], [n.__hash__() for n in graph.get_all_neighbors(neighbour_node)])
                path_edge_coverages = [
                    graph.get_edges_between_nodes(
                        graph.get_node_by_hash(linear_path[i]),
                        graph.get_node_by_hash(linear_path[i + 1]),
                    )[0].get_edge_coverage()
                    for i in range(len(linear_path) - 1)
                ]
                unitig = Unitig(
                    [node_mapping[h] for h in linear_path],
                    list_of_genes,
                    component,
                    current_unitig_id,
                    statistics.mean(path_edge_coverages),
                )
                for h in linear_path:
                    node_to_unitig_mapping[h] = unitig
                unitigs_in_component.append(unitig)
                current_unitig_id += 1
                seen_unitigs.add(tuple(linear_path))
    return unitigs_in_component, node_to_unitig_mapping, current_unitig_id


def unitig_list_of_genes(unitigs_in_component):
    fw_genes_per_unitig = {}
    rv_genes_per_unitig = {}
    for i in range(len(unitigs_in_component)):
        fw_genes_per_unitig[i] = unitigs_in_component[i].get_genes()
        rv_genes_per_unitig[i] = unitigs_in_component[i].get_reverse_genes()
    return fw_genes_per_unitig, rv_genes_per_unitig


def calculate_unitig_distance_matrix(unitigs_in_component):
    from collections import Counter

    from tqdm import tqdm

    fw_genes_per_unitig, rv_genes_per_unitig = unitig_list_of_genes(unitigs_in_component)

    unitig_correction = {}
    for i in tqdm(fw_genes_per_unitig):
        distance_matrix = [0 for _ in range(len(unitigs_in_component))]
        for j in fw_genes_per_unitig:
            if i != j:
                sorted_matches = sorted(
                    [fw_genes_per_unitig[j][:], rv_genes_per_unitig[j][:]],
                    key=lambda x: len(set(fw_genes_per_unitig[i]).intersection(set(x))),
                    reverse=True,
                )
                counter1 = Counter(fw_genes_per_unitig[i])
                counter2 = Counter(sorted_matches[0])
                shared_count = len(list((counter1 & counter2).elements()))
                union_count = len(fw_genes_per_unitig[i])
                distance_matrix[j] = shared_count / union_count
        with open("dist_mat.txt", "a+") as o:
            o.write("\t".join([str(round(c, 1)) for c in distance_matrix]) + "\n")
        max_index, max_value = max(enumerate(distance_matrix), key=lambda x: x[1])
        # only correct unitigs if there is another one that is more than 80% similar
        if (
            max_value > 0.8
            and unitigs_in_component[max_index].get_coverage()
            >= unitigs_in_component[i].get_coverage()
        ):
            unitig_correction[unitigs_in_component[i]] = unitigs_in_component[max_index]
        else:
            unitig_correction[unitigs_in_component[i]] = unitigs_in_component[i]
    return unitig_correction


def split_unitigs(unitigs_in_component):
    true_unitigs = []
    false_unitigs = []
    for unitig in unitigs_in_component:
        if unitig.get_mean_edge_coverage() < 10:
            false_unitigs.append(unitig)
        else:
            true_unitigs.append(unitig)
    return true_unitigs, false_unitigs


def get_unitigs(graph, output_dir):
    from tqdm import tqdm

    """Main function to extract unitigs from a graph and handle corrections."""
    current_unitig_id = 0
    reads_as_unitigs = {}
    for component in tqdm(graph.components()):
        nodes_in_component = graph.get_nodes_in_component(component)
        unitig_boundaries0, node_mapping = get_unitig_boundaries(nodes_in_component, graph)
        # get the unitigs in the component
        unitigs_in_component, node_to_unitig_mapping, current_unitig_id = generate_unitigs(
            node_mapping, unitig_boundaries, graph, component, current_unitig_id
        )
        if component == 5:
            true_unitigs, false_unitigs = split_unitigs(unitigs_in_component)
            # make a dictionary to correct obviously false unitigs
            unitig_correction = correct_unitigs(true_unitigs, false_unitigs)
            print(unitig_correction)


def correct_unitigs(true_unitigs, false_unitigs):
    from collections import Counter

    from tqdm import tqdm

    fw_genes_true_unitigs, rv_genes_true_unitigs = unitig_list_of_genes(true_unitigs)
    for unitig in tqdm(false_unitigs):
        fw_genes = unitig.get_genes()
        distances = []
        for i in range(len(true_unitigs)):
            sorted_matches = sorted(
                [fw_genes_true_unitigs[i][:], rv_genes_true_unitigs[i][:]],
                key=lambda x: len(set(fw_genes).intersection(set(x))),
                reverse=True,
            )
            counter1 = Counter(fw_genes)
            counter2 = Counter(sorted_matches[0])
            shared_count = len(list((counter1 & counter2).elements()))
            union_count = len(fw_genes_true_unitigs[i])
            distances.append(shared_count / len(list((counter1 | counter2).elements())))
        max_index, max_value = max(enumerate(distances), key=lambda x: x[1])
        print(max_value, fw_genes, fw_genes_true_unitigs[max_index])


def main() -> None:
    # get command line options
    args = get_options()
    # make the output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # import the list of genes of interest
    with open(args.path_to_interesting_genes, "r") as i:
        genesOfInterest = i.read().splitlines()
    # remove invalid characters from the AMR gene headers
    cleaned_genesOfInterest = set()
    chars_to_remove = ".|()-*+#:=/,'"
    for g in genesOfInterest:
        cleaned_gene = "".join(char for char in g if char not in chars_to_remove)
        cleaned_genesOfInterest.add(cleaned_gene)
    genesOfInterest = cleaned_genesOfInterest
    # import a JSON of genes on reads
    if args.pandoraJSON:
        annotatedReads, genesOfInterest = process_pandora_json(args.pandoraJSON, genesOfInterest)
    # load the pandora consensus and convert to a dictionary
    if args.pandoraSam:
        pandora_consensus = parse_fastq(args.pandoraConsensus)
        annotatedReads, genesOfInterest, distances = convert_pandora_output(
            args.pandoraSam, pandora_consensus, genesOfInterest, args.gene_min_coverage
        )
        # plot_log_histogram(distances, os.path.join(args.output_dir, "distance_between_genes.png"))
    print(list(sorted(list(genesOfInterest))))
    # write out debug files if specified
    if args.debug:
        write_debug_files(annotatedReads, args.geneMer_size, genesOfInterest, args.output_dir)
    # build the gene-mer graph
    graph = GeneMerGraph(annotatedReads, args.geneMer_size)
    # get rid of low coverage components
    graph.remove_low_coverage_components(5)
    # correct on a per component basis
    new_annotatedReads = {}
    for component in graph.components():
        sys.stderr.write(
            f"\nAmira: correcting reads for component {component}/{len(graph.components())}\n"
        )
        # get the nodes in the component
        nodes_in_component = graph.get_nodes_in_component(component)
        # get the read annotations for the reads in this component
        component_annotated_reads = {}
        for node in nodes_in_component:
            for read in node.get_reads():
                if not read in component_annotated_reads:
                    component_annotated_reads[read] = annotatedReads[read]
        # build a graph of just this component
        component_graph = GeneMerGraph(component_annotated_reads, args.geneMer_size)
        # remove short linear paths
        component_graph.remove_short_linear_paths(20)
        component_annotated_reads = component_graph.correct_reads()
        component_graph = GeneMerGraph(component_annotated_reads, args.geneMer_size)
        # dynamically determine the node threshold for filtering
        if not args.node_min_coverage:
            try:
                node_min_coverage = find_trough(
                    component_graph.get_all_node_coverages(),
                    os.path.join(args.output_dir, f"node_coverages_{component}.png"),
                )
            except:
                node_min_coverage = 5
        else:
            node_min_coverage = args.node_min_coverage
        component_graph.filter_graph(node_min_coverage, 1)
        component_annotated_reads = component_graph.correct_reads()
        # remove low coverage components
        component_graph = GeneMerGraph(component_annotated_reads, args.geneMer_size)
        component_graph.remove_low_coverage_components(5)
        component_annotated_reads = component_graph.correct_reads()
        # collect the component read annotations
        for read_id in component_annotated_reads:
            if read_id in new_annotatedReads:
                print(read_id)
                print(new_annotatedReads[read_id])
                print(component_annotated_reads[read_id])
            new_annotatedReads[read_id] = component_annotated_reads[read_id]
    sys.stderr.write("\nAmira: building corrected gene-mer graph\n")
    graph = GeneMerGraph(new_annotatedReads, args.geneMer_size)
    # color nodes in the graph if --debug is used
    if args.debug:
        for node in graph.all_nodes():
            node.color_node(genesOfInterest)
    # write out the graph as a GML
    sys.stderr.write("\nAmira: writing gene-mer graph\n")
    graph.generate_gml(
        os.path.join(args.output_dir, "gene_mer_graph"),
        args.geneMer_size,
        None,
        args.edge_min_coverage,
    )
    # assign reads to AMR genes by path
    sys.stderr.write("\nAmira: clustering reads\n")
    clusters = graph.assign_reads_to_genes(genesOfInterest)
    # write out the clustered reads
    with open(os.path.join(args.output_dir, "reads_per_amr_gene.json"), "w") as o:
        o.write(json.dumps(clusters))
    sys.exit(0)


if __name__ == "__main__":
    main()
