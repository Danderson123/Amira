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
    #graph.remove_low_coverage_components(5)
    # remove short linear paths
    #graph.remove_short_linear_paths(20)
    new_annotatedReads = graph.correct_reads()
    #graph = GeneMerGraph(new_annotatedReads, args.geneMer_size)
    # dynamically determine the node threshold for filtering
    # if not args.node_min_coverage:
    #     try:
    #         node_min_coverage = find_trough(
    #             graph.get_all_node_coverages(),
    #             os.path.join(args.output_dir, f"node_coverages.png"),
    #         )
    #     except:
    #         node_min_coverage = 5
    # else:
    #     node_min_coverage = args.node_min_coverage
    #graph.filter_graph(node_min_coverage, 1)
    #new_annotatedReads = graph.correct_reads()
    # remove low coverage components
    # graph = GeneMerGraph(new_annotatedReads, args.geneMer_size)
    # graph.remove_low_coverage_components(5)
    # new_annotatedReads = graph.correct_reads()
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
