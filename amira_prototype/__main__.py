import argparse
from joblib import Parallel, delayed
import json
import os
import shutil
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks, savgol_filter
from tqdm import tqdm

from amira_prototype.construct_graph import GeneMerGraph, build_graph, merge_graphs
from amira_prototype.construct_unitig import parse_fastq, write_fastq
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
        help="Path to Pandora consensus fastq.",
        required=False,
    )
    parser.add_argument(
        "--readfile", dest="readfile", help="path of gzipped long read fastq.", required=True
    )
    parser.add_argument(
        "--output",
        dest="output_dir",
        type=str,
        default="gene_de_Bruijn_graph",
        help="Directory for Amira outputs.",
    )
    parser.add_argument(
        "-k",
        dest="geneMer_size",
        type=int,
        default=3,
        help="k-mer length for the gene de Bruijn graph.",
    )
    parser.add_argument(
        "-n",
        dest="node_min_coverage",
        type=int,
        default=None,
        help="Minimum threshold for gene-mer coverage.",
    )
    parser.add_argument(
        "-e",
        dest="edge_min_coverage",
        type=int,
        default=1,
        help="Minimum threshold for edge coverage between gene-mers.",
    )
    parser.add_argument(
        "-g",
        dest="gene_min_coverage",
        type=int,
        default=1,
        help="Minimum threshold for gene filtering.",
    )
    parser.add_argument(
        "--minimum-length-proportion",
        dest="lower_gene_length_threshold",
        type=float,
        default=0.5,
        help="Minimum length threshold for gene filtering.",
    )
    parser.add_argument(
        "--maximum-length-proportion",
        dest="upper_gene_length_threshold",
        type=float,
        default=1.5,
        help="Maximum length threshold for gene filtering.",
    )
    parser.add_argument(
        "--gene-path",
        dest="path_to_interesting_genes",
        help="Path to a newline delimited file of genes of interest.",
        required=True,
    )
    parser.add_argument(
        "--cores",
        dest="cores",
        type=int,
        help="Number of CPUs.",
        default=1,
    )
    parser.add_argument(
        "--racon-path",
        dest="racon_path",
        help="Path to racon installation.",
        default="racon",
    )
    parser.add_argument(
        "--sample-reads",
        dest="sample_reads",
        action="store_true",
        default=False,
        help="Randomly sample to a maximum of 100,000 input reads.",
    )
    parser.add_argument(
        "--quiet",
        dest="quiet",
        action="store_true",
        default=False,
        help="Supress progress updates.",
    )
    parser.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        default=False,
        help="Output Amira debugging files.",
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
    window_length, poly_order = 30, 3  # Example values; need to be chosen based on your data
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

def get_final_filter_threshold(node_coverages,
                            filename):
    # Calculate the frequency of each coverage value
    max_coverage = max(node_coverages)
    coverages = {i: 0 for i in range(max_coverage + 1)}
    for cov in node_coverages:
        if cov in coverages:  # This check ensures we only count values within the range we initialized
            coverages[cov] += 1
    # Sort the coverage values and their frequencies
    vals = sorted(coverages.items())
    x_values = [v[0] for v in vals]
    y_values = [v[1] for v in vals]
    # Apply log transformation to the counts, adding 1 to avoid log(0)
    log_counts = np.log(np.array(y_values) + 1)
    # Smooth the log-transformed histogram counts using a Savitzky-Golay filter
    window_length, poly_order = 20, 3  # Example values; need to be chosen based on your data
    if len(log_counts) < window_length:  # Ensure we have enough data points for the chosen window
        window_length = (
            len(log_counts) // 2 * 2 + 1
        )  # Make the window length the next odd number less than the data length
    smoothed_log_counts = list(savgol_filter(log_counts, window_length, poly_order))
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

    peaks, _ = find_peaks(smoothed_log_counts, prominence=0.2)
    first_peak_index = x_values.index(peaks[0])
    trough_index = smoothed_log_counts.index(min(smoothed_log_counts[:first_peak_index]))
    trough_value = x_values[trough_index]

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

    #distances = random.sample(distances, 10000)
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
    plt.xlim([0, 200])
    # Set plot title and labels
    plt.title("Histogram of distances between genes")
    plt.xlabel("Distance")
    plt.ylabel("Log of absolute frequency")
    # Save the plot as a PNG file
    plt.savefig(filename, dpi=300)  # Adjust dpi for higher resolution if needed
    plt.close()

def build_multiprocessed_graph(annotatedReads, geneMer_size, cores, gene_positions=None):
    batches = [set(list(annotatedReads.keys())[i::cores]) for i in range(cores)]
    if gene_positions is not None:
        sub_graphs = Parallel(n_jobs=cores)(
                        delayed(build_graph)(
                                {k: annotatedReads[k] for k in annotatedReads if k in batch},
                                geneMer_size,
                                {k: gene_positions[k] for k in gene_positions if k in batch}) for batch in batches
                        )
    else:
        sub_graphs = Parallel(n_jobs=cores)(
                        delayed(build_graph)(
                                {k: annotatedReads[k] for k in annotatedReads if k in batch},
                                geneMer_size) for batch in batches
                        )
    merged_graph = merge_graphs(sub_graphs)
    return merged_graph

def write_debug_files(
    annotatedReads: dict[str, list[str]],
    geneMer_size: int,
    genesOfInterest: list[str],
    output_dir: str,
    cores: int
) -> GeneMerGraph:
    sys.stderr.write("\nAmira: building pre-correction gene-mer graph...\n")
    raw_graph = build_multiprocessed_graph(annotatedReads, geneMer_size, cores)
    #raw_graph = GeneMerGraph(annotatedReads, geneMer_size)
    # color nodes in the graph
    for node in raw_graph.all_nodes():
        node.color_node(genesOfInterest)
    import json

    with open(os.path.join(output_dir, "genesAnnotatedOnReads.json"), "w") as o:
        o.write(json.dumps(annotatedReads))
    raw_graph.generate_gml(
        os.path.join(output_dir, "pre_correction_gene_mer_graph"), geneMer_size, 1, 1
    )
    return raw_graph

def plot_unitig_coverages(coverages, filename):
    counts = {}
    for c in range(min(500, max(max(coverages), 100)) + 1):
        if c in coverages:
            counts[c] = coverages.count(c)
        else:
            counts[c] = 0
    x_vals = []
    y_vals = []
    for c in counts:
        x_vals.append(c)
        y_vals.append(counts[c])
    y_vals = np.log10(np.array(y_vals)+1)
    plt.figure(figsize=(10, 6))
    plt.bar(x_vals, y_vals)
    plt.title("Histogram of coverages")
    plt.xlabel("Coverage")
    plt.ylabel("Log of absolute frequency")
    plt.xlim([0, min(500, max(max(coverages), 100))])
    plt.ylim([0, max(y_vals)])
    window_length, poly_order = 50, 3  # Example values; need to be chosen based on your data
    if len(y_vals) < window_length:  # Ensure we have enough data points for the chosen window
        window_length = (
            len(y_vals) // 2 * 2 + 1
        )  # Make the window length the next odd number less than the data length
    smoothed_log_counts = list(savgol_filter(y_vals, window_length, poly_order))
    plt.plot(
        x_vals, smoothed_log_counts, color="red", label="Smoothed counts"
    )
    plt.savefig(filename)
    plt.close()

def main() -> None:
    # get command line options
    args = get_options()
    # make the output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # import the list of genes of interest
    with open(args.path_to_interesting_genes, "r") as i:
        genesOfInterest = i.read().splitlines()
    # import a JSON of genes on reads
    if args.pandoraJSON:
        # output sample information
        sys.stderr.write(
            f"Sample name: {os.path.basename(os.path.dirname(args.pandoraSam))}\nJSON file: {os.path.basename(args.pandoraJSON)}\nFASTQ path: {os.path.basename(args.readfile)}\nSubsample reads: {args.sample_reads}\n"
        )
        annotatedReads, sample_genesOfInterest = process_pandora_json(
            args.pandoraJSON, genesOfInterest
        )
    # load the pandora consensus and convert to a dictionary
    if args.pandoraSam:
        # output sample information
        sys.stderr.write(
            f"Sample name: {os.path.basename(os.path.dirname(args.pandoraSam))}\nSAM file: {os.path.basename(args.pandoraSam)}\nFASTQ path: {os.path.basename(args.readfile)}\nSubsample reads: {args.sample_reads}\n"
        )
        sys.stderr.write("\nAmira: loading Pandora SAM file...\n")
        pandora_consensus = parse_fastq(args.pandoraConsensus)
        try:
            annotatedReads, sample_genesOfInterest, gene_position_dict = convert_pandora_output(
                args.pandoraSam,
                pandora_consensus,
                genesOfInterest,
                args.gene_min_coverage,
                args.lower_gene_length_threshold,
                args.upper_gene_length_threshold,
            )
        except:
            os.remove(args.pandoraConsensus)
            os.remove(args.pandoraSam)
            sys.exit(0)
        # randomly sample 150,000 reads
        if args.sample_reads:
            import random

            annotatedReads = dict(
                random.sample(annotatedReads.items(), min(len(annotatedReads), 100000))
            )
        # plot_log_histogram(distances, os.path.join(args.output_dir, "distance_between_genes.png"))
    print(list(sorted(list(sample_genesOfInterest))))
    # with open(os.path.join(args.output_dir, "gene_gap_sizes.json"), "w") as o:
    #    o.write(json.dumps(distances))
    # plot_log_histogram(proportion_gene_length, os.path.join(args.output_dir, "proportion_genes_covered.png"))
    # write out debug files if specified
    if args.debug:
        write_debug_files(
            annotatedReads, args.geneMer_size, sample_genesOfInterest, args.output_dir, args.cores
        )
    # load the fastq data
    fastq_content = parse_fastq(args.readfile)
    # build the gene-mer graph
    sys.stderr.write("\nAmira: building intitial gene-mer graph...\n")
    #graph = GeneMerGraph(annotatedReads, args.geneMer_size, gene_position_dict)
    graph = build_multiprocessed_graph(annotatedReads, args.geneMer_size, args.cores, gene_position_dict)
    #get_final_filter_threshold(graph.get_all_node_coverages(),
    #                                os.path.join(args.output_dir, f"final_correction_node_coverages.png"))
    # filter junk reads
    graph.filter_graph(2, 1)
    new_annotatedReads, new_gene_position_dict = graph.remove_junk_reads(0.80)
    # dynamically determine the node threshold for filtering
    if not args.node_min_coverage:
        try:
            node_min_coverage = find_trough(
                graph.get_all_node_coverages(),
                os.path.join(args.output_dir, f"node_coverages.png"),
            )
        except:
            node_min_coverage = 5
    else:
        node_min_coverage = args.node_min_coverage
    sys.stderr.write(f"\nAmira: removing low coverage components and nodes with coverage < {node_min_coverage}...\n")
    #graph = GeneMerGraph(new_annotatedReads, args.geneMer_size, new_gene_position_dict)
    graph = build_multiprocessed_graph(new_annotatedReads, args.geneMer_size, args.cores, new_gene_position_dict)
    graph.remove_low_coverage_components(5)
    graph.filter_graph(node_min_coverage, 1)
    new_annotatedReads, new_gene_position_dict = graph.correct_reads(fastq_content)
    #graph = GeneMerGraph(new_annotatedReads, args.geneMer_size, new_gene_position_dict)
    #graph = GeneMerGraph(new_annotatedReads, args.geneMer_size, new_gene_position_dict)
    graph = build_multiprocessed_graph(new_annotatedReads, args.geneMer_size, args.cores, new_gene_position_dict)
    graph.filter_graph(node_min_coverage, 1)
    new_annotatedReads = graph.get_valid_reads_only()
    # parse the original fastq file
    cleaning_iterations = 5
    for this_iteration in range(cleaning_iterations):
        sys.stderr.write(
            f"\nAmira: running graph cleaning iteration {this_iteration+1}/{cleaning_iterations}...\n"
        )
        sys.stderr.write(f"\n\tAmira: removing dead ends...\n")
        #graph = GeneMerGraph(new_annotatedReads, args.geneMer_size, new_gene_position_dict)
        graph = build_multiprocessed_graph(new_annotatedReads, args.geneMer_size, args.cores, new_gene_position_dict)
        graph.remove_short_linear_paths(args.geneMer_size)
        new_annotatedReads, new_gene_position_dict = graph.correct_reads(fastq_content)
        sys.stderr.write(f"\n\tAmira: popping bubbles using {args.cores} CPUs...\n")
        #graph = GeneMerGraph(new_annotatedReads, args.geneMer_size, new_gene_position_dict)
        graph = build_multiprocessed_graph(new_annotatedReads, args.geneMer_size, args.cores, new_gene_position_dict)
        #new_annotatedReads, new_gene_position_dict = graph.correct_low_coverage_paths(fastq_content, sample_genesOfInterest, args.cores)
    # merge paths that are very similar in terms of minimizers
    #graph = GeneMerGraph(new_annotatedReads, args.geneMer_size, new_gene_position_dict)
    graph = build_multiprocessed_graph(new_annotatedReads, args.geneMer_size, args.cores, new_gene_position_dict)
    new_annotatedReads, new_gene_position_dict = graph.correct_low_coverage_paths(fastq_content, sample_genesOfInterest, args.cores, True)
    final_filtering_threshold = get_final_filter_threshold(graph.get_all_node_coverages(),
                                    os.path.join(args.output_dir, f"final_correction_node_coverages.png"))
    # do a final round of filtering
    #graph = GeneMerGraph(new_annotatedReads, args.geneMer_size, new_gene_position_dict)
    graph = build_multiprocessed_graph(new_annotatedReads, args.geneMer_size, args.cores, new_gene_position_dict)
    # decide the threshold for filtering
    graph.filter_graph(final_filtering_threshold, 1)
    new_annotatedReads, new_gene_position_dict = graph.correct_reads(fastq_content)
    #graph = GeneMerGraph(new_annotatedReads, args.geneMer_size, new_gene_position_dict)
    #graph.filter_graph(final_filtering_threshold, 1)
    #new_annotatedReads = graph.get_valid_reads_only()
    # build the corrected gene-mer graph
    sys.stderr.write("\nAmira: building corrected gene-mer graph...\n")
    with open(os.path.join(args.output_dir, "corrected_genesAnnotatedOnReads.json"), "w") as o:
        o.write(json.dumps(new_annotatedReads))
    #graph = GeneMerGraph(new_annotatedReads, args.geneMer_size, new_gene_position_dict)
    graph = build_multiprocessed_graph(new_annotatedReads, args.geneMer_size, args.cores, new_gene_position_dict)
    # remove low coverage components
    graph.remove_low_coverage_components(5)
    # color nodes in the graph if --debug is used
    if args.debug:
        for node in graph.all_nodes():
            node.color_node(sample_genesOfInterest)
    # write out the graph as a GML
    sys.stderr.write("\nAmira: writing gene-mer graph...\n")
    graph.generate_gml(
        os.path.join(args.output_dir, "gene_mer_graph"),
        args.geneMer_size,
        node_min_coverage,
        args.edge_min_coverage,
    )
    # assign reads to AMR genes by path
    sys.stderr.write("\nAmira: clustering reads...\n")
    clusters_of_interest = graph.new_assign_reads_to_genes(sample_genesOfInterest, fastq_content)
    # write out the fastq files
    if not os.path.exists(os.path.join(args.output_dir, "AMR_allele_fastqs")):
        os.mkdir(os.path.join(args.output_dir, "AMR_allele_fastqs"))
    # subset the fastq data based on the cluster assignments
    files_to_assemble = []
    sys.stderr.write("\nAmira: writing fastqs...\n")
    for allele in tqdm(clusters_of_interest):
        read_subset = {}
        for r in clusters_of_interest[allele]:
            underscore_split = r.split("_")
            fastq_data = fastq_content[underscore_split[0]].copy()
#            print(underscore_split, len(fastq_data["sequence"]))
            fastq_data["sequence"] = fastq_data["sequence"][max([0, int(underscore_split[1]) - 100]): min([len(fastq_data["sequence"])-1, int(underscore_split[2]) + 101])]
            fastq_data["quality"] = fastq_data["quality"][max([0, int(underscore_split[1]) - 100]): min([len(fastq_data["quality"])-1, int(underscore_split[2]) + 101])]
            read_subset[underscore_split[0]] = fastq_data
            assert read_subset[underscore_split[0]]["sequence"] != ""
        write_fastq(
            os.path.join(args.output_dir, "AMR_allele_fastqs", allele + ".fastq.gz"), read_subset
        )
        files_to_assemble.append(
            os.path.join(args.output_dir, "AMR_allele_fastqs", allele + ".fastq.gz")
        )
    # run racon to polish the pandora consensus
    sys.stderr.write("\nAmira: obtaining nucleotide sequences...\n")
    genes_to_remove = graph.polish_pandora_consensus(
        files_to_assemble,
        args.racon_path,
        pandora_consensus,
        1,
        os.path.join(args.output_dir, "AMR_allele_fastqs"),
    )
    # remove genes that do not have sufficient mapping coverage
    for g in genes_to_remove:
        if g is not None:
            sys.stderr.write(
                f"\nAmira: allele {g[0]} removed due to insufficient coverage ({g[1]}).\n"
            )
            del clusters_of_interest[g[0]]
            #shutil.rmtree(os.path.join(args.output_dir, "AMR_allele_fastqs", g[0]))
            #os.remove(os.path.join(args.output_dir, "AMR_allele_fastqs", g[0] + ".fastq.gz"))
    # write out the clustered reads
    for allele in clusters_of_interest:
        new_reads = set()
        for r in clusters_of_interest[allele]:
            new_reads.add(r.split("_")[0])
        clusters_of_interest[allele] = list(new_reads)
    with open(os.path.join(args.output_dir, "reads_per_amr_gene.json"), "w") as o:
        o.write(json.dumps(clusters_of_interest))
    sys.exit(0)


if __name__ == "__main__":
    main()
