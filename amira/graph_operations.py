import os
import statistics
import sys

import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
from scipy.signal import find_peaks, savgol_filter

from amira.construct_graph import GeneMerGraph


def build_graph(read_dict, kmer_size, gene_positions=None):
    graph = GeneMerGraph(read_dict, kmer_size, gene_positions)
    return graph


def merge_nodes(sub_graphs, fastq_data=None):
    reference_graph = sub_graphs[0]
    # iterate through the subgraphs
    for graph in sub_graphs[1:]:
        # iterate through the reads in the subgraph
        for read_id in graph.get_readNodes():
            # get the nodes on each read
            node_hashes_on_read = graph.get_readNodes()[read_id]
            # get the direction of each node on each read
            node_directions_on_read = graph.get_readNodeDirections()[read_id]
            # get the positions of each node on each read
            node_positions_on_read = graph.get_readNodePositions()[read_id]
            # iterate through the nodes
            for i in range(len(node_hashes_on_read)):
                # get the node object
                node_in_subgraph = graph.get_node_by_hash(node_hashes_on_read[i])
                # add the node to the reference graph
                node_in_reference = reference_graph.add_node(
                    node_in_subgraph.get_geneMer(), node_in_subgraph.get_reads()
                )
                # add to the minhash
                if fastq_data is not None:
                    if node_in_reference.get_minhash() is None:
                        node_in_reference.set_minhash(node_in_subgraph.get_minhash())
                    node_in_reference.get_minhash().add_many(node_in_subgraph.get_minhash())
                # increment the node coverage
                node_in_reference.increment_node_coverage()
                # add the node to the read in the reference graph
                reference_graph.add_node_to_read(
                    node_in_reference,
                    read_id,
                    node_directions_on_read[i],
                    node_positions_on_read[i],
                )
    return reference_graph


def merge_edges(sub_graphs, reference_graph):
    for graph in sub_graphs[1:]:
        for edge_hash in graph.get_edges():
            if edge_hash not in reference_graph.get_edges():
                # get the edge object in the subgraph
                edge_in_subgraph = graph.get_edge_by_hash(edge_hash)
                # change the source and target node to those in the reference graph
                reference_source_node = edge_in_subgraph.set_sourceNode(
                    reference_graph.get_node_by_hash(edge_in_subgraph.get_sourceNode().__hash__())
                )
                edge_in_subgraph.set_targetNode(
                    reference_graph.get_node_by_hash(edge_in_subgraph.get_targetNode().__hash__())
                )
                # add the edge to the node
                reference_graph.add_edge_to_node(reference_source_node, edge_in_subgraph)
                # add the edge object to the graph
                reference_graph.get_edges()[edge_hash] = edge_in_subgraph
            else:
                # get the reference edge
                reference_edge = reference_graph.get_edge_by_hash(edge_hash)
                # extend the coverage
                reference_edge.extend_edge_coverage(reference_edge.get_edge_coverage())


def merge_reads(sub_graphs, reference_graph):
    for graph in sub_graphs[1:]:
        for read in graph.get_reads():
            reference_graph.get_reads()[read] = graph.get_reads()[read]
            if reference_graph.get_gene_positions() is not None:
                reference_graph.get_gene_positions()[read] = graph.get_gene_positions()[read]
        for read in graph.get_short_read_annotations():
            reference_graph.get_short_read_annotations()[read] = graph.get_short_read_annotations()[
                read
            ]
            if graph.get_gene_positions() is not None:
                reference_graph.get_short_read_gene_positions()[
                    read
                ] = graph.get_short_read_gene_positions()[read]


def merge_graphs(sub_graphs):
    # merge the nodes
    reference_graph = merge_nodes(sub_graphs)
    # merge the edges
    merge_edges(sub_graphs, reference_graph)
    # merge the reads
    merge_reads(sub_graphs, reference_graph)
    reference_graph.assign_component_ids()
    return reference_graph


def build_multiprocessed_graph(annotatedReads, geneMer_size, cores, gene_positions=None):
    batches = [set(list(annotatedReads.keys())[i::cores]) for i in range(cores)]
    if gene_positions is not None:
        sub_graphs = Parallel(n_jobs=cores)(
            delayed(build_graph)(
                {k: annotatedReads[k] for k in annotatedReads if k in batch},
                geneMer_size,
                {k: gene_positions[k] for k in gene_positions if k in batch},
            )
            for batch in batches
        )
    else:
        sub_graphs = Parallel(n_jobs=cores)(
            delayed(build_graph)(
                {k: annotatedReads[k] for k in annotatedReads if k in batch}, geneMer_size
            )
            for batch in batches
        )
    merged_graph = merge_graphs(sub_graphs)
    return merged_graph


def iterative_bubble_popping(
    new_annotatedReads,
    new_gene_position_dict,
    cleaning_iterations,
    geneMer_size,
    cores,
    short_reads,
    short_read_gene_positions,
    fastq_content,
    output_dir,
    node_min_coverage,
    sample_genesOfInterest,
    min_path_coverage,
):
    prev_nodes = 0
    components_to_skip = set()
    for this_iteration in range(cleaning_iterations):
        sys.stderr.write(
            f"\nAmira: running graph cleaning iteration {this_iteration+1}/{cleaning_iterations}\n"
        )
        graph = build_multiprocessed_graph(
            new_annotatedReads, geneMer_size, cores, new_gene_position_dict
        )
        graph.filter_graph(node_min_coverage, 1)
        new_annotatedReads, new_gene_position_dict = graph.correct_reads(fastq_content)
        graph = build_multiprocessed_graph(
            new_annotatedReads, geneMer_size, cores, new_gene_position_dict
        )
        # check if the current number of nodes is equal to the previous number of nodes
        if len(graph.get_nodes()) == prev_nodes:
            sys.stderr.write(f"\n\tAmira: terminating cleaning at iteration {this_iteration+1}\n")
            break
        prev_nodes = len(graph.get_nodes())
        sys.stderr.write("\n\tAmira: removing dead ends\n")
        # collect the reads that have fewer than k genes
        short_reads.update(graph.get_short_read_annotations())
        short_read_gene_positions.update(graph.get_short_read_gene_positions())
        graph.remove_short_linear_paths(geneMer_size)
        new_annotatedReads, new_gene_position_dict = graph.correct_reads(fastq_content)
        sys.stderr.write("\n\tAmira: popping bubbles using 1 CPU\n")
        graph = build_multiprocessed_graph(
            new_annotatedReads, geneMer_size, cores, new_gene_position_dict
        )
        # collect the reads that have fewer than k genes
        short_reads.update(graph.get_short_read_annotations())
        short_read_gene_positions.update(graph.get_short_read_gene_positions())
        graph.generate_gml(
            os.path.join(output_dir, "intermediate_graph"),
            geneMer_size,
            node_min_coverage,
            1,
        )
        new_annotatedReads, new_gene_position_dict, path_coverages, min_path_coverage = (
            graph.correct_low_coverage_paths(
                fastq_content,
                sample_genesOfInterest,
                cores,
                min_path_coverage,
                components_to_skip,
                True,
            )
        )
    return new_annotatedReads, new_gene_position_dict


def plot_node_coverages(unitig_coverages, filename):
    # Calculate the frequency of each coverage value with bins of width 5
    max_coverage = max(unitig_coverages)
    # Create bins with a step size of 5
    bins = np.arange(0, max_coverage + 5, 5)  # Bins of width 5
    hist, bin_edges = np.histogram(unitig_coverages, bins=bins)
    # Midpoints of bins for plotting
    x_values = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    # Apply log transformation to the counts, adding 1 to avoid log(0)
    log_counts = np.log(hist + 1)
    # Smooth the log-transformed histogram counts using a Savitzky-Golay filter
    window_length, poly_order = 31, 5  # Window length must be odd
    if len(log_counts) < window_length:
        window_length = max(5, len(log_counts) // 2 * 2 - 1)  # Smallest odd number >= 3
    smoothed_log_counts = savgol_filter(log_counts, window_length, poly_order)

    # Plot histogram
    plt.figure(figsize=(10, 6))
    plt.bar(
        x_values,
        log_counts,
        width=5,
        label="Counts",
        color="white",
        edgecolor="black",
        align="center",
    )
    plt.plot(x_values, smoothed_log_counts, color="red", label="Smoothed counts")
    plt.title("Histogram of mean unitig coverages with Smoothed Curve")
    plt.xlabel("Unitig Coverage")
    plt.ylabel("Log of absolute frequency")
    plt.xlim([0, max(x_values) + 5])
    plt.legend()
    plt.savefig(filename)
    plt.close()

    # Identify peaks and troughs
    peaks, _ = find_peaks(
        [min(smoothed_log_counts)] + list(smoothed_log_counts),
        [min(smoothed_log_counts)],
        prominence=0.1,
    )
    peaks = peaks - 1
    first_peak_index = np.where(x_values == x_values[peaks[0]])[0][0]
    second_peak_index = np.where(x_values == x_values[peaks[1]])[0][0]
    trough_index = (
        np.argmin(smoothed_log_counts[first_peak_index : second_peak_index + 1]) + first_peak_index
    )
    trough_value = x_values[trough_index]

    # Plot the histogram and the trough
    plt.figure(figsize=(10, 6))
    plt.bar(
        x_values,
        log_counts,
        width=5,
        label="Counts",
        color="white",
        edgecolor="black",
        align="center",
    )
    plt.plot(x_values, smoothed_log_counts, color="red", label="Smoothed counts")
    plt.axvline(x=trough_value, color="r", linestyle="--", label=f"Trough at x={trough_value:.2f}")
    plt.title("Histogram of node coverages with Smoothed Curve")
    plt.xlabel("Node Coverage")
    plt.ylabel("Log of absolute frequency")
    plt.xlim([0, max(x_values) + 5])
    plt.legend()
    plt.savefig(filename)
    plt.close()

    return trough_value


def choose_kmer_size(
    overall_mean_node_coverage,
    new_annotatedReads,
    cores,
    new_gene_position_dict,
    sample_genesOfInterest,
):
    geneMer_size = 3
    if overall_mean_node_coverage >= 20:
        for k in range(3, 16, 2):
            # Build the graph with the current k value
            graph = build_multiprocessed_graph(
                new_annotatedReads.copy(), k, cores, new_gene_position_dict.copy()
            )

            def is_component_valid(component):
                amr_nodes = {
                    n.__hash__()
                    for g in sample_genesOfInterest
                    for n in graph.get_nodes_containing(g)
                }
                nodes_in_component = [n.__hash__() for n in graph.get_nodes_in_component(component)]
                reads = graph.collect_reads_in_path(
                    [n for n in nodes_in_component if n in amr_nodes]
                )
                lengths = [len(graph.get_reads()[r]) for r in reads]
                if len(lengths) != 0:
                    return (
                        len([length for length in lengths if length >= (2 * k - 1)]) / len(lengths)
                        >= 0.8
                    )
                else:
                    return True

            if all(is_component_valid(c) for c in graph.components()):
                geneMer_size = k
            else:
                break
    return geneMer_size


def get_overall_mean_node_coverages(graph):
    overall_mean_node_coverages = {}
    for k in range(3, 16, 2):
        coverages = []
        for node in graph.all_nodes():
            cov = 0
            for read in node.get_reads():
                if len(graph.get_reads()[read]) >= k:
                    cov += 1
            coverages.append(cov)
        if not len(coverages) == 0:
            overall_mean_node_coverages[k] = statistics.mean(coverages)
        else:
            overall_mean_node_coverages[k] = 0
    return overall_mean_node_coverages
