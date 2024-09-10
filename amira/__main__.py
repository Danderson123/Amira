import argparse
import gzip
import json
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
from scipy.signal import find_peaks, savgol_filter
from tqdm import tqdm

from amira.construct_graph import GeneMerGraph, build_graph, merge_graphs
from amira.pre_process import convert_pandora_output, process_pandora_json

matplotlib.use("Agg")


def get_options() -> argparse.Namespace:
    """define args from the command line"""
    parser = argparse.ArgumentParser(description="Build a prototype gene de Bruijn graph.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pandoraSam", dest="pandoraSam", help="Pandora map SAM file path.")
    group.add_argument("--pandoraJSON", dest="pandoraJSON", help="Pandora map JSON file path.")
    group.add_argument(
        "--gene-positions", dest="gene_positions", help="Gene position JSON file path."
    )
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
        default=3,
        help="Minimum threshold for gene-mer coverage.",
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
        help="Path to a multi_FASTA file of the AMR gene alleles of interest.",
        required=True,
    )
    parser.add_argument(
        "--phenotypes",
        dest="phenotypes",
        help="Path to a JSON of phenotypes for each AMR allele.",
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


def plot_read_length_distribution(annotatedReads, output_dir):
    read_lengths = []
    for read in annotatedReads:
        read_lengths.append(len(annotatedReads[read]))
    plt.figure(figsize=(10, 6))
    plt.hist(read_lengths, bins=50, edgecolor="black")
    plt.title("Number of genes per read")
    plt.xlabel("Number of genes")
    plt.ylabel("Absolute frequency")
    plt.savefig(os.path.join(output_dir, "read_lengths.png"), dpi=600)
    plt.close()


def write_debug_files(
    annotatedReads: dict[str, list[str]],
    geneMer_size: int,
    genesOfInterest: list[str],
    output_dir: str,
    cores: int,
) -> GeneMerGraph:
    sys.stderr.write("\nAmira: building pre-correction gene-mer graph\n")
    raw_graph = build_multiprocessed_graph(annotatedReads, geneMer_size, cores)
    # color nodes in the graph
    for node in raw_graph.all_nodes():
        node.color_node(genesOfInterest)
    raw_graph.generate_gml(
        os.path.join(output_dir, "pre_correction_gene_mer_graph"), geneMer_size, 1, 1
    )
    return raw_graph


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


def estimate_copy_number(amira_allele, copy_numbers_per_component, additional_copy_numbers):
    copy_numbers = {}
    for component in copy_numbers_per_component:
        for gene in copy_numbers_per_component[component]:
            for allele in copy_numbers_per_component[component][gene]:
                copy_numbers[allele] = round(copy_numbers_per_component[component][gene][allele], 2)
    if amira_allele in copy_numbers:
        return copy_numbers[amira_allele]
    return additional_copy_numbers[amira_allele]


def write_allele_fastq(reads_for_allele, fastq_content, output_dir, allele_name):
    read_subset = {}
    for r in reads_for_allele:
        underscore_split = r.split("_")
        fastq_data = fastq_content[underscore_split[0]].copy()
        fastq_data["sequence"] = fastq_data["sequence"][
            max([0, int(underscore_split[1]) - 250]) : min(
                [len(fastq_data["sequence"]) - 1, int(underscore_split[2]) + 250]
            )
        ]
        fastq_data["quality"] = fastq_data["quality"][
            max([0, int(underscore_split[1]) - 250]) : min(
                [len(fastq_data["quality"]) - 1, int(underscore_split[2]) + 250]
            )
        ]
        if fastq_data["sequence"] != "":
            read_subset[underscore_split[0]] = fastq_data
    if not os.path.exists(os.path.join(output_dir, "AMR_allele_fastqs", allele_name)):
        os.mkdir(os.path.join(output_dir, "AMR_allele_fastqs", allele_name))
    write_fastq(
        os.path.join(output_dir, "AMR_allele_fastqs", allele_name, allele_name + ".fastq.gz"),
        read_subset,
    )
    return os.path.join(output_dir, "AMR_allele_fastqs", allele_name, allele_name + ".fastq.gz")


def parse_fastq_lines(fh):
    # Initialize a counter to keep track of the current line number
    line_number = 0
    # Iterate over the lines in the file
    for line in fh:
        # Increment the line number
        line_number += 1
        # If the line number is divisible by 4, it's a sequence identifier line
        if line_number % 4 == 1:
            # Extract the identifier from the line
            identifier = line.split(" ")[0][1:]
        # If the line number is divisible by 4, it's a sequence line
        elif line_number % 4 == 2:
            sequence = line.strip()
        elif line_number % 4 == 0:
            # Yield the identifier, sequence and quality
            yield identifier, sequence, line.strip()


def parse_fastq(fastq_file):
    # Initialize an empty dictionary to store the results
    results = {}
    # Open the fastq file
    if ".gz" in fastq_file:
        try:
            with gzip.open(fastq_file, "rt") as fh:
                # Iterate over the lines in the file
                for identifier, sequence, quality in parse_fastq_lines(fh):
                    # Add the identifier and sequence to the results dictionary
                    results[identifier.replace("\n", "")] = {
                        "sequence": sequence,
                        "quality": quality,
                    }
            return results
        except OSError:
            pass
    with open(fastq_file, "r") as fh:
        # Iterate over the lines in the file
        for identifier, sequence, quality in parse_fastq_lines(fh):
            # Add the identifier and sequence to the results dictionary
            results[identifier.replace("\n", "")] = {"sequence": sequence, "quality": quality}
    # Return the dictionary of results
    return results


def write_fastq(fastq_file, data):
    # Open the fastq file
    with gzip.open(fastq_file, "wt") as fh:
        # Iterate over the data
        for identifier, value in data.items():
            # Write the identifier line
            fh.write(f"@{identifier}\n")
            # Write the sequence line
            fh.write(f'{value["sequence"]}\n')
            # Write the placeholder quality lines
            fh.write("+\n")
            fh.write(f'{value["quality"]}\n')


def main() -> None:
    # get command line options
    args = get_options()
    # make the output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # import the list of genes of interest
    with open(args.path_to_interesting_genes, "r") as i:
        reference_content = i.read().split(">")[1:]
    genesOfInterest = set()
    reference_alleles = {}
    for allele in reference_content:
        newline_split = allele.split("\n")
        assert (
            newline_split[0].count(";") == 1
        ), "Reference FASTA headers can only contain 1 semicolon"
        gene_name, allele_name = newline_split[0].split(";")
        genesOfInterest.add(gene_name)
        sequence = "".join(newline_split[1:])
        if gene_name not in reference_alleles:
            reference_alleles[gene_name] = {}
        reference_alleles[gene_name][allele_name] = sequence
    # import a JSON of genes on reads
    if args.pandoraJSON:
        # output sample information
        sys.stderr.write(
            f"Sample name: {os.path.basename(os.path.dirname(args.pandoraSam))}\n"
            f"JSON file: {os.path.basename(args.pandoraJSON)}\n"
            f"FASTQ path: {os.path.basename(args.readfile)}\n"
            f"Subsample reads: {args.sample_reads}\n"
        )
        annotatedReads, sample_genesOfInterest = process_pandora_json(
            args.pandoraJSON, genesOfInterest
        )
    # load the pandora consensus and convert to a dictionary
    if args.pandoraSam:
        # output sample information
        sys.stderr.write(
            f"Sample name: {os.path.basename(os.path.dirname(args.pandoraSam))}\n"
            f"SAM file: {os.path.basename(args.pandoraSam)}\n"
            f"FASTQ path: {os.path.basename(args.readfile)}\n"
            f"Subsample reads: {args.sample_reads}\n"
        )
        sys.stderr.write("\nAmira: loading Pandora SAM file\n")
        pandora_consensus = parse_fastq(args.pandoraConsensus)
        annotatedReads, sample_genesOfInterest, gene_position_dict = convert_pandora_output(
            args.pandoraSam,
            pandora_consensus,
            genesOfInterest,
            args.gene_min_coverage,
            args.lower_gene_length_threshold,
            args.upper_gene_length_threshold,
            args.readfile,
            args.cores,
            args.output_dir,
        )
        with open(
            os.path.join(args.output_dir, "gene_positions_with_gene_filtering.json"), "w"
        ) as o:
            o.write(json.dumps(gene_position_dict))
        with open(os.path.join(args.output_dir, "gene_calls_with_gene_filtering.json"), "w") as o:
            o.write(json.dumps(annotatedReads))
        # randomly sample 150,000 reads
        if args.sample_reads:
            import random

            annotatedReads = dict(
                random.sample(annotatedReads.items(), min(len(annotatedReads), 100000))
            )
    print(list(sorted(list(sample_genesOfInterest))))
    # terminate if no AMR genes were found
    if len(sample_genesOfInterest) == 0:
        # write an empty dataframe
        results = "Gene name\tSequence name\tClosest reference\tReference length\t"
        results += (
            "Identity (%)\tCoverage (%)\tAmira allele\tNumber of reads\tApproximate copy number\n"
        )
        with open(os.path.join(args.output_dir, "amira_results.tsv"), "w") as o:
            o.write(results)
        # exit
        sys.exit(0)
    # load the fastq data
    fastq_content = parse_fastq(args.readfile)
    # write out debug files if specified
    if args.debug:
        plot_read_length_distribution(annotatedReads, args.output_dir)
        write_debug_files(
            annotatedReads, args.geneMer_size, sample_genesOfInterest, args.output_dir, args.cores
        )
    # build the gene-mer graph
    sys.stderr.write("\nAmira: building intitial gene-mer graph\n")
    # graph = GeneMerGraph(annotatedReads, geneMer_size, gene_position_dict)
    graph = build_multiprocessed_graph(
        annotatedReads, args.geneMer_size, args.cores, gene_position_dict
    )
    try:
        min_path_coverage = plot_node_coverages(
            graph.get_all_node_coverages(),
            os.path.join(args.output_dir, "initial_node_coverages.png"),
        )
    except (ValueError, IndexError):
        min_path_coverage = 10
    # collect the reads that have fewer than k genes
    short_reads = graph.get_short_read_annotations()
    short_read_gene_positions = graph.get_short_read_gene_positions()
    # filter junk reads
    graph.filter_graph(2, 1)
    new_annotatedReads, new_gene_position_dict, rejected_reads, rejected_read_positions = graph.remove_junk_reads(0.80)
    node_min_coverage = args.node_min_coverage
    sys.stderr.write(
        f"\nAmira: removing low coverage components and nodes with coverage < {node_min_coverage}\n"
    )
    graph = build_multiprocessed_graph(
        new_annotatedReads, args.geneMer_size, args.cores, new_gene_position_dict
    )
    # collect the reads that have fewer than k genes
    short_reads.update(graph.get_short_read_annotations())
    short_read_gene_positions.update(graph.get_short_read_gene_positions())
    graph.remove_low_coverage_components(5)
    graph.filter_graph(node_min_coverage, 1)
    new_annotatedReads, new_gene_position_dict = graph.correct_reads(fastq_content)
    graph = build_multiprocessed_graph(
        new_annotatedReads, args.geneMer_size, args.cores, new_gene_position_dict
    )
    # collect the reads that have fewer than k genes
    short_reads.update(graph.get_short_read_annotations())
    short_read_gene_positions.update(graph.get_short_read_gene_positions())
    graph.filter_graph(node_min_coverage, 1)
    new_annotatedReads = graph.get_valid_reads_only()
    # choose a value for k
    sys.stderr.write("\nAmira: selecting a gene-mer size (k)\n")
    geneMer_size = 3
    for k in range(3, 16, 2):
        # Build the graph with the current k value
        graph = build_multiprocessed_graph(
            new_annotatedReads.copy(), k, args.cores, new_gene_position_dict.copy()
        )

        def is_component_valid(component):
            amr_nodes = {
                n.__hash__() for g in sample_genesOfInterest for n in graph.get_nodes_containing(g)
            }
            nodes_in_component = [n.__hash__() for n in graph.get_nodes_in_component(component)]
            reads = graph.collect_reads_in_path([n for n in nodes_in_component if n in amr_nodes])
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
    sys.stderr.write(f"\nAmira: selected k={geneMer_size}\n")
    # parse the original fastq file
    cleaning_iterations = 10
    prev_nodes = 0
    for this_iteration in range(cleaning_iterations):
        sys.stderr.write(
            f"\nAmira: running graph cleaning iteration {this_iteration+1}/{cleaning_iterations}\n"
        )
        sys.stderr.write("\n\tAmira: removing dead ends\n")
        graph = build_multiprocessed_graph(
            new_annotatedReads, geneMer_size, args.cores, new_gene_position_dict
        )
        # check if the current number of nodes is equal to the previous number of nodes
        if len(graph.get_nodes()) == prev_nodes:
            sys.stderr.write(f"\n\tAmira: terminating cleaning at iteration {this_iteration+1}\n")
            break
        prev_nodes = len(graph.get_nodes())
        # collect the reads that have fewer than k genes
        short_reads.update(graph.get_short_read_annotations())
        short_read_gene_positions.update(graph.get_short_read_gene_positions())
        graph.remove_short_linear_paths(geneMer_size)
        new_annotatedReads, new_gene_position_dict = graph.correct_reads(fastq_content)
        sys.stderr.write(f"\n\tAmira: popping bubbles using {args.cores} CPUs\n")
        graph = build_multiprocessed_graph(
            new_annotatedReads, geneMer_size, args.cores, new_gene_position_dict
        )
        # collect the reads that have fewer than k genes
        short_reads.update(graph.get_short_read_annotations())
        short_read_gene_positions.update(graph.get_short_read_gene_positions())
        new_annotatedReads, new_gene_position_dict, path_coverages, min_path_coverage = (
            graph.correct_low_coverage_paths(
                fastq_content, sample_genesOfInterest, args.cores, min_path_coverage, True
            )
        )
    # build the corrected gene-mer graph
    sys.stderr.write("\nAmira: building corrected gene-mer graph\n")
    with open(os.path.join(args.output_dir, "corrected_gene_calls_after_filtering.json"), "w") as o:
        o.write(json.dumps(new_annotatedReads))
    with open(
        os.path.join(args.output_dir, "corrected_gene_positions_after_filtering.json"), "w"
    ) as o:
        o.write(json.dumps(new_gene_position_dict))
    graph = build_multiprocessed_graph(
        new_annotatedReads, geneMer_size, args.cores, new_gene_position_dict
    )
    # collect the reads that have fewer than k genes
    short_reads.update(graph.get_short_read_annotations())
    short_read_gene_positions.update(graph.get_short_read_gene_positions())
    # remove low coverage components
    graph.remove_low_coverage_components(5)
    # color nodes in the graph if --debug is used
    if args.debug:
        for node in graph.all_nodes():
            node.color_node(sample_genesOfInterest)
    # write out the graph as a GML
    sys.stderr.write("\nAmira: writing gene-mer graph\n")
    graph.generate_gml(
        os.path.join(args.output_dir, "gene_mer_graph"),
        geneMer_size,
        node_min_coverage,
        1,
    )
    # assign reads to AMR genes by path
    sys.stderr.write("\nAmira: clustering reads\n")
    clusters_of_interest, cluster_copy_numbers, allele_counts = graph.assign_reads_to_genes(
        sample_genesOfInterest, fastq_content
    )
    # get the unique genes we have found
    found_genes_of_interest = set()
    # iterate through the components
    for component_id in clusters_of_interest:
        for gene in clusters_of_interest[component_id]:
            found_genes_of_interest.add(gene)
    # add AMR alleles that have come from the reads shorter than k
    clusters_to_add = {}
    cluster_copy_numbers_to_add = {}
    for read_id in short_reads:
        for g in range(len(short_reads[read_id])):
            strandless_gene = short_reads[read_id][g][1:]
            if (
                strandless_gene in sample_genesOfInterest
                and strandless_gene not in found_genes_of_interest
            ):
                if f"{strandless_gene}_1" not in clusters_to_add:
                    clusters_to_add[f"{strandless_gene}_1"] = []
                gene_start = short_read_gene_positions[read_id][g][0]
                gene_end = short_read_gene_positions[read_id][g][1]
                clusters_to_add[f"{strandless_gene}_1"].append(f"{read_id}_{gene_start}_{gene_end}")
    for amira_allele in clusters_to_add:
        cluster_copy_numbers_to_add[amira_allele] = max(
            1.0, len(clusters_to_add[amira_allele]) / graph.get_mean_node_coverage()
        )
    # write out the fastq files
    if not os.path.exists(os.path.join(args.output_dir, "AMR_allele_fastqs")):
        os.mkdir(os.path.join(args.output_dir, "AMR_allele_fastqs"))
    # subset the fastq data based on the cluster assignments
    files_to_assemble = []
    sys.stderr.write("\nAmira: writing fastqs\n")
    supplemented_clusters_of_interest = {}
    for component in tqdm(clusters_of_interest):
        for gene in clusters_of_interest[component]:
            for allele in clusters_of_interest[component][gene]:
                if len(clusters_of_interest[component][gene][allele]) > 4:
                    files_to_assemble.append(
                        write_allele_fastq(
                            clusters_of_interest[component][gene][allele],
                            fastq_content,
                            args.output_dir,
                            allele,
                        )
                    )
                    supplemented_clusters_of_interest[allele] = clusters_of_interest[component][
                        gene
                    ][allele]
    # add the genes from the short reads
    for allele in clusters_to_add:
        if len(clusters_to_add[allele]) > 4:
            files_to_assemble.append(
                write_allele_fastq(clusters_to_add[allele], fastq_content, args.output_dir, allele)
            )
            supplemented_clusters_of_interest[allele] = clusters_to_add[allele]
        else:
            sys.stderr.write(
                f"Amira: allele {allele} filtered due to an insufficient number of reads.\n"
            )
    # run racon to polish the pandora consensus
    sys.stderr.write("\nAmira: obtaining nucleotide sequences\n")
    result_df = graph.get_alleles(
        files_to_assemble,
        args.racon_path,
        args.cores,
        os.path.join(args.output_dir, "AMR_allele_fastqs"),
        reference_alleles,
        pandora_consensus,
        args.phenotypes,
        args.debug,
    )

    # get the copy number estimates of each allele
    result_df["Approximate copy number"] = result_df.apply(
        lambda row: estimate_copy_number(
            row["Amira allele"], cluster_copy_numbers, cluster_copy_numbers_to_add
        ),
        axis=1,
    )
    # remove genes that do not have sufficient mapping coverage
    alleles_to_delete = []
    for index, row in result_df.iterrows():
        if row["Identity (%)"] < 90:
            message = f"\nAmira: allele {row['Amira allele']} removed "
            message += f"due to insufficient similarity ({row['Identity (%)']}).\n"
            sys.stderr.write(message)
            alleles_to_delete.append(row["Amira allele"])
        else:
            if row["Coverage (%)"] < 85:
                message = f"\nAmira: allele {row['Amira allele']} removed "
                message += f"due to insufficient coverage ({row['Coverage (%)']}).\n"
                sys.stderr.write(message)
                alleles_to_delete.append(row["Amira allele"])
    # remove genes as necessary
    for amira_allele in alleles_to_delete:
        del supplemented_clusters_of_interest[amira_allele]
        result_df = result_df[result_df["Amira allele"] != amira_allele]
    # write out the clustered reads
    final_clusters_of_interest = {}
    for allele in supplemented_clusters_of_interest:
        # get the gene name with the allele name appended
        if os.path.exists(
            os.path.join(args.output_dir, "AMR_allele_fastqs", allele, "06.final_sequence.fasta")
        ):
            with open(
                os.path.join(
                    args.output_dir, "AMR_allele_fastqs", allele, "06.final_sequence.fasta"
                )
            ) as i:
                reference_allele_name = i.read().split(" ")[0].replace(">", "")
        else:
            with open(
                os.path.join(
                    args.output_dir, "AMR_allele_fastqs", allele, "03.sequence_to_polish.fasta"
                )
            ) as i:
                reference_allele_name = i.read().split(" ")[0].replace(">", "")
        if "\n" in reference_allele_name:
            reference_allele_name = reference_allele_name.split("\n")[0]
        new_name = f"{allele};{reference_allele_name}"
        new_reads = set()
        for r in supplemented_clusters_of_interest[allele]:
            new_reads.add(r.split("_")[0])
        final_clusters_of_interest[new_name] = list(new_reads)
    with open(os.path.join(args.output_dir, "reads_per_amr_gene.json"), "w") as o:
        o.write(json.dumps(final_clusters_of_interest))
    # write the result tsv
    result_df.to_csv(os.path.join(args.output_dir, "amira_results.tsv"), sep="\t", index=False)
    # make the result table
    sys.exit(0)


if __name__ == "__main__":
    main()
