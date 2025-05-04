import argparse
import os
import random
import sys
import time

import matplotlib

from amira.__init__ import __version__
from amira.construct_graph import GeneMerGraph
from amira.graph_utils import (
    build_multiprocessed_graph,
    choose_kmer_size,
    get_overall_mean_node_coverages,
    iterative_bubble_popping,
    plot_node_coverages,
)
from amira.pre_processing import (
    convert_pandora_output,
    estimate_mean_core_gene_counts,
    load_species_specific_files,
    process_pandora_json,
    process_reference_alleles,
    run_pandora_map,
    subsample_reads_and_estimate_read_depth,
)
from amira.read_utils import (
    parse_fastq,
    plot_read_length_distribution,
    write_modified_fastq,
)
from amira.result_utils import (
    estimate_copy_numbers,
    filter_results,
    genotype_promoters,
    get_alleles,
    output_component_fastqs,
    process_reads,
    supplement_result_df,
    write_empty_result,
    write_fastqs_for_genes,
    write_fastqs_for_genes_with_short_reads,
    write_pandora_gene_calls,
    write_reads_per_AMR_gene,
)

matplotlib.use("Agg")


def get_options() -> argparse.Namespace:
    """define args from the command line"""
    parser = argparse.ArgumentParser(
        description="Identify acquired AMR genes from bacterial long read sequences."
    )
    parser.add_argument("--pandoraSam", dest="pandoraSam", help=argparse.SUPPRESS, default=None)
    parser.add_argument("--pandoraJSON", dest="pandoraJSON", help=argparse.SUPPRESS, default=None)
    parser.add_argument("--gene-positions", help=argparse.SUPPRESS, default=None)
    parser.add_argument(
        "--pandoraConsensus",
        dest="pandoraConsensus",
        help=argparse.SUPPRESS,
        required=False,
    )
    parser.add_argument(
        "--reads", dest="reads", help="path to FASTQ file of long reads.", required=True
    )
    parser.add_argument(
        "--species",
        dest="species",
        choices=[
            "Escherichia_coli",
            "Klebsiella_pneumoniae",
            "Enterococcus_faecium",
            "Streptococcus_pneumoniae",
            "Staphylococcus_aureus",
        ],
        help="The species you want to run Amira on.",
        required=True,
    )
    parser.add_argument(
        "--panRG-path",
        dest="panRG_path",
        help="Path to pandora panRG ending .panidx.zip.",
        required=True,
    )
    parser.add_argument(
        "--output",
        dest="output_dir",
        type=str,
        default="amira_output",
        help="Directory for Amira outputs (default=amira_output).",
    )
    parser.add_argument(
        "-n",
        dest="node_min_coverage",
        type=int,
        default=3,
        help="Minimum threshold for gene-mer coverage in the graph (default=3).",
    )
    parser.add_argument(
        "-g",
        dest="gene_min_coverage",
        type=float,
        default=0.2,
        help="Minimum relative threshold to remove all instances of a gene (default=0.2).",
    )
    parser.add_argument(
        "--minimum-length-proportion",
        dest="lower_gene_length_threshold",
        type=float,
        default=0.5,
        help="Minimum length threshold to filter a gene from a read (default=0.5).",
    )
    parser.add_argument(
        "--maximum-length-proportion",
        dest="upper_gene_length_threshold",
        type=float,
        default=1.5,
        help="Maximum length threshold to filter a gene from a read (default=1.5).",
    )
    parser.add_argument(
        "--sample-size",
        dest="sample_size",
        type=int,
        default=100000,
        help="Number of reads to subsample to (default=100,000).",
    )
    parser.add_argument(
        "--promoter-mutations",
        dest="promoters",
        action="store_true",
        default=False,
        help="Genotype the promoter sequences of certain AMR genes (only for Escherichia_coli).",
    )
    parser.add_argument(
        "--identity",
        dest="identity",
        help="Minimum identity to a reference allele needed to report an AMR gene (default=0.9).",
        required=False,
        type=float,
        default=0.9,
    )
    parser.add_argument(
        "--coverage",
        dest="coverage",
        help="Minimum alignment coverage of a reference allele to keep an AMR gene (default=0.9).",
        required=False,
        type=float,
        default=0.9,
    )
    parser.add_argument(
        "--min-relative-depth",
        dest="min_relative_depth",
        help="Minimum relative read depth to keep an AMR gene (default=0.2).",
        required=False,
        type=float,
        default=0.2,
    )
    parser.add_argument(
        "--cores",
        dest="cores",
        type=int,
        help="Number of CPUs (default=1).",
        default=1,
    )
    parser.add_argument(
        "--pandora-path",
        dest="pandora_path",
        help="Path to pandora binary (default=pandora).",
        default="pandora",
    )
    parser.add_argument(
        "--minimap2-path",
        dest="minimap2_path",
        help="Path to minimap2 binary (default=minimap2).",
        default="minimap2",
    )
    parser.add_argument(
        "--samtools-path",
        dest="samtools_path",
        help="Path to samtools binary (default=samtools).",
        default="samtools",
    )
    parser.add_argument(
        "--racon-path",
        dest="racon_path",
        help="Path to racon binary (default=racon).",
        default="racon",
    )
    parser.add_argument(
        "--seed",
        dest="seed",
        type=int,
        help="Set the seed (default=2025).",
        default=2025,
    )
    parser.add_argument(
        "--no-sampling",
        dest="sample_reads",
        action="store_false",
        default=True,
        help="Do not randomly sample to a maximum of 100,000 input reads (default=False).",
    )
    parser.add_argument(
        "--quiet",
        dest="quiet",
        action="store_true",
        default=False,
        help="Suppress progress updates (default=False).",
    )
    parser.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        default=False,
        help="Output Amira debugging files (default=False).",
    )
    parser.add_argument(
        "--no-trim",
        dest="no_trim",
        action="store_true",
        default=False,
        help="Prevent trimming of the graph (default=False).",
    )
    parser.add_argument(
        "--output-component-fastqs",
        dest="output_components",
        action="store_true",
        default=False,
        help="Output FASTQs of the reads for each connected component (default=False).",
    )
    parser.add_argument(
        "--amr-fasta", dest="amr_fasta", help=argparse.SUPPRESS, required=False, default=None
    )
    parser.add_argument(
        "--amr-calls", dest="amr_calls", help=argparse.SUPPRESS, required=False, default=None
    )
    parser.add_argument(
        "--core-genes", dest="core_genes", help=argparse.SUPPRESS, required=False, default=None
    )
    parser.add_argument("--version", action="version", version="%(prog)s v" + __version__)
    args = parser.parse_args()
    if args.pandoraJSON and not args.gene_positions:
        parser.error("--gene-positions is required when --pandoraJSON is used.")
    return args


def write_debug_files(
    annotatedReads: dict[str, list[str]],
    geneMer_size: int,
    genesOfInterest: list[str],
    output_dir: str,
    cores: int,
) -> GeneMerGraph:
    sys.stderr.write("\nAmira: building pre-correction gene-mer graph.\n")
    raw_graph = build_multiprocessed_graph(annotatedReads, geneMer_size, cores)
    # color nodes in the graph
    for node in raw_graph.all_nodes():
        node.color_node(genesOfInterest)
    raw_graph.generate_gml(
        os.path.join(output_dir, "pre_correction_gene_mer_graph"), geneMer_size, 1, 1
    )
    return raw_graph


def main() -> None:
    # get the runtime
    start_time = time.time()
    # get command line options
    args = get_options()
    # set the seed
    random.seed(args.seed)
    # make the output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # get the relevant species-specific files
    AMR_gene_reference_FASTA, sequence_names, core_genes = load_species_specific_files(
        args.species, args.amr_fasta, args.amr_calls, args.core_genes
    )
    # import the list of genes of interest
    reference_alleles, genesOfInterest = process_reference_alleles(
        AMR_gene_reference_FASTA, args.promoters
    )
    # load the fastq data
    if not args.quiet:
        sys.stderr.write("\nAmira: loading FASTQ file.\n")
    fastq_content = parse_fastq(args.reads)
    # write out the fastq with new read IDs
    read_fastq_path, fastq_content = write_modified_fastq(
        fastq_content, args.reads, args.output_dir
    )
    # run pandora
    if args.pandoraSam is None and args.pandoraJSON is None:
        if not args.quiet:
            sys.stderr.write("\nAmira: running Pandora map.\n")
        pandoraSam, pandoraConsensus = run_pandora_map(
            args.pandora_path,
            args.panRG_path,
            read_fastq_path,
            args.output_dir,
            args.cores,
            args.seed,
        )
    else:
        pandoraSam = args.pandoraSam
        pandoraConsensus = args.pandoraConsensus
    # import a JSON of genes on reads
    if args.pandoraJSON:
        if not args.quiet:
            # output sample information
            sys.stderr.write(
                f"\nJSON file: {os.path.basename(args.pandoraJSON)}\n"
                f"FASTQ file: {os.path.basename(args.reads)}\n"
                f"Subsample reads: {args.sample_reads}\n"
                f"Trim graph: {False if args.no_trim else True}\n"
            )
        annotatedReads, sample_genesOfInterest, gene_position_dict = process_pandora_json(
            args.pandoraJSON, genesOfInterest, args.gene_positions
        )
        # sort the keys
        annotatedReads = dict(sorted(annotatedReads.items()))
        pandora_consensus = parse_fastq(args.pandoraConsensus)
        # get the mean depth across core genes
        mean_read_depth = estimate_mean_core_gene_counts(annotatedReads, core_genes)
        sys.stderr.write(f"\nAmira: mean read depth = {mean_read_depth}.\n")
    # load the pandora consensus and convert to a dictionary
    if pandoraSam:
        if not args.quiet:
            # output sample information
            sys.stderr.write(
                f"\nSAM file: {os.path.basename(pandoraSam)}\n"
                f"FASTQ file: {os.path.basename(args.reads)}\n"
                f"Subsample reads: {args.sample_reads}\n"
                f"Trim graph: {False if args.no_trim else True}\n"
            )
            sys.stderr.write("\nAmira: loading Pandora SAM file.\n")
        # load the pandora consensus
        pandora_consensus = parse_fastq(pandoraConsensus)
        # process the output of pandora
        annotatedReads, sample_genesOfInterest, gene_position_dict = convert_pandora_output(
            pandoraSam,
            pandora_consensus,
            genesOfInterest,
            args.gene_min_coverage,
            args.lower_gene_length_threshold,
            args.upper_gene_length_threshold,
            read_fastq_path,
            args.cores,
            args.output_dir,
            args.minimap2_path,
            args.samtools_path,
        )
        # sort the keys
        annotatedReads = dict(sorted(annotatedReads.items()))
        # subsample the reads
        if args.sample_reads is True:
            annotatedReads, mean_read_depth = subsample_reads_and_estimate_read_depth(
                annotatedReads, args.sample_size, args.output_dir, args.samtools_path, core_genes
            )
        # write out the gene calls
        write_pandora_gene_calls(
            args.output_dir,
            gene_position_dict,
            annotatedReads,
            os.path.join(args.output_dir, "gene_positions_with_gene_filtering.json"),
            os.path.join(args.output_dir, "gene_calls_with_gene_filtering.json"),
        )
        sys.stderr.write(f"\nAmira: mean read depth = {mean_read_depth}.\n")
    # terminate if no AMR genes were found
    if len(sample_genesOfInterest) == 0:
        # write an empty dataframe
        write_empty_result(args.output_dir)
        # exit
        sys.exit(0)
    # write out debug files if specified
    if args.debug:
        plot_read_length_distribution(annotatedReads, args.output_dir)
        write_debug_files(annotatedReads, 3, sample_genesOfInterest, args.output_dir, args.cores)
    # build the gene-mer graph
    if not args.quiet:
        sys.stderr.write("\nAmira: building intitial gene-mer graph.\n")
    graph = build_multiprocessed_graph(annotatedReads, 3, args.cores, gene_position_dict)
    # get the mean node coverages at different k-mer lengths
    overall_mean_node_coverages = get_overall_mean_node_coverages(graph)
    # collect the reads that have fewer than k genes
    short_reads = graph.get_short_read_annotations()
    short_read_gene_positions = graph.get_short_read_gene_positions()
    # mark nodes in the graph that contain at least 1 AMR gene
    if not args.no_trim:
        graph.remove_non_AMR_associated_nodes(sample_genesOfInterest)
        new_annotatedReads, new_gene_position_dict = graph.correct_reads(fastq_content)
        graph = build_multiprocessed_graph(
            new_annotatedReads, 3, args.cores, new_gene_position_dict
        )
    try:
        min_path_coverage = plot_node_coverages(
            graph.get_all_node_coverages(),
            os.path.join(args.output_dir, "initial_node_coverages.png"),
        )
    except (ValueError, IndexError):
        min_path_coverage = 10
    # filter junk reads
    graph.filter_graph(2, 1)
    new_annotatedReads, new_gene_position_dict, rejected_reads, rejected_read_positions = (
        graph.remove_junk_reads(0.80)
    )

    node_min_coverage = args.node_min_coverage
    if not args.quiet:
        message = "\nAmira: removing low coverage components "
        message += f"and nodes with coverage < {node_min_coverage}.\n"
        sys.stderr.write(message)
    # rebuild the graph
    graph = build_multiprocessed_graph(new_annotatedReads, 3, args.cores, new_gene_position_dict)
    # collect the reads that have fewer than k genes
    short_reads.update(graph.get_short_read_annotations())
    short_read_gene_positions.update(graph.get_short_read_gene_positions())
    graph.remove_low_coverage_components(5)
    graph.filter_graph(node_min_coverage, 1)
    new_annotatedReads, new_gene_position_dict = graph.correct_reads(fastq_content)

    # rebuild the graph
    graph = build_multiprocessed_graph(new_annotatedReads, 3, args.cores, new_gene_position_dict)
    # collect the reads that have fewer than k genes
    short_reads.update(graph.get_short_read_annotations())
    short_read_gene_positions.update(graph.get_short_read_gene_positions())
    graph.filter_graph(node_min_coverage, 1)
    new_annotatedReads = graph.get_valid_reads_only()
    # return an empty DF if there are no AMR genes
    if len(new_annotatedReads) == 0:
        # write an empty dataframe
        write_empty_result(args.output_dir)
        # exit
        sys.exit(0)
    # choose a value for k
    if not args.quiet:
        sys.stderr.write("\nAmira: selecting a gene-mer size (k).\n")
    geneMer_size = choose_kmer_size(
        overall_mean_node_coverages[3],
        new_annotatedReads,
        args.cores,
        new_gene_position_dict,
        sample_genesOfInterest,
    )
    overall_mean_node_coverage = overall_mean_node_coverages[geneMer_size]
    if not args.quiet:
        sys.stderr.write(f"\nAmira: selected k={geneMer_size}\n.")
        sys.stderr.write(f"\nAmira: mean node depth = {overall_mean_node_coverage}.\n")
    # correct the graph
    cleaning_iterations = 30
    new_annotatedReads, new_gene_position_dict = iterative_bubble_popping(
        new_annotatedReads,
        new_gene_position_dict,
        cleaning_iterations,
        geneMer_size,
        args.cores,
        short_reads,
        short_read_gene_positions,
        fastq_content,
        args.output_dir,
        node_min_coverage,
        sample_genesOfInterest,
        min_path_coverage,
    )
    # build the corrected gene-mer graph
    if not args.quiet:
        sys.stderr.write("\nAmira: building corrected gene-mer graph.\n")
    graph = build_multiprocessed_graph(
        new_annotatedReads, geneMer_size, args.cores, new_gene_position_dict
    )
    # write out the corrected gene calls
    write_pandora_gene_calls(
        args.output_dir,
        new_gene_position_dict,
        new_annotatedReads,
        os.path.join(args.output_dir, "corrected_gene_calls_after_filtering.json"),
        os.path.join(args.output_dir, "corrected_gene_positions_after_filtering.json"),
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
    if not args.quiet:
        sys.stderr.write("\nAmira: writing gene-mer graph.\n")
    graph.generate_gml(
        os.path.join(args.output_dir, "gene_mer_graph"),
        geneMer_size,
        node_min_coverage,
        1,
    )
    # write out a fastq of the reads in each connected component
    if args.output_components is True:
        output_component_fastqs(args.output_dir, graph, fastq_content)
    # assign reads to AMR genes by path
    if not args.quiet:
        sys.stderr.write("\nAmira: clustering reads.\n")
    (clusters_to_add, clusters_of_interest, path_reads) = process_reads(
        graph,
        sample_genesOfInterest,
        args.cores,
        short_reads,
        short_read_gene_positions,
        overall_mean_node_coverage,
    )
    # write out the fastq files
    if not os.path.exists(os.path.join(args.output_dir, "AMR_allele_fastqs")):
        os.mkdir(os.path.join(args.output_dir, "AMR_allele_fastqs"))
    # subset the fastq data based on the cluster assignments
    if not args.quiet:
        sys.stderr.write("\nAmira: writing fastqs.\n")
    (
        longest_reads_for_genes,
        supplemented_clusters_of_interest,
        allele_component_mapping,
        files_to_assemble,
    ) = write_fastqs_for_genes(
        clusters_of_interest, overall_mean_node_coverage, fastq_content, args.output_dir
    )
    # add the genes from the short reads
    longest_reads_for_genes, files_to_assemble = write_fastqs_for_genes_with_short_reads(
        clusters_to_add,
        overall_mean_node_coverage,
        longest_reads_for_genes,
        args.output_dir,
        files_to_assemble,
        fastq_content,
        supplemented_clusters_of_interest,
        allele_component_mapping,
    )
    # write out the longest reads
    longest_read_lengths = {}
    for row in longest_reads_for_genes:
        longest_read_lengths[row.split("\n")[0].replace(">", "")] = len(
            "".join(row.split("\n")[1:])
        )
    # run racon to polish the pandora consensus
    if not args.quiet:
        sys.stderr.write("\nAmira: obtaining nucleotide sequences.\n")
    result_df = get_alleles(
        files_to_assemble,
        args.racon_path,
        args.cores,
        os.path.join(args.output_dir, "AMR_allele_fastqs"),
        reference_alleles,
        pandora_consensus,
        sequence_names,
        args.debug,
        args.minimap2_path,
        args.identity,
        args.coverage,
        args.samtools_path,
    )
    # return an empty DF if there are no AMR genes
    if len(result_df) == 0:
        # write an empty dataframe
        write_empty_result(args.output_dir)
        # exit
        sys.exit(0)
    # estimate the copy numbers
    if not args.quiet:
        sys.stderr.write("\nAmira: estimating cellular copy numbers.\n")
    copy_numbers, mean_depth_per_reference = estimate_copy_numbers(
        fastq_content,
        path_reads,
        set(result_df["Amira allele"]),
        read_fastq_path,
        args.output_dir,
        args.cores,
        args.samtools_path,
        mean_read_depth,
        args.debug,
    )
    # supplement the result dataframe
    result_df = supplement_result_df(
        result_df, copy_numbers, mean_depth_per_reference, longest_read_lengths, args.debug
    )
    # get the component of each allele
    if args.output_components is True:
        result_df["Component ID"] = result_df.apply(
            lambda row: allele_component_mapping[row["Amira allele"]],
            axis=1,
        )
    # filter hits from the result df
    result_df = filter_results(
        result_df,
        args.min_relative_depth,
        supplemented_clusters_of_interest,
        annotatedReads,
        sample_genesOfInterest,
        args.identity,
        args.coverage,
        mean_read_depth,
    )
    # genotype promoters if specified
    if args.promoters:
        if not args.quiet:
            sys.stderr.write("\nAmira: genotyping promoters.\n")
        result_df = genotype_promoters(
            result_df,
            reference_alleles,
            os.path.join(args.output_dir, "AMR_allele_fastqs"),
            args.minimap2_path,
            args.racon_path,
            sequence_names,
            args.debug,
            args.output_components,
            args.samtools_path,
        )
    # write out the clustered reads
    if args.debug:
        write_reads_per_AMR_gene(args.output_dir, supplemented_clusters_of_interest)
    # sort the results
    result_df = result_df.sort_values(by="Determinant name")
    # write the result tsv
    result_df.to_csv(os.path.join(args.output_dir, "amira_results.tsv"), sep="\t", index=False)
    if not args.quiet:
        # display the runtime
        sys.stderr.write(f"\nAmira: Total runtime {round(time.time() - start_time)} seconds.\n")
    sys.exit(0)


if __name__ == "__main__":
    main()
