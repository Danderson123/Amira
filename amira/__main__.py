import argparse
import json
import os
import random
import sys
import time

import matplotlib

from amira.__init__ import __version__
from amira.construct_graph import GeneMerGraph
from amira.graph_operations import (
    build_multiprocessed_graph,
    choose_kmer_size,
    get_overall_mean_node_coverages,
    iterative_bubble_popping,
    plot_node_coverages,
)
from amira.pre_process_operations import (
    convert_pandora_output,
    get_core_gene_mean_depth,
    process_pandora_json,
    process_reference_alleles,
)
from amira.read_operations import downsample_reads, parse_fastq, plot_read_length_distribution
from amira.result_operations import (
    estimate_copy_numbers,
    filter_results,
    genotype_promoters,
    get_alleles,
    output_component_fastqs,
    process_reads,
    write_fastqs_for_genes,
    write_fastqs_for_genes_with_short_reads,
    write_reads_per_AMR_gene,
)

matplotlib.use("Agg")


def get_options() -> argparse.Namespace:
    """define args from the command line"""
    parser = argparse.ArgumentParser(
        description="Identify acquired AMR genes from bacterial long read sequences."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pandoraSam", dest="pandoraSam", help="Pandora map SAM file path.")
    group.add_argument("--pandoraJSON", dest="pandoraJSON", help="Pandora map JSON file path.")
    parser.add_argument("--gene-positions", help="Gene position JSON file path.")
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
        "--promoter-mutations",
        dest="promoters",
        action="store_true",
        default=False,
        help="Report point mutations in promoters that are associated with AMR.",
    )
    parser.add_argument(
        "--annotation",
        dest="phenotypes",
        help="Path to a JSON of phenotypes for each AMR allele.",
        required=True,
    )
    parser.add_argument(
        "--core-genes",
        dest="core_genes",
        help="Path to a newline delimited list of core genes in this species.",
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
        help="Path to racon binary.",
        default="racon",
    )
    parser.add_argument(
        "--minimap2-path",
        dest="minimap2_path",
        help="Path to minimap2 binary.",
        default="minimap2",
    )
    parser.add_argument(
        "--seed",
        dest="seed",
        type=int,
        help="Set the seed.",
        default=2024,
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
        "--filter-contaminants",
        dest="filter_contamination",
        action="store_true",
        default=False,
        help="Filter AMR alleles that are suspected contaminants.",
    )
    parser.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        default=False,
        help="Output Amira debugging files.",
    )
    parser.add_argument(
        "--no-trim",
        dest="no_trim",
        action="store_true",
        default=False,
        help="Prevent trimming of the graph.",
    )
    parser.add_argument(
        "--component-fastqs",
        dest="output_components",
        action="store_true",
        default=False,
        help="Output FASTQs of the reads for each connected component in the graph.",
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
    sys.stderr.write("\nAmira: building pre-correction gene-mer graph\n")
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
    # import the list of genes of interest
    reference_alleles, genesOfInterest = process_reference_alleles(
        args.path_to_interesting_genes, args.promoters
    )
    # import a JSON of genes on reads
    if args.pandoraJSON:
        if not args.quiet:
            # output sample information
            sys.stderr.write(
                f"Sample name: {os.path.basename(os.path.dirname(args.pandoraJSON))}\n"
                f"JSON file: {os.path.basename(args.pandoraJSON)}\n"
                f"FASTQ file: {os.path.basename(args.readfile)}\n"
                f"Subsample reads: {args.sample_reads}\n"
                f"Trim graph: {False if args.no_trim else True}\n"
                f"Detect promoter SNPs: {args.promoters}\n"
            )
        annotatedReads, sample_genesOfInterest, gene_position_dict = process_pandora_json(
            args.pandoraJSON, genesOfInterest, args.gene_positions
        )
        pandora_consensus = parse_fastq(args.pandoraConsensus)
        # randomly sample 100,000 reads
        if args.sample_reads:
            annotatedReads = downsample_reads(annotatedReads, 10000)
        # initialise mean read depth
        mean_read_depth = None
    # load the pandora consensus and convert to a dictionary
    if args.pandoraSam:
        if not args.quiet:
            # output sample information
            sys.stderr.write(
                f"Sample name: {os.path.basename(os.path.dirname(args.pandoraSam))}\n"
                f"SAM file: {os.path.basename(args.pandoraSam)}\n"
                f"FASTQ file: {os.path.basename(args.readfile)}\n"
                f"Subsample reads: {args.sample_reads}\n"
                f"Trim graph: {False if args.no_trim else True}\n"
                f"Detect promoter SNPs: {args.promoters}\n"
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
            args.minimap2_path,
        )
        with open(
            os.path.join(args.output_dir, "gene_positions_with_gene_filtering.json"), "w"
        ) as o:
            o.write(json.dumps(gene_position_dict))
        with open(os.path.join(args.output_dir, "gene_calls_with_gene_filtering.json"), "w") as o:
            o.write(json.dumps(annotatedReads))
        # randomly sample 100,000 reads
        if args.sample_reads:
            annotatedReads = downsample_reads(annotatedReads, 100000)
        # get the mean depth across core genes
        mean_read_depth = get_core_gene_mean_depth(
            os.path.join(args.output_dir, "mapped_to_consensus.bam"), args.core_genes
        )
        sys.stderr.write(f"\nAmira: mean read depth = {mean_read_depth}\n")
    # terminate if no AMR genes were found
    if len(sample_genesOfInterest) == 0:
        # write an empty dataframe
        results = "Determinant name\tSequence name\tClosest reference\tReference length\t"
        results += (
            "Identity (%)\tCoverage (%)\tAmira allele\tNumber of reads\tApproximate copy number\n"
        )
        with open(os.path.join(args.output_dir, "amira_results.tsv"), "w") as o:
            o.write(results)
        # exit
        sys.exit(0)
    # load the fastq data
    fastq_content = parse_fastq(args.readfile)
    fastq_content = {r.replace("_", ""): fastq_content[r] for r in fastq_content}
    # write out debug files if specified
    if args.debug:
        plot_read_length_distribution(annotatedReads, args.output_dir)
        write_debug_files(annotatedReads, 3, sample_genesOfInterest, args.output_dir, args.cores)
    # build the gene-mer graph
    if not args.quiet:
        sys.stderr.write("\nAmira: building intitial gene-mer graph\n")
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
        message += f"and nodes with coverage < {node_min_coverage}\n"
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
        results = "Determinant name\tSequence name\tClosest reference\tReference length\t"
        results += (
            "Identity (%)\tCoverage (%)\tAmira allele\tNumber of reads\tApproximate copy number\n"
        )
        with open(os.path.join(args.output_dir, "amira_results.tsv"), "w") as o:
            o.write(results)
        # exit
        sys.exit(0)
    # choose a value for k
    if not args.quiet:
        sys.stderr.write("\nAmira: selecting a gene-mer size (k)\n")
    geneMer_size = choose_kmer_size(
        overall_mean_node_coverages[3],
        new_annotatedReads,
        args.cores,
        new_gene_position_dict,
        sample_genesOfInterest,
    )
    overall_mean_node_coverage = overall_mean_node_coverages[geneMer_size]
    if not args.quiet:
        sys.stderr.write(f"\nAmira: selected k={geneMer_size}\n")
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
    if not args.quiet:
        sys.stderr.write("\nAmira: writing gene-mer graph\n")
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
        sys.stderr.write("\nAmira: clustering reads\n")
    (
        clusters_to_add,
        clusters_of_interest,
    ) = process_reads(
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
        sys.stderr.write("\nAmira: writing fastqs\n")
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
    with open(os.path.join(args.output_dir, "AMR_allele_fastqs", "longest_reads.fasta"), "w") as o:
        o.write("\n".join(longest_reads_for_genes))
    # run racon to polish the pandora consensus
    if not args.quiet:
        sys.stderr.write("\nAmira: obtaining nucleotide sequences\n")
    result_df = get_alleles(
        files_to_assemble,
        args.racon_path,
        args.cores,
        os.path.join(args.output_dir, "AMR_allele_fastqs"),
        reference_alleles,
        pandora_consensus,
        args.phenotypes,
        args.debug,
        args.minimap2_path,
    )
    # return an empty DF if there are no AMR genes
    if len(result_df) == 0:
        # write an empty dataframe
        results = "Determinant name\tSequence name\tClosest reference\tReference length\t"
        results += (
            "Identity (%)\tCoverage (%)\tAmira allele\tNumber of reads\tApproximate copy number\n"
        )
        with open(os.path.join(args.output_dir, "amira_results.tsv"), "w") as o:
            o.write(results)
        # exit
        sys.exit(0)
    # estimate the copy numbers
    if not args.quiet:
        sys.stderr.write("\nAmira: estimating cellular copy numbers\n")
    copy_numbers = estimate_copy_numbers(
        mean_read_depth,
        os.path.join(args.output_dir, "AMR_allele_fastqs", "longest_reads.fasta"),
        args.readfile,
        args.cores,
    )
    estimates = []
    for index, row in result_df.iterrows():
        estimates.append(copy_numbers[row["Amira allele"]])
    result_df["Approximate copy number"] = estimates
    # get the component of each allele
    if args.output_components is True:
        result_df["Component ID"] = result_df.apply(
            lambda row: allele_component_mapping[row["Amira allele"]],
            axis=1,
        )
    # filter hits from the result df
    result_df = filter_results(
        result_df,
        args.filter_contamination,
        supplemented_clusters_of_interest,
        annotatedReads,
        sample_genesOfInterest,
    )
    # genotype promoters if specified
    if args.promoters:
        if not args.quiet:
            sys.stderr.write("\nAmira: genotyping promoters\n")
        result_df = genotype_promoters(
            result_df,
            reference_alleles,
            os.path.join(args.output_dir, "AMR_allele_fastqs"),
            args.minimap2_path,
            args.racon_path,
            args.phenotypes,
            args.debug,
            args.output_components,
        )
    # write out the clustered reads
    write_reads_per_AMR_gene(args.output_dir, supplemented_clusters_of_interest)
    # sort the results
    result_df = result_df.sort_values(by="Determinant name")
    # write the result tsv
    result_df.to_csv(os.path.join(args.output_dir, "amira_results.tsv"), sep="\t", index=False)
    if not args.quiet:
        # display the runtime
        sys.stderr.write(f"\nAmira: Total runtime {round(time.time() - start_time)} seconds\n")
    sys.exit(0)


if __name__ == "__main__":
    main()
