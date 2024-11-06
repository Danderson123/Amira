import argparse
import json
import os
import random
import sys
import time

import matplotlib
from tqdm import tqdm

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
    process_pandora_json,
    process_reference_alleles,
)
from amira.read_operations import (
    downsample_reads,
    parse_fastq,
    plot_read_length_distribution,
    write_fastq,
)
from amira.result_operations import estimate_copy_number, process_reads, write_allele_fastq

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
    reference_alleles, genesOfInterest = process_reference_alleles(args.path_to_interesting_genes)
    # import a JSON of genes on reads
    if args.pandoraJSON:
        if not args.quiet:
            # output sample information
            sys.stderr.write(
                f"Sample name: {os.path.basename(os.path.dirname(args.pandoraJSON))}\n"
                f"JSON file: {os.path.basename(args.pandoraJSON)}\n"
                f"FASTQ path: {os.path.basename(args.readfile)}\n"
                f"Subsample reads: {args.sample_reads}\n"
                f"Trim graph: {False if args.no_trim else True}\n"
            )
        annotatedReads, sample_genesOfInterest, gene_position_dict = process_pandora_json(
            args.pandoraJSON, genesOfInterest, args.gene_positions
        )
        pandora_consensus = parse_fastq(args.pandoraConsensus)
        # randomly sample 100,000 reads
        if args.sample_reads:
            annotatedReads = downsample_reads(annotatedReads, 10000)
    # load the pandora consensus and convert to a dictionary
    if args.pandoraSam:
        if not args.quiet:
            # output sample information
            sys.stderr.write(
                f"Sample name: {os.path.basename(os.path.dirname(args.pandoraSam))}\n"
                f"SAM file: {os.path.basename(args.pandoraSam)}\n"
                f"FASTQ path: {os.path.basename(args.readfile)}\n"
                f"Subsample reads: {args.sample_reads}\n"
                f"Trim graph: {False if args.no_trim else True}\n"
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
    if not args.quiet:
        sys.stderr.write("\nAmira: building intitial gene-mer graph\n")
    graph = build_multiprocessed_graph(
        annotatedReads, args.geneMer_size, args.cores, gene_position_dict
    )
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
            new_annotatedReads, args.geneMer_size, args.cores, new_gene_position_dict
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
    graph = build_multiprocessed_graph(
        new_annotatedReads, args.geneMer_size, args.cores, new_gene_position_dict
    )
    # collect the reads that have fewer than k genes
    short_reads.update(graph.get_short_read_annotations())
    short_read_gene_positions.update(graph.get_short_read_gene_positions())
    graph.remove_low_coverage_components(5)
    graph.filter_graph(node_min_coverage, 1)
    new_annotatedReads, new_gene_position_dict = graph.correct_reads(fastq_content)
    # rebuild the graph
    graph = build_multiprocessed_graph(
        new_annotatedReads, args.geneMer_size, args.cores, new_gene_position_dict
    )
    # collect the reads that have fewer than k genes
    short_reads.update(graph.get_short_read_annotations())
    short_read_gene_positions.update(graph.get_short_read_gene_positions())
    graph.filter_graph(node_min_coverage, 1)
    new_annotatedReads = graph.get_valid_reads_only()
    # return an empty DF if there are no AMR genes
    if len(new_annotatedReads) == 0:
        # write an empty dataframe
        results = "Gene name\tSequence name\tClosest reference\tReference length\t"
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
        overall_mean_node_coverages[args.geneMer_size],
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
    if not os.path.exists(os.path.join(args.output_dir, "component_fastqs")):
        os.mkdir(os.path.join(args.output_dir, "component_fastqs"))
    for component in graph.components():
        node_hashes_in_component = [n.__hash__() for n in graph.get_nodes_in_component(component)]
        reads_in_component = graph.collect_reads_in_path(node_hashes_in_component)
        component_fastq_data = {r: fastq_content[r] for r in reads_in_component}
        write_fastq(
            os.path.join(args.output_dir, "component_fastqs", f"{component}.fastq.gz"),
            component_fastq_data,
        )
    # assign reads to AMR genes by path
    if not args.quiet:
        sys.stderr.write("\nAmira: clustering reads\n")
    (
        clusters_to_add,
        cluster_copy_numbers_to_add,
        clusters_of_interest,
        cluster_copy_numbers,
        allele_counts,
    ) = process_reads(
        graph,
        sample_genesOfInterest,
        fastq_content,
        short_reads,
        short_read_gene_positions,
        overall_mean_node_coverage,
    )
    # write out the fastq files
    if not os.path.exists(os.path.join(args.output_dir, "AMR_allele_fastqs")):
        os.mkdir(os.path.join(args.output_dir, "AMR_allele_fastqs"))
    # subset the fastq data based on the cluster assignments
    files_to_assemble = []
    if not args.quiet:
        sys.stderr.write("\nAmira: writing fastqs\n")
    supplemented_clusters_of_interest = {}
    allele_component_mapping = {}
    for component in tqdm(clusters_of_interest):
        for gene in clusters_of_interest[component]:
            for allele in clusters_of_interest[component][gene]:
                if (
                    len(clusters_of_interest[component][gene][allele])
                    > overall_mean_node_coverage / 20
                ):
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
                    # store the component of the allele
                    allele_component_mapping[allele] = component
                else:
                    message = f"\nAmira: allele {allele} removed "
                    message += "due to an insufficient number of reads "
                    message += f"({len(clusters_of_interest[component][gene][allele])}).\n"
                    sys.stderr.write(message)
    # add the genes from the short reads
    for allele in clusters_to_add:
        if len(clusters_to_add[allele]) > overall_mean_node_coverage / 20:
            files_to_assemble.append(
                write_allele_fastq(clusters_to_add[allele], fastq_content, args.output_dir, allele)
            )
            supplemented_clusters_of_interest[allele] = clusters_to_add[allele]
            allele_component_mapping[allele] = None
        else:
            message = f"\nAmira: allele {allele} removed "
            message += f"due to an insufficient number of reads ({len(clusters_to_add[allele])}).\n"
            sys.stderr.write(message)
    # run racon to polish the pandora consensus
    if not args.quiet:
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
        args.minimap2_path,
    )
    # return an empty DF if there are no AMR genes
    if len(result_df) == 0:
        # write an empty dataframe
        results = "Gene name\tSequence name\tClosest reference\tReference length\t"
        results += (
            "Identity (%)\tCoverage (%)\tAmira allele\tNumber of reads\tApproximate copy number\n"
        )
        with open(os.path.join(args.output_dir, "amira_results.tsv"), "w") as o:
            o.write(results)
        # exit
        sys.exit(0)
    # get the copy number estimates of each allele
    result_df["Approximate copy number"] = result_df.apply(
        lambda row: estimate_copy_number(
            row["Amira allele"], cluster_copy_numbers, cluster_copy_numbers_to_add
        ),
        axis=1,
    )
    # get the component of each allele
    result_df["Component ID"] = result_df.apply(
        lambda row: allele_component_mapping[row["Amira allele"]],
        axis=1,
    )
    # remove genes that do not have sufficient mapping coverage
    alleles_to_delete = []
    for index, row in result_df.iterrows():
        if isinstance(row["Identity (%)"], str) and "/" in row["Identity (%)"]:
            identity = float(row["Identity (%)"].split("/")[0])
        else:
            identity = row["Identity (%)"]
        if identity < 90:
            message = f"\nAmira: allele {row['Amira allele']} removed "
            message += f"due to insufficient similarity ({identity}).\n"
            sys.stderr.write(message)
            alleles_to_delete.append(row["Amira allele"])
            continue
        else:
            if isinstance(row["Coverage (%)"], str) and "/" in row["Coverage (%)"]:
                coverage = float(row["Coverage (%)"].split("/")[0])
            else:
                coverage = row["Coverage (%)"]
            if coverage < 90:
                message = f"\nAmira: allele {row['Amira allele']} removed "
                message += f"due to insufficient coverage ({coverage}).\n"
                sys.stderr.write(message)
                alleles_to_delete.append(row["Amira allele"])
                continue
        # check if filter contaminants is on
        if args.filter_contamination is True:
            # remove alleles where all of the reads just contain AMR genes
            reads = supplemented_clusters_of_interest[row["Amira allele"]]
            if all(
                all(g[1:] in sample_genesOfInterest for g in annotatedReads[r.split("_")[0]])
                for r in reads
            ):
                alleles_to_delete.append(row["Amira allele"])
                message = f"\nAmira: allele {row['Amira allele']} removed "
                message += "due to suspected contamination.\n"
                sys.stderr.write(message)
                continue
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
    if not args.quiet:
        # display the runtime
        sys.stderr.write(f"\nAmira: Total runtime {round(time.time() - start_time)} seconds\n")
    sys.exit(0)


if __name__ == "__main__":
    main()
