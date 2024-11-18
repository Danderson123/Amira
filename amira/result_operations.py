import os

from amira.read_operations import write_fastq


def get_found_genes(clusters_of_interest):
    found_genes = set()
    for component_id in clusters_of_interest:
        for gene in clusters_of_interest[component_id]:
            found_genes.add(gene)
    return found_genes


def add_amr_alleles(short_reads, short_read_gene_positions, sample_genesOfInterest, found_genes):
    clusters_to_add = {}
    for read_id in short_reads:
        for g in range(len(short_reads[read_id])):
            strandless_gene = short_reads[read_id][g][1:]
            if strandless_gene in sample_genesOfInterest and strandless_gene not in found_genes:
                if f"{strandless_gene}_1" not in clusters_to_add:
                    clusters_to_add[f"{strandless_gene}_1"] = []
                gene_start, gene_end = short_read_gene_positions[read_id][g]
                clusters_to_add[f"{strandless_gene}_1"].append(f"{read_id}_{gene_start}_{gene_end}")
    return clusters_to_add


def calculate_cluster_copy_numbers(clusters_to_add, overall_mean_node_coverage):
    cluster_copy_numbers = {}
    for allele in clusters_to_add:
        cluster_copy_numbers[allele] = max(
            1.0, len(clusters_to_add[allele]) / overall_mean_node_coverage
        )
    return cluster_copy_numbers


def process_reads(
    graph,
    sample_genesOfInterest,
    fastq_content,
    short_reads,
    short_read_gene_positions,
    overall_mean_node_coverage,
):
    # Step 1: Assign reads to genes
    clusters_of_interest, cluster_copy_numbers, allele_counts = graph.assign_reads_to_genes(
        sample_genesOfInterest, fastq_content, {}, overall_mean_node_coverage
    )
    # Step 2: Get the unique genes found in clusters of interest
    found_genes_of_interest = get_found_genes(clusters_of_interest)
    # Step 3: Add AMR alleles that come from short reads
    clusters_to_add = add_amr_alleles(
        short_reads, short_read_gene_positions, sample_genesOfInterest, found_genes_of_interest
    )
    # Step 4: Calculate cluster copy numbers for newly added clusters
    cluster_copy_numbers_to_add = calculate_cluster_copy_numbers(
        clusters_to_add, overall_mean_node_coverage
    )
    # Return the processed results
    return (
        clusters_to_add,
        cluster_copy_numbers_to_add,
        clusters_of_interest,
        cluster_copy_numbers,
        allele_counts,
    )


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
