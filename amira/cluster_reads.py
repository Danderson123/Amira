import statistics
import math
from amira.construct_graph import GeneMerGraph
import numpy as np
from scipy.special import comb


def get_relevant_reads(calls, genesofinterest):
    return {read for read in calls if any(g[1:] in genesofinterest for g in calls[read])}

def remove_irrelevant_nodes(graph, relevant_reads):
    nodes_to_delete = [graph.get_node_by_hash(node) for node in graph.get_nodes() if not any(r in relevant_reads for r in graph.get_node_by_hash(node).get_reads())]
    for node in nodes_to_delete:
        graph.remove_node(node)

def get_call_subset(graph, calls):
    return {read: calls[read] for read in graph.get_readNodes() if not all(n == None for n in graph.get_readNodes()[read])}

def calculate_coverage_probabilities(component_reads, gene_coverages):
    # calculate the probability that at least 10 reads cover each node
    mean_coverage = (len(component_reads) * statistics.mean(component_reads.values())) / len(gene_coverages)
    prob_X_ge_10 = 1 - (math.exp(-mean_coverage) * sum(mean_coverage ** k / math.factorial(k) for k in range(10)))
    return prob_X_ge_10 ** len(gene_coverages)

def get_mean_node_coverage(input_graph):
    return input_graph.get_mean_node_coverage()

def process_component(component, input_graph, new_annotatedReads, mean_node_coverage, copy_numbers, reads_in_component, junction_tracking):
    nodes_in_component = input_graph.get_nodes_in_component(component)
    node_coverages, junction = analyze_component_nodes(nodes_in_component, input_graph)
    copy_numbers[component] = estimate_copy_number(node_coverages, mean_node_coverage)
    reads_in_component[component] = collect_component_reads(nodes_in_component, new_annotatedReads)
    junction_tracking[component] = junction

def analyze_component_nodes(nodes_in_component, input_graph):
    junction = False
    node_coverages = []
    for node in nodes_in_component:
        if not junction and (len(input_graph.get_forward_neighbors(node)) > 1 or len(input_graph.get_backward_neighbors(node)) > 1):
            junction = True
        node_coverages.append(node.get_node_coverage())
    return node_coverages, junction

def estimate_copy_number(node_coverages, mean_node_coverage):
    if node_coverages:
        return max(1.0, statistics.mean(node_coverages) / mean_node_coverage)
    return 1.0

def collect_component_reads(nodes_in_component, new_annotatedReads):
    component_reads = {}
    for node in nodes_in_component:
        for read in node.get_reads():
            component_reads[read] = new_annotatedReads.get(read, {})
    return component_reads

# Generate the probability distribution for the number of successful trials (i.e., gene appearing in a read of length at least k + k - 1)
def binom_pmf(k, n, p):
    """
    Calculate the probability mass function for the binomial distribution.
    
    Args:
    k (int): Number of successes.
    n (int): Number of trials.
    p (float): Probability of success in each trial.
    
    Returns:
    float: Probability of exactly k successes.
    """
    return comb(n, k) * (p ** k) * ((1 - p) ** (n - k))

def prob_at_least_10_successes(n, p):
    """
    Calculate the probability of getting at least 10 successes in a binomial distribution.
    
    Args:
    n (int): Number of trials.
    p (float): Probability of success in each trial.
    
    Returns:
    float: Probability of getting at least 10 successes.
    """
    # Calculate the cumulative probability of getting fewer than 10 successes
    cumulative_prob = np.sum([binom_pmf(k, n, p) for k in range(15)])
    # Probability of getting at least 10 successes is 1 minus this cumulative probability
    prob_at_least_10 = 1 - cumulative_prob
    return prob_at_least_10

def gene_read_probability_np(R, genes_per_read, reads_per_gene, k):
    """
    Calculate the probability that each unique gene occurs in at least 10 reads of length at least k using NumPy.
    Args:
    R (int): Total number of reads.
    genes_per_read (dict): Dictionary of read to list of genes.
    reads_per_gene (dict): Dictionary of gene to list of reads.
    k (int): Threshold for the minimum read length.
    Returns:
    float: Probability of a gene occurring in at least 10 reads of sufficient length.
    """
    # Number of unique genes
    unique_genes = set(reads_per_gene.keys())
    # Probability that a read is of length at least k
    p_read_length_ge_k = sum(1 for r in genes_per_read if len(genes_per_read[r]) >= (k + k - 1)) / R
    probabilities = {}
    for gene in unique_genes:
        num_reads_with_gene = reads_per_gene[gene]
        # Calculate the probability that the gene appears in at least 10 reads of sufficient length
        prob_at_least_10 = prob_at_least_10_successes(num_reads_with_gene, p_read_length_ge_k)
        probabilities[gene] = prob_at_least_10
    return probabilities

def build_subset_graph_and_assign_reads(aggregated_clusters_of_interest, component, component_reads, k, gene_positions, mean_node_coverage, sample_genesOfInterest, fastq_content, allele_counts):
    new_graph = GeneMerGraph(component_reads, k, gene_positions)
    clusters_of_interest, cluster_copy_numbers, allele_counts = new_graph.new_assign_reads_to_genes(
                                                sample_genesOfInterest, fastq_content, allele_counts, mean_node_coverage
                                            )
    for c in clusters_of_interest:
        aggregated_clusters_of_interest[component].update(clusters_of_interest[c])

def cluster_reads(sample_genesOfInterest, fastq_content, input_graph, new_annotatedReads, new_gene_position_dict):
    mean_node_coverage = get_mean_node_coverage(input_graph)
    copy_numbers = {}
    reads_in_component = {}
    junction_tracking = {}
    for component in input_graph.components():
        process_component(component, input_graph, new_annotatedReads, mean_node_coverage, copy_numbers, reads_in_component, junction_tracking)
    # iterate through the components
    current_k = 3
    aggregated_clusters_of_interest = {}
    allele_counts = {}
    # iterate through the components in the input graph
    for component in junction_tracking:
        aggregated_clusters_of_interest[component] = {}
        component_reads = reads_in_component[component].copy()
        if junction_tracking[component] is False:
            build_subset_graph_and_assign_reads(aggregated_clusters_of_interest,
                                            component,
                                            component_reads,
                                            current_k,
                                            new_gene_position_dict,
                                            mean_node_coverage,
                                            sample_genesOfInterest,
                                            fastq_content,
                                            allele_counts)
        while junction_tracking[component] is True and current_k < 11:
            current_k += 2
            new_graph = GeneMerGraph(component_reads, current_k, new_gene_position_dict)
            new_component_junction_tracking = {}
            for other_component in new_graph.components():
                junction = False
                other_component_reads = {}
                for node in new_graph.get_nodes_in_component(other_component):
                    if junction is False:
                        if len(new_graph.get_forward_neighbors(node)) > 1 or len(new_graph.get_backward_neighbors(node)) > 1:
                            junction = True
                    for read in node.get_reads():
                        if None not in new_graph.get_readNodes()[read]:
                            other_component_reads[read] = component_reads[read]
                gene_coverages = {}
                for r in other_component_reads:
                    for g in other_component_reads[r]:
                        if g[1:] in sample_genesOfInterest:
                            if g[1:] not in gene_coverages:
                                gene_coverages[g[1:]] = 0
                            gene_coverages[g[1:]] += 1
                if junction is True:
                    print(component, other_component, junction, current_k, len(other_component_reads), len(gene_coverages), gene_read_probability_np(len(other_component_reads), other_component_reads, gene_coverages, current_k))
                new_component_junction_tracking[other_component] = junction
                if junction is False or current_k == 11:
                    build_subset_graph_and_assign_reads(aggregated_clusters_of_interest,
                                                        component,
                                                        other_component_reads,
                                                        current_k,
                                                        new_gene_position_dict,
                                                        mean_node_coverage,
                                                        sample_genesOfInterest,
                                                        fastq_content,
                                                        allele_counts)
                    for r in other_component_reads:
                        del component_reads[r]
            #print(component, current_k, new_component_junction_tracking, len(component_reads))
            #for g in aggregated_clusters_of_interest[component]:
            #    print(g, len(aggregated_clusters_of_interest[component][g]))
            if all(new_component_junction_tracking[other_component] is False for other_component in new_component_junction_tracking):
                junction_tracking[component] = False
    dkdkkd

def cluster_reads_old(sample_genesOfInterest, fastq_content, input_graph, new_annotatedReads, new_gene_position_dict):
    reads_in_component = {}
    mean_node_coverage = input_graph.get_mean_node_coverage()
    copy_numbers = {}
    # iterate through the components
    for component in input_graph.components():
        reads_in_component[component] = {}
        # get the nodes in this component
        nodes_in_component = input_graph.get_nodes_in_component(component)
        k_vals = [3, 5, 7, 9, 11, 13, 15, 17, 31]
        gene_coverages = {k: {} for k in k_vals}
        read_lengths = {k: {} for k in k_vals}
        node_coverages = []
        # iterate through the nodes in this component
        for node in nodes_in_component:
            # collect the coverages of the nodes in this component
            node_coverages.append(node.get_node_coverage())
            # iterate through the reads for this node
            for read in node.get_reads():
                # subset the annotations of reads in this component
                reads_in_component[component][read] = new_annotatedReads[read]
                # iterate through the k values
                for k in k_vals:
                    if len(new_annotatedReads[read]) >= k:
                        # get the length of the reads in this component in genes
                        read_lengths[k][read] = len(new_annotatedReads[read])
                    # iterate through the genes on the read
                    for gene in new_annotatedReads[read]:
                        if gene[1:] not in gene_coverages[k]:
                            gene_coverages[k][gene[1:]] = 0
                        if len(new_annotatedReads[read]) >= k:
                            gene_coverages[k][gene[1:]] += 1
        # estimate the copy number of this component
        copy_numbers[component] = max(1.0, statistics.mean(node_coverages) / mean_node_coverage)
        # calculate the probability that each gene is found in 5 genes
        probabilities = {}
        for k in k_vals:
            prob_all_bases_covered = calculate_coverage_probabilities(read_lengths[k], list(gene_coverages[k].values()))
            probabilities[k] = prob_all_bases_covered
        print(component, probabilities, "\n")
    jduduud
    # initialise dictionaries to track thmin(path_coverages[path]) / mean_node_coveragee alleles
    aggregated_clusters_of_interest = {7: {}, 3: {}}
    aggregated_cluster_copy_numbers = {7: {}, 3: {}}
    allele_counts = {}
    for c in reads_in_component:
        subsetted_graph = GeneMerGraph(reads_in_component[c], k, new_gene_position_dict)
        clusters_of_interest, cluster_copy_numbers, allele_counts = subsetted_graph.new_assign_reads_to_genes(
                                                    sample_genesOfInterest, fastq_content, allele_counts, mean_node_coverage
                                                )
        for c in clusters_of_interest:
            aggregated_clusters_of_interest[k].update(clusters_of_interest[c])
            aggregated_cluster_copy_numbers[k].update(cluster_copy_numbers[c])
    for k in aggregated_clusters_of_interest:
        for g in aggregated_clusters_of_interest[k]:
            for a in aggregated_clusters_of_interest[k][g]:
                print(k, g, a, len(aggregated_clusters_of_interest[k][g][a]))
    print("\n")
    return aggregated_clusters_of_interest, aggregated_cluster_copy_numbers
