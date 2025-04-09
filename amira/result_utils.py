import json
import os
import shutil
import statistics
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from joblib import Parallel, delayed
from scipy.optimize import minimize
from scipy.signal import find_peaks, savgol_filter
from scipy.stats import poisson
from tqdm import tqdm

from amira.read_utils import write_fastq


def get_found_genes(clusters_of_interest):
    found_genes = set()
    for component_id in clusters_of_interest:
        for gene in clusters_of_interest[component_id]:
            found_genes.add(gene)
    return found_genes


def add_amr_alleles(
    short_reads, short_read_gene_positions, sample_genesOfInterest, found_genes, path_reads
):
    clusters_to_add = {}
    for read_id in short_reads:
        for g in range(len(short_reads[read_id])):
            strandless_gene = short_reads[read_id][g][1:]
            if strandless_gene in sample_genesOfInterest and strandless_gene not in found_genes:
                if f"{strandless_gene}_1" not in clusters_to_add:
                    clusters_to_add[f"{strandless_gene}_1"] = []
                gene_start, gene_end = short_read_gene_positions[read_id][g]
                clusters_to_add[f"{strandless_gene}_1"].append(f"{read_id}_{gene_start}_{gene_end}")
                path_tuple = tuple([f"+{strandless_gene}_1"])
                if path_tuple not in path_reads:
                    path_reads[path_tuple] = set()
                path_reads[path_tuple].add(read_id)
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
    cores,
    short_reads,
    short_read_gene_positions,
    overall_mean_node_coverage,
):
    # Step 1: Assign reads to genes
    clusters_of_interest, path_reads = graph.assign_reads_to_genes(
        sample_genesOfInterest, cores, {}, overall_mean_node_coverage
    )
    # Step 2: Get the unique genes found in clusters of interest
    found_genes_of_interest = get_found_genes(clusters_of_interest)
    # Step 3: Add AMR alleles that come from short reads
    clusters_to_add = add_amr_alleles(
        short_reads,
        short_read_gene_positions,
        sample_genesOfInterest,
        found_genes_of_interest,
        path_reads,
    )
    # Return the processed results
    return (clusters_to_add, clusters_of_interest, path_reads)


def write_path_fastq(reads_for_path, fastq_content, output_dir, path_id):
    read_subset = {}
    for r in reads_for_path:
        fastq_data = fastq_content[r].copy()
        if fastq_data["sequence"] != "":
            read_subset[r] = fastq_data
    if not os.path.exists(os.path.join(output_dir)):
        os.mkdir(os.path.join(output_dir))
    write_fastq(
        os.path.join(output_dir, f"{path_id}.fastq.gz"),
        read_subset,
    )
    return os.path.join(output_dir, f"{path_id}.fastq.gz")


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


def filter_results(
    result_df,
    min_relative_depth,
    supplemented_clusters_of_interest,
    annotatedReads,
    sample_genesOfInterest,
    required_identity,
    required_coverage,
    mean_read_depth,
):
    # remove genes that do not have sufficient mapping coverage
    alleles_to_delete = []
    comments = []
    # skip filtering by copy number if mean read depth is <20x
    if mean_read_depth < 20:
        skip_depth_filtering = True
        message = "\nAmira: skipping filtering by depth as read depth <20x.\n"
        sys.stderr.write(message)
    else:
        skip_depth_filtering = False
    # modify required coverage to make a percentage
    required_coverage = required_coverage * 100
    for index, row in result_df.iterrows():
        # store comments about this allele
        flags = []
        if isinstance(row["Identity (%)"], str) and "/" in row["Identity (%)"]:
            identity = float(row["Identity (%)"].split("/")[0])
        else:
            identity = row["Identity (%)"]
        if identity < required_identity:
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
            if coverage < required_coverage:
                message = f"\nAmira: allele {row['Amira allele']} removed "
                message += f"due to insufficient coverage ({coverage}).\n"
                sys.stderr.write(message)
                alleles_to_delete.append(row["Amira allele"])
                continue
            else:
                if skip_depth_filtering is False:
                    relative_depth = row["Mean read depth"] / mean_read_depth
                    if relative_depth < min_relative_depth:
                        message = f"\nAmira: allele {row['Amira allele']} removed "
                        message += "due to insufficient relative read depth "
                        message += f"({relative_depth}).\n"
                        sys.stderr.write(message)
                        alleles_to_delete.append(row["Amira allele"])
                        continue
                if coverage < 90:
                    flags.append("Partially present gene.")
        # remove alleles where all of the reads just contain AMR genes
        reads = supplemented_clusters_of_interest[row["Amira allele"]]
        if all(
            all(g[1:] in sample_genesOfInterest for g in annotatedReads[r.split("_")[0]])
            for r in reads
        ):
            flags.append("Potential contaminant.")
        # collect the flags
        comments.append(" ".join(flags))
    # remove genes as necessary
    for amira_allele in alleles_to_delete:
        del supplemented_clusters_of_interest[amira_allele]
        result_df = result_df[result_df["Amira allele"] != amira_allele]
    # add the comments as a column
    result_df["Comments"] = comments
    return result_df


def output_component_fastqs(output_dir, graph, fastq_content):
    if not os.path.exists(os.path.join(output_dir, "component_fastqs")):
        os.mkdir(os.path.join(output_dir, "component_fastqs"))
    for component in graph.components():
        node_hashes_in_component = [n.__hash__() for n in graph.get_nodes_in_component(component)]
        reads_in_component = graph.collect_reads_in_path(node_hashes_in_component)
        component_fastq_data = {r: fastq_content[r] for r in reads_in_component}
        write_fastq(
            os.path.join(output_dir, "component_fastqs", f"{component}.fastq.gz"),
            component_fastq_data,
        )


def write_reads_per_AMR_gene(output_dir, supplemented_clusters_of_interest):
    final_clusters_of_interest = {}
    for allele in supplemented_clusters_of_interest:
        # get the gene name with the allele name appended
        if os.path.exists(
            os.path.join(output_dir, "AMR_allele_fastqs", allele, "06.final_sequence.fasta")
        ):
            with open(
                os.path.join(output_dir, "AMR_allele_fastqs", allele, "06.final_sequence.fasta")
            ) as i:
                reference_allele_name = i.read().split(" ")[0].replace(">", "")
        else:
            with open(
                os.path.join(output_dir, "AMR_allele_fastqs", allele, "03.sequence_to_polish.fasta")
            ) as i:
                reference_allele_name = i.read().split(" ")[0].replace(">", "")
        if "\n" in reference_allele_name:
            reference_allele_name = reference_allele_name.split("\n")[0]
        new_name = f"{allele};{reference_allele_name}"
        new_reads = set()
        for r in supplemented_clusters_of_interest[allele]:
            new_reads.add(r.split("_")[0])
        final_clusters_of_interest[new_name] = list(new_reads)
    with open(os.path.join(output_dir, "reads_per_amr_gene.json"), "w") as o:
        o.write(json.dumps(final_clusters_of_interest))


def write_fasta(file_path, sequences):
    with open(file_path, "w") as outFasta:
        outFasta.write("\n".join(sequences))


def run_subprocess(command):
    subprocess.run(command, shell=True, check=True)


def map_reads(
    output_dir,
    reference_fasta,
    input_fasta,
    output_prefix,
    minimap_opts,
    minimap2_path,
    samtools_path,
):
    sam_file = os.path.join(output_dir, f"{output_prefix}.mapped.sam")
    bam_file = os.path.join(output_dir, f"{output_prefix}.mapped.bam")
    map_command = f"{minimap2_path} {minimap_opts} {reference_fasta} {input_fasta} > {sam_file}"
    run_subprocess(map_command)
    sort_and_index_command = (
        f"{samtools_path} sort {sam_file} > {bam_file} && {samtools_path} index {bam_file}"
    )
    run_subprocess(sort_and_index_command)
    return bam_file


def prepare_reference_sequences(gene_name, reference_alleles, output_dir):
    sequences = [f">{allele}\n{seq}" for allele, seq in reference_alleles[gene_name].items()]
    write_fasta(os.path.join(output_dir, "01.reference_alleles.fasta"), sequences)
    return [len(seq) for seq in reference_alleles[gene_name].values()]


def racon_polish(
    output_dir,
    racon_path,
    input_fasta,
    sam_file,
    reference_fasta,
    output_fasta,
    window_length,
):
    racon_command = f"{racon_path} -u -t 1 --no-trimming -w {window_length} {input_fasta} "
    racon_command += f"{sam_file} {reference_fasta} > {output_fasta}"
    run_subprocess(racon_command)


def racon_one_iteration(
    output_dir,
    racon_path,
    read_file,
    sam_file,
    sequence_to_polish,
    polished_sequence,
    window_size,
    minimap2_path,
    samtools_path,
):
    # run minimap2
    bam_file = map_reads(
        output_dir,
        os.path.join(output_dir, sequence_to_polish),
        read_file,
        sam_file.replace(".mapped.sam", ""),
        "-a --MD -t 1 -x map-ont --eqx",
        minimap2_path,
        samtools_path,
    )
    # run racon
    racon_polish(
        output_dir,
        racon_path,
        read_file,
        os.path.join(output_dir, sam_file),
        os.path.join(output_dir, sequence_to_polish),
        os.path.join(output_dir, polished_sequence),
        window_size,
    )
    # rename the polished sequence
    shutil.copy(
        os.path.join(output_dir, polished_sequence),
        os.path.join(output_dir, sequence_to_polish),
    )
    return bam_file


def create_output_dir(output_dir, allele_name):
    output_dir = os.path.join(output_dir, allele_name)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    return output_dir


def get_closest_allele(
    bam_file_path, mapping_type, required_identity, required_coverage, ref_cov_proportion=None
):
    valid_references = []
    invalid_references = []
    ref_covered = {}
    ref_matching = {}
    ref_lengths = {}
    ref_cigarstrings = {}
    ref_cigartuples = {}
    unique_reads = set()
    # Open the BAM file
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        for read in bam_file.fetch():
            # check if the read is mapped
            if read.is_unmapped:
                continue
            # add the read name to a set
            unique_reads.add(read.query_name)
            # get the length of the reference
            total_length = bam_file.get_reference_length(read.reference_name)
            # add the reference to the coverage dictionary
            if read.reference_name not in ref_covered:
                ref_covered[read.reference_name] = 0
                ref_matching[read.reference_name] = 0
                ref_lengths[read.reference_name] = total_length
            # get the proportion of bases matching the reference
            matching_bases = 0
            for op, length in read.cigartuples:
                if op == 7:
                    matching_bases += length
            if mapping_type == "reads":
                prop_matching = matching_bases / total_length  # read.query_alignment_length
                prop_covered = ref_cov_proportion[read.reference_name]
            if mapping_type == "allele":
                prop_matching = matching_bases / read.infer_read_length()
                # get the proportion of the reference covered by the query
                prop_covered = read.query_alignment_length / total_length
            # add the stats to the dictionaries
            if prop_matching > ref_matching[read.reference_name]:
                ref_matching[read.reference_name] = prop_matching
                ref_cigarstrings[read.reference_name] = read.cigarstring
                ref_cigartuples[read.reference_name] = read.cigartuples
            if prop_covered > ref_covered[read.reference_name]:
                ref_covered[read.reference_name] = prop_covered
    for ref in ref_matching:
        if ref_covered[ref] >= required_coverage - 0.05:
            valid_references.append(
                (
                    ref,
                    ref_matching[ref],
                    ref_lengths[ref],
                    ref_covered[ref],
                    ref_cigarstrings[ref],
                    ref_cigartuples[ref],
                )
            )
        else:
            invalid_references.append(
                (
                    ref,
                    ref_matching[ref],
                    ref_lengths[ref],
                    ref_covered[ref],
                    ref_cigarstrings[ref],
                    ref_cigartuples[ref],
                )
            )
    valid_references = sorted(
        valid_references, key=lambda x: (min(1, x[3]), x[1], x[2]), reverse=True
    )
    if len(valid_references) != 0:
        return True, valid_references, unique_reads
    else:
        invalid_references = sorted(invalid_references, key=lambda x: (x[3], x[1]), reverse=True)
        return False, invalid_references, unique_reads


def get_ref_allele_pileups(bam_file, output_dir):
    read_depths = []
    ref_allele_positions = {}
    cov_proportion = {}
    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Iterate over each reference (chromosome/contig)
        for ref in bam.references:
            ref_length = bam.get_reference_length(ref)
            # Initialize a list to store read depths for each position
            depth_list = [0] * ref_length
            # Fetch all reads mapped to this reference
            for read in bam.fetch(ref):
                if read.is_unmapped:
                    continue
                # Iterate over the aligned positions of the read
                for pos in read.get_reference_positions():
                    if pos < ref_length:
                        depth_list[pos] += 1
            # Convert the depth list to a string format
            depth_list_str = ",".join(map(str, depth_list))
            read_depths.append(f">{ref}\n{depth_list_str}")
            # get the first and last non-zero element
            try:
                first_index = depth_list.index(next(x for x in depth_list if x != 0))
                last_index = (
                    len(depth_list)
                    - 1
                    - depth_list[::-1].index(next(x for x in reversed(depth_list) if x != 0))
                )
            except StopIteration:
                first_index, last_index = None, None
            ref_allele_positions[ref] = (first_index, last_index)
            cov_proportion[ref] = len([d for d in depth_list if d != 0]) / len(depth_list)
    # Write the read depths to the output file
    output_file_path = os.path.join(output_dir, "reference_allele_coverages.txt")
    with open(output_file_path, "w") as output_file:
        output_file.write("\n".join(read_depths))
    return ref_allele_positions, cov_proportion


def get_allele_sequence(gene_name, allele, reference_genes):
    return reference_genes[gene_name][allele]


def compare_reads_to_references(
    allele_file,
    output_dir,
    reference_genes,
    racon_path,
    fastqContent,
    phenotypes,
    debug,
    minimap2_path,
    required_identity,
    required_coverage,
    samtools_path,
):
    allele_name = os.path.basename(allele_file).replace(".fastq.gz", "")
    gene_name = "_".join(allele_name.split("_")[:-1])
    output_dir = create_output_dir(output_dir, allele_name)
    prepare_reference_sequences(gene_name, reference_genes, output_dir)
    bam_file = map_reads(
        output_dir,
        os.path.join(output_dir, "01.reference_alleles.fasta"),
        allele_file,
        "02.read",
        "-a --MD -t 1 -x map-ont --eqx",
        minimap2_path,
        samtools_path,
    )
    ref_allele_positions, ref_cov_proportion = get_ref_allele_pileups(bam_file, output_dir)
    if debug is False:
        os.remove(os.path.join(output_dir, "reference_allele_coverages.txt"))
    validity, references, unique_reads = get_closest_allele(
        bam_file, "reads", required_identity, required_coverage, ref_cov_proportion
    )
    if validity is True:
        (
            valid_allele,
            match_proportion,
            match_length,
            coverage_proportion,
            cigarstring,
            cigartuple,
        ) = references[0]
        valid_allele_sequence = get_allele_sequence(gene_name, valid_allele, reference_genes)
        first_base, last_base = ref_allele_positions[valid_allele]
        write_fasta(
            os.path.join(output_dir, "03.sequence_to_polish.fasta"),
            [f">{valid_allele}\n{valid_allele_sequence[first_base: last_base+1]}"],
        )
        iterations = 5
        for _ in range(iterations):
            try:
                racon_one_iteration(
                    output_dir,
                    racon_path,
                    allele_file,
                    "02.read.mapped.sam",
                    "03.sequence_to_polish.fasta",
                    "04.polished_sequence.fasta",
                    str(len(valid_allele_sequence) + 200),
                    minimap2_path,
                    samtools_path,
                )
            except subprocess.CalledProcessError:
                try:
                    gene_name = valid_allele.split(".")[0]
                    closest_ref = valid_allele.split(".")[1]
                except IndexError:
                    gene_name = "_".join(allele_name.split("_")[:-1])
                    closest_ref = valid_allele
                # report the phenotype if there is one
                if valid_allele in phenotypes:
                    phenotype = phenotypes[valid_allele]
                else:
                    phenotype = ""
                # return the result
                return {
                    "Determinant name": gene_name,
                    "Sequence name": phenotype,
                    "Closest reference": closest_ref,
                    "Reference length": match_length,
                    "Identity (%)": 0,
                    "Coverage (%)": 0,
                    "Cigar string": "",
                    "Amira allele": allele_name,
                    "Number of reads used for polishing": len(unique_reads),
                }
        bam_file = map_reads(
            output_dir,
            os.path.join(output_dir, "01.reference_alleles.fasta"),
            os.path.join(output_dir, "04.polished_sequence.fasta"),
            "05.read",
            "-a --MD -t 1 --eqx",
            minimap2_path,
            samtools_path,
        )
        validity, references, _ = get_closest_allele(
            bam_file, "allele", required_identity, required_coverage
        )
        max_similarity = references[0][1]
        references = [r for r in references if r[1] == max_similarity]
        if len(references) == 1:
            (
                closest_allele,
                match_proportion,
                match_length,
                coverage_proportion,
                cigarstring,
                cigartuple,
            ) = references[0]
            with open(os.path.join(output_dir, "04.polished_sequence.fasta")) as i:
                final_allele_sequence = "".join(i.read().split("\n")[1:])
            write_fasta(
                os.path.join(output_dir, "06.final_sequence.fasta"),
                [f">{closest_allele}\n{final_allele_sequence}"],
            )
            try:
                gene_name = closest_allele.split(".")[0]
                closest_ref = closest_allele.split(".")[1]
            except IndexError:
                gene_name = "_".join(allele_name.split("_")[:-1])
                closest_ref = closest_allele
            # report the phenotype if there is one
            if closest_allele in phenotypes:
                phenotype = phenotypes[closest_allele]
            else:
                phenotype = ""
            # calculate the identity, ignoring soft- & hard-clipped bases
            matching = 0
            total = 0
            for op, length in cigartuple:
                if op == 7:
                    matching += length
                if op != 4 and op != 5:
                    total += length
            identity = matching / total
            # return the result
            return {
                "Determinant name": gene_name,
                "Sequence name": phenotype,
                "Closest reference": closest_ref,
                "Reference length": match_length,
                "Identity (%)": round(identity * 100, 1),
                "Coverage (%)": min(100.0, round(coverage_proportion * 100, 1)),
                "Cigar string": cigarstring,
                "Amira allele": allele_name,
                "Number of reads used for polishing": len(unique_reads),
            }
        if len(references) > 1:
            (closest_allele, match_proportion, match_length, coverage_proportion, cigarstrings) = (
                [],
                [],
                [],
                [],
                [],
            )
            for r in references:
                closest_allele.append(r[0])
                match_length.append(r[2])
                coverage_proportion.append(r[3])
                cigarstrings.append(r[4])
                # calculate the identity, ignoring soft- & hard-clipped bases
                matching = 0
                total = 0
                for op, length in r[5]:
                    if op == 7:
                        matching += length
                    if op != 4 and op != 5:
                        total += length
                identity = matching / total
                match_proportion.append(identity)
            with open(os.path.join(output_dir, "04.polished_sequence.fasta")) as i:
                final_allele_sequence = "".join(i.read().split("\n")[1:])
            write_fasta(
                os.path.join(output_dir, "06.final_sequence.fasta"),
                [f">{'/'.join(closest_allele)}\n{final_allele_sequence}"],
            )
            try:
                gene_names = "/".join(sorted(list(set([c.split(".")[0] for c in closest_allele]))))
                closest_refs = "/".join([c.split(".")[1] for c in closest_allele])
                phenotypes = "/".join(
                    [phenotypes[c] if c in phenotypes else "" for c in closest_allele]
                )
            except IndexError:
                gene_names = "_".join(allele_name.split("_")[:-1])
                closest_refs = "/".join(closest_allele)
                phenotypes = "/".join(
                    [phenotypes[c] if c in phenotypes else "" for c in closest_allele]
                )
            return {
                "Determinant name": gene_names,
                "Sequence name": phenotypes,
                "Closest reference": closest_refs,
                "Reference length": "/".join([str(m) for m in match_length]),
                "Identity (%)": "/".join([str(round(p * 100, 1)) for p in match_proportion]),
                "Coverage (%)": "/".join(
                    [str(min(100.0, round(p * 100, 1))) for p in coverage_proportion]
                ),
                "Cigar string": "/".join(cigarstrings),
                "Amira allele": allele_name,
                "Number of reads used for polishing": len(unique_reads),
            }
    else:
        if len(references) != 0:
            (
                invalid_allele,
                match_proportion,
                match_length,
                coverage_proportion,
                cigarstring,
                cigartuple,
            ) = references[0]
            try:
                gene_name = invalid_allele.split(".")[0]
                closest_ref = invalid_allele.split(".")[1]
            except IndexError:
                gene_name = "_".join(allele_name.split("_")[:-1])
                closest_ref = invalid_allele
            # report the phenotype if there is one
            if invalid_allele in phenotypes:
                phenotype = phenotypes[invalid_allele]
            else:
                phenotype = ""
            # calculate the identity, ignoring soft- & hard-clipped bases
            matching = 0
            total = 0
            for op, length in cigartuple:
                if op == 7:
                    matching += length
                if op != 4 and op != 5:
                    total += length
            identity = matching / total
            # return the results
            return {
                "Determinant name": gene_name,
                "Sequence name": phenotype,
                "Closest reference": closest_ref,
                "Reference length": match_length,
                "Identity (%)": round(identity * 100, 1),
                "Coverage (%)": min(100.0, round(coverage_proportion * 100, 1)),
                "Cigar string": cigarstring,
                "Amira allele": allele_name,
                "Number of reads used for polishing": len(unique_reads),
            }
        else:
            return {
                "Determinant name": "",
                "Sequence name": "",
                "Closest reference": "",
                "Reference length": 0,
                "Identity (%)": 0,
                "Coverage (%)": 0,
                "Cigar string": "",
                "Amira allele": allele_name,
                "Number of reads used for polishing": len(unique_reads),
            }


def get_alleles(
    readFiles,
    racon_path,
    threads,
    output_dir,
    reference_genes,
    fastqContent,
    phenotypes_path,
    debug,
    minimap2_path,
    required_identity,
    required_coverage,
    samtools_path,
):
    # import the phenotypes
    with open(phenotypes_path) as i:
        phenotypes = json.load(i)
    # batch the read files for multi processing
    job_list = [readFiles[i : i + threads] for i in range(0, len(readFiles), threads)]
    gene_data = []
    for subset in tqdm(job_list):
        gene_data += Parallel(n_jobs=threads)(
            delayed(compare_reads_to_references)(
                r,
                output_dir,
                reference_genes,
                racon_path,
                fastqContent.get("_".join(os.path.basename(r).split("_")[:-1]), None),
                phenotypes,
                debug,
                minimap2_path,
                required_identity,
                required_coverage,
                samtools_path,
            )
            for r in subset
        )
    return pd.DataFrame(gene_data)


def genotype_promoters(
    result_df,
    reference_alleles,
    output_dir,
    minimap2_path,
    racon_path,
    phenotypes,
    debug,
    output_components,
    samtools_path,
):
    # skip if there are no promoters in the reference
    if not any("_promoter" in a for a in reference_alleles):
        sys.stderr.write("\nAmira: No promoters found in reference FASTA.\n")
        return result_df
    # iterate through the rows in the amira output
    for index, row in result_df.iterrows():
        # get the amira gene name for this row
        amira_gene = "_".join(row["Amira allele"].split("_")[:-1])
        # get the hypothetical promoter name
        promoter_name = amira_gene + "_promoter"
        # check if the promoter name exists in the reference
        if promoter_name in reference_alleles:
            # get the allele index
            gene_index = row["Amira allele"].split("_")[-1]
            # get the new allele name
            promoter_allele_name = f"{promoter_name}_{gene_index}"
            # make the output directory
            output_allele_dir = create_output_dir(output_dir, promoter_allele_name)
            # copy the fastq
            allele_file = os.path.join(output_allele_dir, f"{promoter_name}_{gene_index}.fastq.gz")
            shutil.copyfile(
                os.path.join(
                    output_dir,
                    row["Amira allele"],
                    f"{row['Amira allele']}.fastq.gz",
                ),
                allele_file,
            )
            # get the closest promoter allele
            closest_reference = compare_reads_to_references(
                allele_file,
                output_dir,
                reference_alleles,
                racon_path,
                {},
                phenotypes,
                debug,
                minimap2_path,
                samtools_path,
            )
            # skip if no allele was similar
            if not os.path.exists(
                os.path.join(output_dir, promoter_allele_name, "06.final_sequence.fasta")
            ):
                continue
            # import the final bam
            SNPs_present = {
                "Determinant name": [],
                "Sequence name": [],
                "Closest reference": [],
                "Reference length": [],
                "Identity (%)": [],
                "Coverage (%)": [],
                "Cigar string": [],
                "Amira allele": [],
                "Number of reads used for polishing": [],
                "Approximate cellular copy number": [],
            }
            if output_components is True:
                SNPs_present["Component ID"] = []
            # Load reference sequence
            ref_sequences = pysam.FastaFile(
                os.path.join(output_dir, promoter_allele_name, "01.reference_alleles.fasta")
            )
            read_sequence = pysam.FastaFile(
                os.path.join(output_dir, promoter_allele_name, "06.final_sequence.fasta")
            )
            bam = pysam.AlignmentFile(
                os.path.join(output_dir, promoter_allele_name, "05.read.mapped.sam"), "rb"
            )
            # skip if the allele is the same as the reference
            if not closest_reference["Identity (%)"] < 100:
                continue
            # Store changes per read
            changes = {}
            for read in bam.fetch():
                if read.is_unmapped:
                    continue

                read_changes = []
                read_seq = read_sequence.fetch(read.query_name)
                ref_positions = read.get_reference_positions(full_length=True)
                ref_seq = ref_sequences.fetch(read.reference_name)
                read_index = 0

                for cigar_op, length in read.cigartuples:
                    if cigar_op == 8:  # mismatch
                        for i in range(length):
                            ref_pos = ref_positions[read_index + i]
                            if ref_pos is not None:  # Handle potential None in ref_positions
                                ref_base = ref_seq[ref_pos].upper()
                                read_base = read_seq[read_index + i].upper()
                                read_changes.append(f"{ref_base}{ref_pos + 1}{read_base}")
                        read_index += length

                    elif cigar_op == 1:  # insertion
                        # Record the insertion as occurring after the last reference position
                        insertion_seq = read_seq[read_index : read_index + length].upper()
                        last_ref_pos = ref_positions[read_index - 1] if read_index > 0 else None
                        if last_ref_pos is not None:
                            read_changes.append(f"{last_ref_pos + 1}I{insertion_seq}")
                        read_index += length

                    elif cigar_op == 2:  # deletion
                        # Record the deletion range in the reference
                        del_start = ref_positions[read_index - 1] + 1 if read_index > 0 else None
                        del_end = ref_positions[read_index + length - 1]
                        if del_start is not None and del_end is not None:
                            del_seq = ref_seq[
                                del_start - 1 : del_end
                            ].upper()  # Fetch the deleted sequence from reference
                            read_changes.append(f"{del_start}-{del_end}D{del_seq}")
                    else:
                        read_index += (
                            length  # Increment for other operations (e.g., matches, soft clips)
                        )
                changes[read.reference_name] = read_changes
            # split up the gene name
            for ref in changes:
                gene_name = ref.split(".")[0] + "_promoter_" + "_".join(changes[ref])
                accession = ".".join(ref.split(".")[0:2])
                # get the identity
                prop_matching = closest_reference["Identity (%)"]
                # get the coverage
                prop_covered = closest_reference["Coverage (%)"]
                # report the phenotype if there is one
                if ref in phenotypes:
                    phenotype = phenotypes[read.reference_name]
                else:
                    phenotype = ""
                # add the information
                SNPs_present["Determinant name"].append(gene_name)
                SNPs_present["Sequence name"].append(phenotype)
                SNPs_present["Closest reference"].append(accession)
                SNPs_present["Reference length"].append(closest_reference["Reference length"])
                SNPs_present["Identity (%)"].append(prop_matching)
                SNPs_present["Coverage (%)"].append(prop_covered)
                SNPs_present["Cigar string"].append(closest_reference["Cigar string"])
                SNPs_present["Amira allele"].append(promoter_allele_name)
                SNPs_present["Number of reads used for polishing"].append(
                    closest_reference["Number of reads used for polishing"]
                )
                SNPs_present["Approximate cellular copy number"].append(row["Approximate cellular copy number"])
                if output_components is True:
                    SNPs_present["Component ID"].append(row["Component ID"])
            # Close the files
            bam.close()
            ref_sequences.close()
            read_sequence.close()
            if len(SNPs_present["Determinant name"]) > 0:
                new_row = pd.DataFrame(SNPs_present)
                result_df = pd.concat([result_df, new_row], ignore_index=True)
    return result_df


def get_mean_read_depth_per_contig(bam_file, samtools_path, core_genes=None):
    # Run samtools depth command and capture output
    command = [
        samtools_path,
        "mpileup",
        "-aa",
        "-Q",
        "0",
        "--ff",
        "UNMAP,QCFAIL,DUP,SUPPLEMENTARY",
        bam_file,
    ]
    # Run the command and capture the output
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    # Parse the mpileup output
    contig_depths = {}
    contig_positions = {}
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) > 3:  # Ensure sufficient fields exist
            contig = fields[0]  # Contig name is the 1st column
            if core_genes is not None:
                if contig not in core_genes:
                    continue
            depth = int(fields[3])  # 4th column is the depth
            if contig not in contig_depths:
                contig_depths[contig] = 0
                contig_positions[contig] = 0
            contig_depths[contig] += depth
            contig_positions[contig] += 1
    # Calculate mean depth for each contig
    mean_depths = {
        contig: (contig_depths[contig] / contig_positions[contig]) for contig in contig_depths
    }
    return mean_depths


def kmer_cutoff_estimation(kmer_counts):
    i_values = np.array(list(kmer_counts.keys()))
    xi_values = np.array(list(kmer_counts.values()))

    # Define the likelihood function
    def neg_log_likelihood(params):
        w, c = params  # w = error proportion, c = coverage mean
        if w < 0 or w > 1 or c <= 0:
            return np.inf
        # compute Poisson probabilities
        error_prob = poisson.pmf(i_values, mu=1)
        real_prob = poisson.pmf(i_values, mu=c)
        # mixture model
        mix_prob = w * error_prob + (1 - w) * real_prob
        mix_prob[mix_prob == 0] = 1e-10
        # compute negative log-likelihood
        return -np.sum(xi_values * np.log(mix_prob))

    # initial parameters
    init_params = [0.1, 10]
    # optimise parameters using BFGS
    result = minimize(
        neg_log_likelihood, init_params, method="L-BFGS-B", bounds=[(0, 1), (1e-5, None)]
    )
    # get optimised parameters
    w_opt, c_opt = result.x
    for i in i_values:
        error_prob = poisson.pmf(i, mu=1) * w_opt
        real_prob = poisson.pmf(i, mu=c_opt) * (1 - w_opt)
        if real_prob > error_prob:
            return i
    return 0


def estimate_kmer_depth(kmer_counts, histogram_path, debug=False):
    x_values, y_values = zip(*sorted(kmer_counts.items()))
    log_counts = np.log(np.array(y_values) + 1)
    window_length = min(30, len(log_counts) // 2 * 2 + 1)
    smoothed_log_counts = savgol_filter(log_counts, window_length, 3)
    peak_indices, _ = find_peaks(smoothed_log_counts)
    # Get the highest peak
    max_peak_index = peak_indices[np.argmax(smoothed_log_counts[peak_indices])]
    depth = x_values[max_peak_index]  # Get k-mer depth based on x-values
    if debug is True:
        plt.bar(x_values, log_counts)
        plt.plot(x_values, smoothed_log_counts, color="red", label="Smoothed counts")
        plt.axvline(depth, color="red")
        plt.xlim(0, 500)
        plt.savefig(histogram_path.replace(".filtered.histo", ".png"), dpi=600)
    return depth


def import_jellyfish_histo(histo_path):
    with open(histo_path, "r") as f:
        lines = f.read().split("\n")
    kmer_counts = {}
    for line in lines:
        if line != "":
            val = int(line.split(" ")[0])
            count = int(line.split(" ")[1])
            kmer_counts[val] = count
    return kmer_counts


def load_kmer_counts(counts_file):
    # load the counts
    abundances = []
    with open(counts_file) as i:
        for line in i:
            parts = line.split()
            if len(parts) == 2:  # Ensure line has both k-mer and count
                count = int(parts[1])
                if count != 0:
                    abundances.append(count)
    return abundances


def estimate_overall_read_depth(full_reads, k, threads, debug):
    # count the k-mers in the full read set
    jf_out = full_reads.replace(".fastq.gz", ".jf")
    histo_out = full_reads.replace(".fastq.gz", ".histo")
    sys.stderr.write("\nAmira: counting k-mers using Jellyfish.\n")
    command = (
        f"bash -c 'jellyfish count -m {k} -s 20M -t {threads} -C <(zcat {full_reads}) -o {jf_out}'"
    )
    command += f" && jellyfish histo {jf_out} > {histo_out}"
    subprocess.run(command, shell=True, check=True)
    # import the counts
    kmer_counts = import_jellyfish_histo(histo_out)
    # estimate the kmer cutoff
    cutoff = kmer_cutoff_estimation(kmer_counts)
    sys.stderr.write(f"\nAmira: filtering k-mers with count below cutoff ({cutoff}).\n")
    jf_out = full_reads.replace(".fastq.gz", ".filtered.jf")
    histo_out = full_reads.replace(".fastq.gz", ".filtered.histo")
    kmers_out = full_reads.replace(".fastq.gz", ".filtered.kmers.txt")
    command = f"bash -c 'jellyfish count -m {k} -s 20M -L {cutoff} "
    command += f"-t {threads} -C <(zcat {full_reads}) -o {jf_out}'"
    command += f" && jellyfish dump -c {jf_out} > {kmers_out}"
    command += f" && jellyfish histo {jf_out} > {histo_out}"
    subprocess.run(command, shell=True, check=True)
    # import the counts
    filtered_kmer_counts = import_jellyfish_histo(histo_out)
    # estimate the read depth from the kmerss
    kmer_depth = estimate_kmer_depth(filtered_kmer_counts, histo_out, debug)
    return kmer_depth, jf_out


def estimate_depth(counts_file):
    # load the counts
    abundances = load_kmer_counts(counts_file)
    return statistics.median(abundances)


def estimate_copy_numbers(
    fastq_content,
    path_reads,
    amira_alleles,
    fastq_file,
    threads,
    samtools_path,
    raw_read_depth,
    debug,
):
    # make the output directory
    outdir = os.path.join(os.path.dirname(fastq_file), "AMR_allele_fastqs", "path_reads")
    # write out the reads for each path
    sys.stderr.write("\nAmira: writing out the reads for each path.\n")
    path_mapping = {}
    path_list = list(path_reads.keys())
    path_fastqs = []
    for i in range(len(path_list)):
        # asign an id to each path
        path_id = i + 1
        path = path_list[i]
        path_mapping[path_id] = list(path)
        # write out the reads of each path
        path_fastqs.append(write_path_fastq(path_reads[path], fastq_content, outdir, path_id))
    # write out the mapping
    with open(os.path.join(outdir, "path_id_mapping.json"), "w") as o:
        o.write(json.dumps(path_mapping))
    # define k
    k = 15
    # estimate the overall depth
    read_depth, full_jf_out = estimate_overall_read_depth(fastq_file, k, threads, debug)
    sys.stderr.write(f"\nAmira: estimated k-mer depth = {read_depth}.\n")
    # get the genomic copy number of each gene
    gene_counts = {}
    for i in path_mapping:
        gene_counts[i] = {}
        for g in path_mapping[i]:
            strandless = g[1:]
            if strandless in amira_alleles:
                gene = "_".join(strandless.split("_")[:-1])
                gene_counts[i][gene] = gene_counts[i].get(gene, 0) + 1
    # estimate cellular cpy numbers for each subset fastq
    normalised_depths, mean_depth_per_reference = {}, {}
    for locus_reads in tqdm(path_fastqs):
        # query the kmers
        counts_file = locus_reads.replace(".fastq.gz", ".kmer_counts.txt")
        command = (
            f"bash -c 'jellyfish query -s <(zcat {locus_reads}) {full_jf_out} > {counts_file}'"
        )
        subprocess.run(command, shell=True, check=True)
        # get the counts of the requested kmers
        depth_estimate = estimate_depth(counts_file)
        # get the id of the path
        path_id = int(os.path.basename(locus_reads).replace(".fastq.gz", ""))
        # get the allele
        full_path = path_mapping[path_id]
        for g in full_path:
            allele_name = g[1:]
            if allele_name not in amira_alleles:
                continue
            # get the gene
            gene = "_".join(allele_name.split("_")[:-1])
            # store the depth estimate
            normalised_depths[allele_name] = depth_estimate / (
                read_depth * gene_counts[path_id][gene]
            )
            mean_depth_per_reference[allele_name] = depth_estimate
    return normalised_depths, mean_depth_per_reference


def write_fastqs_for_genes_with_short_reads(
    clusters_to_add,
    overall_mean_node_coverage,
    longest_reads_for_genes,
    output_dir,
    files_to_assemble,
    fastq_content,
    supplemented_clusters_of_interest,
    allele_component_mapping,
):
    for allele in clusters_to_add:
        files_to_assemble.append(
            write_allele_fastq(clusters_to_add[allele], fastq_content, output_dir, allele)
        )
        supplemented_clusters_of_interest[allele] = clusters_to_add[allele]
        allele_component_mapping[allele] = None
        # iterate through the reads
        longest_read = None
        longest_read_length = 0
        for read in clusters_to_add[allele]:
            read_name = "_".join(read.split("_")[:-2])
            read_length = len(fastq_content[read_name]["sequence"])
            if read_length > longest_read_length:
                longest_read_length = read_length
                longest_read = read_name
        longest_reads_for_genes.append(f">{allele}\n{fastq_content[longest_read]['sequence']}")
    return longest_reads_for_genes, files_to_assemble


def write_fastqs_for_genes(
    clusters_of_interest, overall_mean_node_coverage, fastq_content, output_dir
):
    # track the longest read for each AMR gene copy
    longest_reads_for_genes = []
    supplemented_clusters_of_interest = {}
    allele_component_mapping = {}
    files_to_assemble = []
    for component in tqdm(clusters_of_interest):
        for gene in clusters_of_interest[component]:
            for allele in clusters_of_interest[component][gene]:
                files_to_assemble.append(
                    write_allele_fastq(
                        clusters_of_interest[component][gene][allele],
                        fastq_content,
                        output_dir,
                        allele,
                    )
                )
                supplemented_clusters_of_interest[allele] = clusters_of_interest[component][gene][
                    allele
                ]
                # store the component of the allele
                allele_component_mapping[allele] = component
                # iterate through the reads
                longest_read = None
                longest_read_length = 0
                for read in clusters_of_interest[component][gene][allele]:
                    read_name = "_".join(read.split("_")[:-2])
                    read_length = len(fastq_content[read_name]["sequence"])
                    if read_length > longest_read_length:
                        longest_read_length = read_length
                        longest_read = read_name
                longest_reads_for_genes.append(
                    f">{allele}\n{fastq_content[longest_read]['sequence']}"
                )
    return (
        longest_reads_for_genes,
        supplemented_clusters_of_interest,
        allele_component_mapping,
        files_to_assemble,
    )


def write_empty_result(output_dir):
    results = "Determinant name\tSequence name\tClosest reference\tReference length\t"
    results += "Identity (%)\tCoverage (%)\tAmira allele\t"
    results += "Number of reads used for polishing\tApproximate cellular copy number\n"
    with open(os.path.join(output_dir, "amira_results.tsv"), "w") as o:
        o.write(results)


def supplement_result_df(
    result_df, copy_numbers, mean_depth_per_reference, longest_read_lengths, debug
):
    estimates = []
    copy_depths = []
    read_lengths = []
    for index, row in result_df.iterrows():
        estimates.append(copy_numbers[row["Amira allele"]])
        copy_depths.append(mean_depth_per_reference[row["Amira allele"]])
        read_lengths.append(longest_read_lengths[row["Amira allele"]])
    result_df["Mean read depth"] = copy_depths
    result_df["Approximate cellular copy number"] = estimates
    if debug:
        result_df["Longest read length"] = read_lengths
    return result_df


def write_pandora_gene_calls(output_dir, gene_position_dict, annotatedReads, outfile_1, outfile_2):
    with open(outfile_1, "w") as o:
        o.write(json.dumps(gene_position_dict))
    with open(outfile_2, "w") as o:
        o.write(json.dumps(annotatedReads))
