import json
import os
import shutil
import subprocess
import sys

import pandas as pd
import pysam
from joblib import Parallel, delayed
from tqdm import tqdm

from amira.read_utils import write_fastq


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
    cores,
    short_reads,
    short_read_gene_positions,
    overall_mean_node_coverage,
):
    # Step 1: Assign reads to genes
    clusters_of_interest = graph.assign_reads_to_genes(
        sample_genesOfInterest, cores, {}, overall_mean_node_coverage
    )
    # Step 2: Get the unique genes found in clusters of interest
    found_genes_of_interest = get_found_genes(clusters_of_interest)
    # Step 3: Add AMR alleles that come from short reads
    clusters_to_add = add_amr_alleles(
        short_reads, short_read_gene_positions, sample_genesOfInterest, found_genes_of_interest
    )
    # Return the processed results
    return (
        clusters_to_add,
        clusters_of_interest,
    )


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
    filter_contamination,
    supplemented_clusters_of_interest,
    annotatedReads,
    sample_genesOfInterest,
    required_identity,
    required_coverage,
):
    # remove genes that do not have sufficient mapping coverage
    alleles_to_delete = []
    comments = []
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
                if coverage < 90:
                    flags.append("Partially present gene.")
        # remove alleles where all of the reads just contain AMR genes
        reads = supplemented_clusters_of_interest[row["Amira allele"]]
        if all(
            all(g[1:] in sample_genesOfInterest for g in annotatedReads[r.split("_")[0]])
            for r in reads
        ):
            # check if filter contaminants is on
            if filter_contamination is True:
                alleles_to_delete.append(row["Amira allele"])
                message = f"\nAmira: allele {row['Amira allele']} removed "
                message += "due to suspected contamination.\n"
                sys.stderr.write(message)
                continue
            else:
                flags.append("Potential contaminant.")
        # collect the flags
        comments.append(" ".join(flags))
    # add the comments as a column
    result_df["Comments"] = comments
    # remove genes as necessary
    for amira_allele in alleles_to_delete:
        del supplemented_clusters_of_interest[amira_allele]
        result_df = result_df[result_df["Amira allele"] != amira_allele]
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
                prop_matching = (
                    matching_bases / read.query_alignment_length
                ) * 100  # total_length) * 100
                prop_covered = ref_cov_proportion[read.reference_name]
            if mapping_type == "allele":
                prop_matching = (matching_bases / read.infer_read_length()) * 100
                # get the proportion of the reference covered by the query
                prop_covered = read.query_alignment_length / total_length
            # add the stats to the dictionaries
            if prop_matching > ref_matching[read.reference_name]:
                ref_matching[read.reference_name] = prop_matching
                ref_cigarstrings[read.reference_name] = read.cigarstring
            if prop_covered > ref_covered[read.reference_name]:
                ref_covered[read.reference_name] = prop_covered
    for ref in ref_matching:
        if ref_covered[ref] >= required_coverage - 0.05:
            valid_references.append(
                (ref, ref_matching[ref], ref_lengths[ref], ref_covered[ref], ref_cigarstrings[ref])
            )
        else:
            invalid_references.append(
                (ref, ref_matching[ref], ref_lengths[ref], ref_covered[ref], ref_cigarstrings[ref])
            )
    valid_references = sorted(valid_references, key=lambda x: (x[1], x[2]), reverse=True)
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
        valid_allele, match_proportion, match_length, coverage_proportion, cigarstring = references[
            0
        ]
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
                    str(len(valid_allele_sequence)),
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
                    "Number of reads": len(unique_reads),
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
            closest_allele, match_proportion, match_length, coverage_proportion, cigarstring = (
                references[0]
            )
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
            # return the result
            return {
                "Determinant name": gene_name,
                "Sequence name": phenotype,
                "Closest reference": closest_ref,
                "Reference length": match_length,
                "Identity (%)": round(match_proportion, 1),
                "Coverage (%)": min(100.0, round(coverage_proportion * 100, 1)),
                "Cigar string": cigarstring,
                "Amira allele": allele_name,
                "Number of reads": len(unique_reads),
            }
        if len(references) > 1:
            closest_allele, match_proportion, match_length, coverage_proportion, cigarstrings = (
                [],
                [],
                [],
                [],
                [],
            )
            for r in references:
                closest_allele.append(r[0])
                match_proportion.append(r[1])
                match_length.append(r[2])
                coverage_proportion.append(r[3])
                cigarstrings.append(r[4])
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
                "Identity (%)": "/".join([str(round(p, 1)) for p in match_proportion]),
                "Coverage (%)": "/".join(
                    [str(min(100.0, round(p * 100, 1))) for p in coverage_proportion]
                ),
                "Cigar string": "/".join(cigarstrings),
                "Amira allele": allele_name,
                "Number of reads": len(unique_reads),
            }
    else:
        if len(references) != 0:
            invalid_allele, match_proportion, match_length, coverage_proportion, cigarstring = (
                references[0]
            )
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
            # return the results
            return {
                "Determinant name": gene_name,
                "Sequence name": phenotype,
                "Closest reference": closest_ref,
                "Reference length": match_length,
                "Identity (%)": round(match_proportion, 1),
                "Coverage (%)": min(100.0, round(coverage_proportion * 100, 1)),
                "Cigar string": cigarstring,
                "Amira allele": allele_name,
                "Number of reads": len(unique_reads),
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
                "Number of reads": len(unique_reads),
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
                "Number of reads": [],
                "Approximate copy number": [],
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
                SNPs_present["Number of reads"].append(closest_reference["Number of reads"])
                SNPs_present["Approximate copy number"].append(row["Approximate copy number"])
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
    command = [samtools_path, "mpileup", "-aa", "-Q", "0", "--ff", "UNMAP,QCFAIL,DUP", bam_file]
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


def estimate_copy_numbers(mean_read_depth, ref_file, fastq_file, threads, samtools_path):
    sam_file = ref_file.replace(".fasta", ".sam")
    bam_file = sam_file.replace(".sam", ".bam")
    command = f"minimap2 -x map-ont -t {threads} "
    command += f"-a --MD --eqx -o {sam_file} {ref_file} {fastq_file} "
    command += f"&& {samtools_path} sort {sam_file} > {bam_file}"
    subprocess.run(command, shell=True, check=True)
    # get mean depths across each reference
    mean_depth_per_reference = get_mean_read_depth_per_contig(bam_file, samtools_path)
    # normalise by the mean depth across core genes
    if mean_read_depth is None:
        mean_read_depth = min(mean_depth_per_reference.values())
    normalised_depths = {
        c: round(mean_depth_per_reference[c] / mean_read_depth, 2) for c in mean_depth_per_reference
    }
    os.remove(sam_file)
    os.remove(bam_file)
    return normalised_depths


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
        if len(clusters_to_add[allele]) > overall_mean_node_coverage / 20:
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
        else:
            message = f"\nAmira: allele {allele} removed "
            message += f"due to an insufficient number of reads ({len(clusters_to_add[allele])}).\n"
            sys.stderr.write(message)
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
                if (
                    len(clusters_of_interest[component][gene][allele])
                    > overall_mean_node_coverage / 20
                ):
                    files_to_assemble.append(
                        write_allele_fastq(
                            clusters_of_interest[component][gene][allele],
                            fastq_content,
                            output_dir,
                            allele,
                        )
                    )
                    supplemented_clusters_of_interest[allele] = clusters_of_interest[component][
                        gene
                    ][allele]
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
                else:
                    message = f"\nAmira: allele {allele} removed "
                    message += "due to an insufficient number of reads "
                    message += f"({len(clusters_of_interest[component][gene][allele])}).\n"
                    sys.stderr.write(message)
    return (
        longest_reads_for_genes,
        supplemented_clusters_of_interest,
        allele_component_mapping,
        files_to_assemble,
    )
