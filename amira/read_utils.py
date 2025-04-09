import gzip
import os
import random
import re
import subprocess

import matplotlib.pyplot as plt


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
    else:
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
        lines = []
        # Iterate over the data
        for identifier, value in data.items():
            lines.append(f"@{identifier}\n")
            lines.append(f"{value['sequence']}\n")
            lines.append("+\n")
            lines.append(f"{value['quality']}\n")
        # write out the data
        fh.writelines(lines)


def downsample_reads(fastq_content, read_path, output_dir, max_reads=100000):
    # If no downsampling is needed, return original annotatedReads
    total_reads = len(fastq_content)
    if total_reads <= max_reads:
        selected_reads = list(fastq_content.keys())
    else:
        # Convert the items to a list before sampling
        selected_reads = set(random.sample(list(fastq_content.keys()), max_reads))
        # sample the fastq dictionary
        fastq_content = {k: v for k, v in fastq_content.items() if k in selected_reads}
    # write out the list of reads
    with open(os.path.join(output_dir, "selected_reads.txt"), "w") as o:
        o.write("\n".join(list(selected_reads)))
    # filter the fastq
    read_fastq_path = os.path.join(output_dir, "subsampled_reads.fq.gz")
    command = f"fastaq filter --ids_file {os.path.join(output_dir, 'selected_reads.txt')} "
    command += f"{read_path} {read_fastq_path}"
    subprocess.run(command, shell=True, check=True)
    return read_fastq_path


def write_modified_fastq(fastq_content, read_path, output_dir):
    # remove underscores in read names
    fastq_content = {re.sub(r"[\W_]+", "", r): fastq_content[r] for r in fastq_content}
    # write the modified fastq data to the output directoy
    if ".gz" not in read_path:
        read_fastq_path = os.path.join(output_dir, os.path.splitext(os.path.basename(read_path))[0])
    else:
        read_fastq_path = os.path.join(
            output_dir, os.path.splitext(os.path.splitext(os.path.basename(read_path))[0])[0]
        )
    read_fastq_path += ".fastq.gz"
    write_fastq(
        read_fastq_path,
        fastq_content,
    )
    return read_fastq_path
