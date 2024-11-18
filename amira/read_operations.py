import gzip
import os
import random

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


def downsample_reads(annotatedReads, max_reads=100000):
    # If no downsampling is needed, return original annotatedReads
    total_reads = len(annotatedReads)
    if total_reads <= max_reads:
        return annotatedReads

    # Convert the items to a list before sampling
    sampled_items = random.sample(list(annotatedReads.items()), max_reads)
    return dict(sampled_items)
