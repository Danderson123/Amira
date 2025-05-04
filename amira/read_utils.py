import gzip
import os
import random
import subprocess
import sys

import matplotlib.pyplot as plt
import pysam


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


def parse_fastq(fastq_file):
    fastq_dict = {}
    with pysam.FastxFile(fastq_file) as fh:
        for entry in fh:
            fastq_dict[entry.name] = {"sequence": entry.sequence, "quality": entry.quality}
    return fastq_dict


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
    # compress the fastq file
    if ".gz" not in read_path:
        sys.stderr.write(f"\nAmira: compressing {read_path}\n")
        subprocess.run(f"gzip -1 {read_path}", shell=True, check=True)
        read_path = read_path + ".gz"
    return read_path, fastq_content
