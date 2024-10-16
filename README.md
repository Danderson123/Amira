# Amira

## Overview

Here are some brief instructions to run Amira on some test data. Amira runs on MacOS and Linux so if you have a windows computer, you will need to first install the Windows Subsystem for Linux from [here](https://learn.microsoft.com/en-us/windows/wsl/install) and run Amira using a Linux terminal.

## Installation

You first need to clone and then cd into the Amira repository by running:
```
git clone https://github.com/Danderson123/Amira && cd Amira
```
We use poetry to manage Amira's python dependencies, which can be installed by running:
```
pip install poetry
```
The python dependencies can then be installed by running:
```
poetry install
```
Amira also has two non-python dependencies, minimap2 and racon. You will need to build a binary for these tools and then give the path to the binary to Amira (see the example command).
* The instructions to build a binary for minimap2 are available from [here](https://github.com/lh3/minimap2).
* The instructions to build a binary for racon are available from [here](https://github.com/isovic/racon).

## Running Amira

Amira can be run directly on the output of Pandora, or it can be run using JSON files that indicate the list of genes on each sequencing read and the list of gene positions on each sequencing read.

To run Amira from the JSON files, you can use this command. You will need to replace `<PATH TO RACON BINARY>` with the absolute path to the racon binary you made earlier and replace `<PATH TO MINIMAP2 BINARY>` with the path to the minimap2 binary.
```
python3 amira/__main__.py --pandoraJSON test_data/gene_calls_with_gene_filtering.json --gene-positions test_data/gene_positions_with_gene_filtering.json --pandoraConsensus test_data/pandora.consensus.fq.gz --readfile test_data/SRR23044220_1.fastq.gz --output amira_output --gene-path test_data/AMR_gene_references.fa --phenotypes test_data/AMR_calls.json --racon-path <PATH TO RACON BINARY> --minimap2-path <PATH TO MINIMAP2 BINARY> --debug --cores 1 --sample-reads --filter-contaminants
```