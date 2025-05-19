# Amira

## Introduction

Amira is an AMR gene detection tool designed to work directly from bacterial long read sequences. Amira makes it easy to reliably identify the AMR genes in a bacterial sample, reduces the time taken to get meaningful results and allows more accurate detection of AMR genes than assembly.

## Overview

Amira leverages the full length of long read sequences to differentiate multi-copy genes by their local genomic context. This is done by first identifying the genes on each sequencing read and using the gene calls to construct a *de Bruijn* graph (DBG) in gene space. Following error correction, the reads containing different copies of multi-copy AMR genes can be clustered together based on their path in the graph, then assembled to obtain the nucleotide sequence.

## Prerequisites

Amira requires Python and three additional non-Python tools for optimal functionality:

- **Python >=3.10,<3.13**.
- **Poetry** to manage the Python dependencies.
- **Pandora >=0.12.0** to identify the genes on each sequencing read.
- **minimap2** for sequence alignment.
- **samtools** for processing alignments.
- **racon** for allele polishing.
- **Jellyfish** for cellular copy number estimation.

## Installation

Follow these steps to install Amira and its dependencies.

### With Singularity (preferred method)

Amira and its dependencies can be run through Singularity. A prebuilt singularity container can be found with each Release. Run this to build the container:

```bash
sudo singularity build amira.img Singularity.def
```

You can then run amira with:
```bash
singularity exec amira.img amira --help
```

### With PyPI

Amira can be installed from PyPI. **NOTE: You will need to install the non-Python dependencies separately if you opt for this method.**
```bash
pip install amira-amr
```
Amira can then be run with:
```bash
amira --help
```

### From source

#### Step 1: Clone the Amira Repository

Open a terminal and run the following command to clone the repository and navigate into it:
```bash
git clone https://github.com/Danderson123/Amira && cd Amira
```
#### Step 2: Install Poetry
Amira’s dependencies are managed with Poetry. Install Poetry by running:
```bash
pip install poetry
```
#### Step 3: Install Python Dependencies
Once Poetry is installed, use it to set up Amira’s python dependencies.
```bash
poetry install
```
You will need to install the non-Python dependencies separately if you opt for this method.

### Installing Non-Python Dependencies
Amira requires Pandora, minimap2, racon and Jellyfish. Follow the links below for instructions on building binaries for each tool. **The easiest way to run Pandora is through a binary that is precompiled for Linux, available from [here](https://github.com/iqbal-lab-org/pandora/releases/download/0.12.0-alpha.0/pandora-linux-precompiled-v0.12.0-alpha.0).**

- [Pandora installation guide](https://github.com/iqbal-lab-org/pandora?tab=readme-ov-file#installation)
- [minimap2 installation guide](https://github.com/lh3/minimap2)
- [samtools installation guide](https://www.htslib.org/download/)
- [racon installation guide](https://github.com/isovic/racon)
- [Jellyfish installation guide](https://github.com/gmarcais/Jellyfish)

After installation, make a note of the paths to these binaries as they will be required when running Amira.

## Pre-built species-specific panRGs
[Pandora](https://github.com/iqbal-lab-org/pandora) uses species-specific reference pan-genomes (panRGs) to identify the genes on each sequencing read (see above for instructions to install Pandora). Click the relevant link below to download a panRG to run Amira on your favorite bacterial species. If we do not currently support a species you are interested in then we are more than happy to build one, please let us know via a GitHub issue!
* [*Escherichia coli*](https://figshare.com/ndownloader/files/54318899)
* [*Klebsiella pneumoniae*](https://figshare.com/ndownloader/files/53398349)
* [*Enterococcus faecium*](https://figshare.com/ndownloader/files/53395052)
* [*Staphylococcus aureus*](https://figshare.com/ndownloader/files/53577833)
* [*Streptococcus pneumoniae*](https://figshare.com/ndownloader/files/53577887)

**Note**: These panRGs can currently detect all acquired AMR genes in the [NCBI Bacterial Antimicrobial Resistance Reference Gene Database](https://www.ncbi.nlm.nih.gov/bioproject/313047) as of 21st August 2024.

## Running Amira

Once you have installed the Python dependencies, Pandora, Racon and Minimap2, and downloaded the panRG for the species you interested in, Amira can be run with this command.

```bash
amira --reads <PATH TO READ FASTQ> --output <OUTPUT_DIR> --species <SPECIES> --panRG-path <PANRG PATH> --pandora-path <PATH TO PANDORA BINARY --racon-path <PATH TO RACON BINARY> --minimap2-path <PATH TO MINIMAP2 BINARY> --samtools-path <PATH TO SAMTOOLS BINARY> --cores <CPUS>
```

## For developers

Amira can also be run on the output of Pandora directly, or from JSON files listing the genes and gene positions on each sequencing read.

### Running with Pandora
After installing Pandora, you can call the genes on your sequencing reads using this command:
```bash
pandora map -t <THREADS> --min-gene-coverage-proportion 0.5 --max-covg 10000 -o pandora_map_output <PANRG PATH> <PATH TO READ FASTQ>
```
Amira can then be run directly on the output of Pandora using this command:
```bash
amira --pandoraSam pandora_map_output/*.sam --pandoraConsensus pandora_map_output/pandora.consensus.fq.gz --panRG-path <PANRG PATH> --reads <PATH TO READ FASTQ> --output amira_output --species <SPECIES> --minimum-length-proportion 0.5 --maximum-length-proportion 1.5 --cores <CPUS> --racon-path <PATH TO RACON BINARY> --minimap2-path <PATH TO MINIMAP2 BINARY> --samtools-path <PATH TO SAMTOOLS BINARY>
 ```

### Running from JSON
To run Amira from the JSON files, you can use this command:
```
amira --pandoraJSON <PATH TO GENE CALL JSON> --gene-positions <PATH TO GENE POSITION JSON> --pandoraConsensus <PATH TO PANDORA CONSENSUS FASTQ> --reads <PATH TO READ FASTQ> --output <OUTPUT DIRECTORY> --panRG-path <PANRG PATH> --species <SPECIES> --racon-path <PATH TO RACON BINARY> --minimap2-path <PATH TO MINIMAP2 BINARY> --cores <CPUS>
```

####  JSON example

Some example JSON data can be downloaded from [here](https://drive.google.com/drive/folders/1mQ8JmzVhFiNkgRy5_1iFQrqV2TLNnlQ4). Amira can then be run using this command:
```
amira --pandoraJSON test_data/gene_calls_with_gene_filtering.json --gene-positions test_data/gene_positions_with_gene_filtering.json --pandoraConsensus test_data/pandora.consensus.fq.gz --reads test_data/SRR23044220_1.fastq.gz --output amira_output --species Escherichia_coli --panRG-path . --racon-path <PATH TO RACON BINARY> --minimap2-path <PATH TO MINIMAP2 BINARY> --samtools-path <PATH TO SAMTOOLS BINARY> --debug --cores <CPUS>
```

### Additional options
For additional options and configurations, run:
```bash
amira --help
```
## Citation
```
@article {Anderson2025.05.16.654303,
	author = {Anderson, Daniel and Lima, Leandro and Le, Trieu and Judd, Louise M and Wick, Ryan R and Iqbal, Zamin},
	title = {Amira: gene-space de Bruijn graphs to improve the detection of AMR genes from bacterial long reads},
	elocation-id = {2025.05.16.654303},
	year = {2025},
	doi = {10.1101/2025.05.16.654303},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Accurate detection of antimicrobial resistance (AMR) genes is essential for the surveillance, epidemiology and genotypic prediction of AMR. This is typically done by generating an assembly from the sequencing reads of a bacterial isolate and running AMR gene detection tools on the assembly. However, despite advances in long-read sequencing that have greatly improved the quality and completeness of bacterial genome assemblies, assembly tools remain prone to large-scale errors caused by repeats in the genome, leading to inaccurate detection of AMR gene content and consequent impact on resistance prediction. In this work we present Amira, a tool to detect AMR genes directly from unassembled long-read sequencing data. Amira leverages the fact that multiple consecutive genes lie within a single read to construct gene-space de Bruijn graphs where the k-mer alphabet is the set of genes in the pan-genome of the species under study. Through this approach, the reads corresponding to different copies of AMR genes can be effectively separated based on the genomic context of the AMR genes, and used to infer the nucleotide sequence of each copy. Amira achieves significant improvements in genomic copy number recall and nucleotide accuracy, demonstrated through objective simulations and comparison with alternative read and assembly-based methods on samples with manually curated truth assemblies. Applied to a dataset of 32 Escherichia coli samples with diverse AMR gene content, Amira achieves a mean genomic-copy-number recall of 98.4\% with precision 97.9\% and nucleotide accuracy 99.9\%. Finally, we compare the number of fully (\&gt;90\%) present genes with good read support by Amira and Flye with AMRFinderPlus across all E. coli, K. pneumoniae and E. faecium nanopore datasets from the ENA (n=8593, 2469 and 420 respectively), quantifying the improvement in recall when handling heterogeneous data.Competing Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2025/05/18/2025.05.16.654303},
	eprint = {https://www.biorxiv.org/content/early/2025/05/18/2025.05.16.654303.full.pdf},
	journal = {bioRxiv}
}
```

## Contributing
If you’d like to contribute to Amira, please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bugfix (git checkout -b feature-name).
3. Commit your changes (git commit -m "Description of feature").
4. Push to the branch (git push origin feature-name).
5. Submit a pull request.

## License
This project is licensed under the Apache License 2.0 License. See the LICENSE file for details.

## Contact
For questions, feedback, or issues, please open an issue on GitHub or contact [Daniel Anderson](<mailto:dander@ebi.ac.uk>).

## Additional Resources
* [Pandora Documentation](https://github.com/iqbal-lab-org/pandora/wiki/Usage)