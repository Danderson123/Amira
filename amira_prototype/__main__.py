import argparse
from cigar import Cigar
import gzip
import matplotlib.pyplot as plt
import os
import pysam
import subprocess
import sys

from construct_graph import GeneMerGraph
from construct_unitig import Unitigs
from construct_gene import convert_int_strand_to_string

def get_options():
    """define args from the command line"""
    parser = argparse.ArgumentParser(description='Build a prototype gene de Bruijn graph.')
    parser.add_argument('--pandora', dest='pandoraDir',
                        help='Pandora map output directory', required=True)
    parser.add_argument('--readfile', dest='readfile',
                        help='path of gzipped long read fastq', required=True)
    parser.add_argument('--output', dest='output_dir', type=str, default="gene_de_Bruijn_graph",
                        help='directory for Amira outputs')
    parser.add_argument('-k', dest='geneMer_size', type=int, default=3,
                        help='kmer length for the gene de Bruijn graph')
    parser.add_argument('-n', dest='node_min_coverage', type=int, default=1,
                        help='minimum threshold for gene-mer coverage')
    parser.add_argument('-e', dest='edge_min_coverage', type=int, default=1,
                        help='minimum threshold for edge coverage between gene-mers')
    parser.add_argument('--gene-path', dest='path_to_interesting_genes',
                        help='path to a newline delimited file of genes of interest', required=True)
    parser.add_argument('--raven-path', dest='raven_path',
                        help='path to raven binary', required=True)
    args = parser.parse_args()
    return args

def get_read_start(cigarString):
    """ return an int of the 0 based position where the read region starts mapping to the gene """
    cigarSplit = list(cigarString.items())
    # check if there are any hard clipped bases at the start of the mapping
    if "H" in cigarSplit[0]:
        regionStart = cigarSplit[0][0]
    else:
        regionStart = 0
    return regionStart

def get_read_end(cigarString,
                regionStart):
    """ return an int of the 0 based position where the read region stops mapping to the gene """
    regionLength = len(cigarString)
    regionEnd = regionStart + regionLength
    return regionEnd, regionLength

def convert_pandora_output(pandoraDir):
    # load the pseudo SAM
    pandora_sam_content = pysam.AlignmentFile(os.path.join(pandoraDir, os.path.basename(pandoraDir) + ".filtered.sam"), "r")
    annotatedReads = {}
    readDict = {}
    # iterate through the read regions
    for read in pandora_sam_content.fetch():
        # check if the read has mapped to any regions
        if read.is_mapped:
            # read the CIGAR string for the mapped read region
            cigarString = Cigar(read.cigarstring)
            # get the read name
            readName = read.query_name
            # get the start base that the region maps to on the read
            regionStart = get_read_start(cigarString)
            # get the end base that the region maps to on the read
            regionEnd, regionLength = get_read_end(cigarString,
                                                regionStart)
            # get the name of the gene we have mapped to
            gene_name = read.reference_name
            # append the strand of the match to the name of the gene
            if not read.is_forward:
                gene_name = "-" + gene_name
            else:
                gene_name = "+" + gene_name
            if not readName in annotatedReads:
                annotatedReads[readName] = []
            if not gene_name[1:] in readDict:
                readDict[gene_name[1:]] = []
            # store the per read gene names, gene starts and gene ends
            readDict[gene_name[1:]].append(regionLength)
            # store the per read gene names
            annotatedReads[readName].append(gene_name)
    return annotatedReads, readDict

def get_consensus_lengths(pandoraDir):
    # Initialize an empty dictionary to store the results
    geneLengths = {}
    # Open the fastq file
    with gzip.open(os.path.join(pandoraDir, "pandora.consensus.fq.gz"), 'rt') as fh:
        # Iterate over the lines in the file
        for identifier, sequence, quality in parse_fastq_lines(fh):
            # Add the identifier and sequence to the results dictionary
            geneLengths[identifier] = len(sequence)
    return geneLengths

def get_gene_length(gene,
                geneLengths):
    geneName = gene.get_name()
    geneStrand = convert_int_strand_to_string(gene.get_strand())
    return max(geneLengths[geneName])

def visualise_unitigs(unitigsOfInterest,
                    geneLengths,
                    output_dir):
    for paralog in unitigsOfInterest:
        fig, ax = plt.subplots()
        count = 1
        for unitig in unitigsOfInterest[paralog]["unitigs"]:
            height = []
            union = set().union(*unitig)
            for gene in union:
                geneLength = get_gene_length(gene,
                                            geneLengths)
                if gene.get_name() == paralog:
                    colour = "r"
                else:
                    colour = "g"
                if gene.get_strand() == 1:
                    ax.arrow(sum(height), count, sum(height) + geneLength, 0, width = 0.05, head_width=0.2, head_length=500, color = colour, length_includes_head=True)
                else:
                    ax.arrow(sum(height) + geneLength, count, -1 * geneLength, 0, width = 0.05, head_width=0.2, head_length=500, color = colour, length_includes_head=True)
                height.append(geneLength)
            count +=1
        plt.yticks([i for i in range(1, count)])
        ax.set_ylabel(paralog)
        ax.set_xlabel("Length (bp)")
        plotFilename = os.path.join(output_dir, paralog, paralog) + ".pdf"
        plt.savefig(plotFilename)

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
    with gzip.open(fastq_file, 'rt') as fh:
        # Iterate over the lines in the file
        for identifier, sequence, quality in parse_fastq_lines(fh):
            # Add the identifier and sequence to the results dictionary
            results[identifier] = {"sequence": sequence,
                                "quality": quality}
    # Return the dictionary of results
    return results

def write_fastq(fastq_file,
                data):
    # Open the fastq file
    with gzip.open(fastq_file, 'wt') as fh:
        # Iterate over the data
        for identifier, value in data.items():
            # Write the identifier line
            fh.write(f'@{identifier}\n')
            # Write the sequence line
            fh.write(f'{value["sequence"]}\n')
            # Write the placeholder quality lines
            fh.write('+\n')
            fh.write(f'{value["quality"]}\n')

def separate_reads(unitigsOfInterest,
                output_dir,
                readFile):
    fastqContent = parse_fastq(readFile)
    readFiles = []
    for geneOfInterest in unitigsOfInterest:
        if not os.path.exists(os.path.join(output_dir, geneOfInterest)):
            os.mkdir(os.path.join(output_dir, geneOfInterest))
        for i in range(len(unitigsOfInterest[geneOfInterest]["unitigs"])):
            unitigReads = unitigsOfInterest[geneOfInterest]["reads"][i]
            subsettedReadData = {}
            for u in unitigReads:
                for r in u:
                    subsettedReadData[r] =  fastqContent[r]
            # write the per unitig fastq data
            filename = os.path.join(output_dir, geneOfInterest, geneOfInterest + "_" + str(i) + ".fastq.gz")
            write_fastq(filename,
                        subsettedReadData)
            readFiles.append(filename)
    return readFiles

def run_raven(inputFastq,
            raven_path):
    outputConsensus = inputFastq.replace(".fastq.gz", ".fasta")
    raven_command = " ".join([raven_path,
                            "-t",
                            "8",
                            inputFastq,
                            ">",
                            outputConsensus])
    subprocess.run(raven_command, shell=True, check=True)

def main():
    # get command line options
    args = get_options()
    # make the output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # convert the Pandora SAM file to a dictionary
    annotatedReads, readDict = convert_pandora_output(args.pandoraDir)
    # build the graph
    graph = GeneMerGraph(annotatedReads,
                        args.geneMer_size)
    graph.filter_graph(args.node_min_coverage,
                    args.edge_min_coverage)
    graph.generate_gml(os.path.join(args.output_dir, "gene_mer_graph"))
    # import the list of genes of interest
    with open(args.path_to_interesting_genes, "r") as i:
        genesOfInterest = i.read().splitlines()
    # build the unitig object
    unitigs = Unitigs(graph,
                    genesOfInterest)
    # get the unitigs of the genes of interest
    unitigsOfInterest = unitigs.get_unitigs_of_interest()
    # generate a visualisation of the unitigs
    readFiles = separate_reads(unitigsOfInterest,
                            args.output_dir,
                            args.readfile)
    # run racon on the subsetted reads
    for r in readFiles:
        run_raven(r,
                args.raven_path)
    # make plots to visualise unitigs
    visualise_unitigs(unitigsOfInterest,
                    readDict,
                    args.output_dir)
    sys.exit(0)

if __name__ == "__main__":
    main()