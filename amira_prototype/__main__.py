import argparse
from cigar import Cigar
import gzip
import matplotlib.pyplot as plt
import os
import pysam
import sys

from construct_graph import GeneMerGraph
from construct_unitig import Unitigs, UnitigBuilder, parse_fastq_lines
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
    # initialise the UnitigBuilder class
    unitigTools = UnitigBuilder(unitigsOfInterest)
    # generate a visualisation of the unitigs
    readFiles = unitigTools.separate_reads(args.output_dir,
                                        args.readfile)
    # run racon on the subsetted reads
    for r in readFiles:
        unitigTools.run_raven(r,
                            args.raven_path)
    # make plots to visualise unitigs
    unitigTools.visualise_unitigs(readDict,
                                args.output_dir)
    sys.exit(0)

if __name__ == "__main__":
    main()