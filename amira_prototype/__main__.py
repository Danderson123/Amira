import argparse
from cigar import Cigar
import os
import pysam
import sys

from construct_graph import GeneMerGraph
from construct_unitig import Unitigs

def get_options():
    """define args from the command line"""
    parser = argparse.ArgumentParser(description='Build a prototype gene de Bruijn graph.')
    parser.add_argument('--input', dest='input_path',
                        help='path of Pandora SAM file', required=True)
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
    return regionEnd

def convert_pandora_output(pandoraSamPath):
    # load the pseudo SAM
    pandora_sam_content = pysam.AlignmentFile(pandoraSamPath, "r")
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
            regionEnd = get_read_end(cigarString,
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
                readDict[readName] = []
            # store the per read gene names, gene starts and gene ends
            readDict[readName].append({"gene name": gene_name,
                                        "region start": regionStart,
                                        "region end": regionEnd})
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
    annotatedReads, readDict = convert_pandora_output(args.input_path)
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
    sys.exit(0)

if __name__ == "__main__":
    main()