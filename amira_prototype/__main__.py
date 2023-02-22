import argparse
from cigar import Cigar
import os
import pysam
import sys
from tqdm import tqdm

from construct_graph import GeneMerGraph
from construct_unitig import UnitigTools

def get_options():
    """define args from the command line"""
    parser = argparse.ArgumentParser(description='Build a prototype gene de Bruijn graph.')
    parser.add_argument('--pandoraSam', dest='pandoraSam',
                        help='Pandora map SAM file path', required=True)
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
    parser.add_argument('--flye-path', dest='flye_path',
                        help='path to Flye binary', default=None, required=False)
    parser.add_argument('--raven-path', dest='raven_path',
                        help='path to Raven binary', default=None, required=False)
    parser.add_argument('--use-consensus', dest='use_pandora_consensus',
                        help='polish the pandora consensus of each gene to recover AMR alleles',
                        action='store_true' , default=None, required=False)
    parser.add_argument("--pandora-consensus", dest="pandoraConsensus",
                        help='path to the pandora consensus fastq for this sample', default=None, required=False)
    parser.add_argument('--racon-path', dest='racon_path',
                        help='path to Racon binary', default=None, required=False)
    parser.add_argument('--threads', dest='threads', type=int, default=1,
                        help='number of threads to use')
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

def determine_gene_strand(read):
    if not read.is_forward:
        gene_name = "-" + read.reference_name.replace("~~~", ";")
    else:
        gene_name = "+" + read.reference_name.replace("~~~", ";")
    return gene_name

def convert_pandora_output(pandoraSam):
    # load the pseudo SAM
    pandora_sam_content = pysam.AlignmentFile(pandoraSam, "rb")
    annotatedReads = {}
    readDict = {}
    # iterate through the read regions
    for read in pandora_sam_content.fetch():
        # convert the cigarsting to a Cigar object
        cigarString = Cigar(read.cigarstring)
        # check if the read has mapped to any regions
        if read.is_mapped:
            # get the start base that the region maps to on the read
            regionStart = get_read_start(cigarString)
            # get the end base that the region maps to on the read
            regionEnd, regionLength = get_read_end(cigarString,
                                                regionStart)
            # append the strand of the match to the name of the gene
            gene_name = determine_gene_strand(read)
            if not read.query_name in annotatedReads:
                annotatedReads[read.query_name] = []
            if not gene_name[1:] in readDict:
                readDict[gene_name[1:]] = []
            # store the per read gene names, gene starts and gene ends
            readDict[gene_name[1:]].append(regionLength)
            # store the per read gene names
            annotatedReads[read.query_name].append(gene_name)
    return annotatedReads, readDict

def main():
    # get command line options
    args = get_options()
    # make the output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # convert the Pandora SAM file to a dictionary
    sys.stderr.write("\nAmira: loading Pandora SAM\n")
    annotatedReads, readDict = convert_pandora_output(args.pandoraSam)
    # build the graph
    sys.stderr.write("\nAmira: building gene-mer graph\n")
    graph = GeneMerGraph(annotatedReads,
                        args.geneMer_size)
    sys.stderr.write("\nAmira: filtering gene-mer graph\n")
    graph.filter_graph(args.node_min_coverage,
                    args.edge_min_coverage)
    sys.stderr.write("\nAmira: writing gene-mer graph\n")
    graph.generate_gml(os.path.join(args.output_dir, "gene_mer_graph"),
                    args.geneMer_size,
                    args.node_min_coverage,
                    args.edge_min_coverage)
    # import the list of genes of interest
    with open(args.path_to_interesting_genes, "r") as i:
        genesOfInterest = i.read().splitlines()
    # initialise the UnitigBuilder class
    unitigTools = UnitigTools(graph,
                            genesOfInterest)
    # generate a visualisation of the unitigs
    sys.stderr.write("\nAmira: separating paralog reads\n")
    readFiles, unitig_mapping = unitigTools.separate_reads(args.output_dir,
                                                        args.readfile)
    # run flye on the subsetted reads
    if args.flye_path:
        sys.stderr.write("\nAmira: assembling paralog reads with Flye\n")
        unitigTools.multithread_flye(readFiles,
                                    args.flye_path,
                                    args.threads)
    # run raven on subsetted reads and pandora consensus
    if args.raven_path:
        sys.stderr.write("\nAmira: assembling paralog reads with Raven\n")
        unitigTools.multithread_raven(readFiles,
                                    args.raven_path,
                                    args.threads)
    # use each gene's pandora consensus as the basis for polishing
    if args.use_pandora_consensus:
        sys.stderr.write("\nAmira: polishing Pandora AMR gene consensus sequences\n")
        unitigTools.polish_pandora_consensus(readFiles,
                                            args.racon_path,
                                            args.pandoraConsensus,
                                            genesOfInterest,
                                            args.threads)
    # make plots to visualise unitigs
    sys.stderr.write("\nAmira: generating unitig plots\n")
    unitigTools.visualise_unitigs(readDict,
                                unitig_mapping,
                                args.output_dir)
    sys.exit(0)

if __name__ == "__main__":
    main()