import argparse

def get_options():
    """define args from the command line"""
    parser = argparse.ArgumentParser(description='Build a prototype gene de Bruijn graph.')
    parser.add_argument('--input', dest='input_path',
                        help='gene-annotated read json', required=True)
    parser.add_argument('--output', dest='output_file', type=str, default="gene_de_Bruijn_graph",
                        help='file name for generated de Bruijn graph')
    parser.add_argument('-k', dest='kmer_length', type=int, default=3,
                        help='kmer length for the gene de Bruijn graph')
    parser.add_argument('-n', dest='node_min_coverage', type=int, default=1,
                        help='minimum threshold for k-mer coverage')
    parser.add_argument('-e', dest='edge_min_coverage', type=int, default=1,
                        help='minimum threshold for edge coverage between gene-mers')
    parser.add_argument('--rc', dest='rc', action='store_true')
    parser.add_argument('--refine', dest='refine', type=str, default=None,
                        help='path of file with a newline separated list of genes we are interested in')
    args = parser.parse_args()
    return args
