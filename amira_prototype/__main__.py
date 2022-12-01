import sys

from def_options import get_options
from process_json_data import load_json_data, classify_read_json
from build_graph import build_graph
from visualise_graph import generate_gml, get_filename, convert_node_to_string

from tqdm import tqdm

def main():
    args = get_options()
    # load the json as a dictionary
    readDict = load_json_data(args.input_path)
    # turn the annotated read json into more flexible classes
    readData, readIdNameDict, geneData, geneIdNameDict, refined_genes = classify_read_json(readDict,
                                                                                        args.refine)
    # build the graph
    graphData, geneMerIdNameDict = build_graph(readData,
                                            args.kmer_length,
                                            refined_genes)
    # generate the gml visualisation
    generate_gml(graphData,
                geneIdNameDict,
                geneMerIdNameDict,
                args.node_min_coverage,
                args.edge_min_coverage,
                get_filename(args.output_file,
                            args.kmer_length,
                            args.node_min_coverage,
                            args.edge_min_coverage,
                            args.rc))
    return sys.exit(0)

if __name__=="__main__":
    main()
