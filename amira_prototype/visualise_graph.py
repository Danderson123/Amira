import networkx as nx

def convert_node_to_string(nodeName,
                        geneIdNameDict):
    node_ints = [int(Id) for Id in nodeName.split(",")]
    node_strands = ["+" if Id > 0 else "-" for Id in node_ints]
    return "~~~".join([node_strands[i] + geneIdNameDict[abs(node_ints[i])] for i in range(len(node_ints))])

def get_filename(output_file,
                k,
                node_min_coverage,
                edge_min_coverage,
                rc):
    filename = [output_file]
    filename.append(str(k) + "_k")
    filename.append(str(node_min_coverage) + "_nc")
    filename.append(str(edge_min_coverage) + "_ec")
    if rc:
        filename.append("rc")
    filename.append("gml")
    return ".".join(filename)

def write_node_entry(node_id,
                    node_string,
                    node_coverage,
                    reads):
    node_entry = "\tnode\t[\n"
    node_entry += "\t\tid\t" + str(node_id) + "\n"
    node_entry += '\t\tlabel\t"' + node_string + '"\n'
    node_entry += "\t\tcoverage\t" + str(node_coverage) + "\n"
    node_entry += '\t\treads\t"' + reads + '"\n'
    node_entry += "\t]\n"
    return node_entry

def write_edge_entry(source_node,
                    target_node,
                    edge_coverage):
    edge_entry = "\tedge\t[\n"
    edge_entry += "\t\tsource\t" + str(source_node) + "\n"
    edge_entry += "\t\ttarget\t" + str(target_node) + "\n"
    edge_entry += "\t\tweight\t" + str(edge_coverage) + "\n"
    edge_entry += "\t]\n"
    return edge_entry

def generate_gml(graphData,
                geneIdNameDict,
                geneMerIdNameDict,
                node_min_coverage,
                edge_min_coverage,
                filename):
#    G=nx.Graph()
    G = "graph\t[\n"
    from tqdm import tqdm
    seen_nodes = set()
    seen_edges = set()
    for node in tqdm(graphData):
        source_gene_mer = convert_node_to_string(node.get_nodeName(),
                                                geneIdNameDict)
        if node.get_nodeCoverage() >= node_min_coverage:
            if not source_gene_mer in seen_nodes:
                this_source_node = write_node_entry(node.get_nodeId(),
                                                    source_gene_mer,
                                                    node.get_nodeCoverage(),
                                                    ",".join([str(i) for i in node.get_readIds()]))
                G += this_source_node
                seen_nodes.add(source_gene_mer)
            filtered_edges = [e for e in node.get_forwardEdgeList() + node.get_reverseEdgeList() if e.get_edgeWeight() > edge_min_coverage]
            for e in filtered_edges:
                absolute_targetId = e.get_targetNodeId()
                target_gene_mer = convert_node_to_string(geneMerIdNameDict[absolute_targetId],
                                                        geneIdNameDict)
                if not target_gene_mer in seen_nodes:
                    this_target_node = write_node_entry(absolute_targetId,
                                                        target_gene_mer,
                                                        graphData[absolute_targetId].get_nodeCoverage(),
                                                        ",".join([str(i) for i in graphData[absolute_targetId].get_readIds()]))
                    if not this_target_node in seen_nodes:
                        G += this_target_node
                        seen_nodes.add(target_gene_mer)
                sorted_edge = tuple(sorted([source_gene_mer, target_gene_mer]))
                if not sorted_edge in seen_edges:
                    edge_coverage = e.get_edgeWeight()
                    this_edge = write_edge_entry(node.get_nodeId(),
                                                absolute_targetId,
                                                edge_coverage)
                    G += this_edge
                    seen_edges.add(sorted_edge)
    G += "]"
    with open(filename, "w") as graphOut:
        graphOut.write(G)
    return