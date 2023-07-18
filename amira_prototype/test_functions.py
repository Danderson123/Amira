import gzip
import os
from tqdm import tqdm

from construct_graph import GeneMerGraph

class TestUnitigTools:
    def __init__(self,
                graph: GeneMerGraph,
                listOfGenes: list,
                readFile: str,
                output_dir: str):
        self._graph = graph
        self._listOfGenes = listOfGenes
        self._fastqContent = parse_fastq(readFile)
        self._output_dir = output_dir
    def get_graph(self):
        """ returns the geneMerGraph """
        return self._graph
    def get_selected_genes(self):
        """ returns the list of selected genes """
        return self._listOfGenes
    def get_fastqContent(self):
        """ returns a dictionary of all read data in the input fastq """
        return self._fastqContent
    def get_output_dir(self):
        """ returns a string of the output directory """
        return self._output_dir
    def get_read_data(self,
                    readId):
        """ return a dictionary of the data for this read """
        return self.get_fastqContent()[readId]
    def get_nodes_of_interest(self,
                            geneOfInterest):
        """ extracts the graph nodes containing the genes of interest and returns them as a list """
        return self.get_graph().get_nodes_containing(geneOfInterest)
    def get_all_nodes_containing_AMR_genes(self):
        """ return a dictionary of nodes containing AMR genes """
        AMRNodes = {}
        # iterate through the list of specified genes
        for geneOfInterest in tqdm(self.get_selected_genes()):
            # get the graph nodes containing this gene
            nodesOfInterest = self.get_nodes_of_interest(geneOfInterest)
            # add nodes of interest to the AMR node set
            for n in nodesOfInterest:
                AMRNodes[n.__hash__()] = n
        return AMRNodes
    def get_junctions(self):
        nodeJunctions = set()
        # get nodes that are anchors for traversing in the forward direction of the graph
        for node in self.get_graph().all_nodes():
            # if the degree of this node is more then this is a junction
            if self.get_graph().get_degree(node) > 2:
                nodeJunctions.add(node.__hash__())
        return nodeJunctions
    def get_AMR_anchors_and_junctions(self,
                                    AMRNodes):
        # get the graph
        graph = self.get_graph()
        # initialise the node anchor and junction sets
        nodeAnchors = set()
        nodeJunctions = set()
        # get nodes that are anchors for traversing in the forward direction of the graph
        for nodeHash in AMRNodes:
            node = graph.get_node_by_hash(nodeHash)
            if not graph.get_degree(node) > 2:
                # get the forward neighbors for this node
                forward_neighbors = graph.get_forward_neighbors(node)
                # get the backward neighbors for this node
                backward_neighbors = graph.get_backward_neighbors(node)
                # get the forward neighbors that contain an AMR gene
                forwardAMRNodes = [n for n in forward_neighbors if n.__hash__() in AMRNodes]
                # get the backward neighbors that contain an AMR gene
                backwardAMRNodes = [n for n in backward_neighbors if n.__hash__() in AMRNodes]
                if len(backwardAMRNodes) == 0 or len(forwardAMRNodes) == 0:
                    nodeAnchors.add(nodeHash)
            # if the number of backward neighbors is not 1 or 0 then this is a forward anchor
            else:
                nodeJunctions.add(nodeHash)
        return nodeAnchors, nodeJunctions
    def split_into_consecutive_chunks(self, input_list):
        chunks = []
        chunk = [input_list[0]]
        for i in range(1, len(input_list)):
            if input_list[i] - input_list[i-1] == 1:
                chunk.append(input_list[i])
            else:
                chunks.append(chunk)
                chunk = [input_list[i]]
        chunks.append(chunk)  # add the last chunk
        return chunks
    def convert_reads_to_junctions_and_anchors(self,
                                            AMRNodes,
                                            nodeAnchors,
                                            nodeJunctions):
        graph = self.get_graph()
        chunkId = 0
        referenceReadsAsChunks = {}
        chunkToIdMapping = {}
        IdToChunkMapping = {}
        readsInChunks = {}
        for readId in tqdm(graph.get_readNodes()):
            readNodes = graph.get_readNodes()[readId]
            anchorIndices = [i for i in range(len(readNodes)) if readNodes[i] in nodeAnchors or readNodes[i] in nodeJunctions]
            if not len(anchorIndices) < 2:
                referenceReadsAsChunks[readId] = []
                chunks = [readNodes[anchorIndices[i]: anchorIndices[i+1] + 1] for i in range(len(anchorIndices) - 1)]
                chunks = [tuple(c) for c in chunks if all(nodeHash in AMRNodes for nodeHash in c)]
                for c in chunks:
                    if not c in chunkToIdMapping:
                        chunkToIdMapping[c] = chunkId
                        IdToChunkMapping[chunkToIdMapping[c]] = c
                        readsInChunks[chunkToIdMapping[c]] = set()
                        chunkId += 1
                    for nodeHash in list(c):
                        for readId in graph.get_node_by_hash(nodeHash).get_reads():
                            readsInChunks[chunkToIdMapping[c]].add(readId)
        return IdToChunkMapping, readsInChunks
    def separate_paralogs(self):
        # isolate nodes containing AMR genes
        AMRNodes = self.get_all_nodes_containing_AMR_genes()
        # get the AMR nodes that are anchors and nodes
        amrAnchors, amrJunctions = self.get_AMR_anchors_and_junctions(AMRNodes)
        # get the anchors and junctions on reads
        IdToChunkMapping, readsInChunks = self.convert_reads_to_junctions_and_anchors(AMRNodes,
                                                                                    amrAnchors,
                                                                                    amrJunctions)
        # make a gml of the chunks
        import networkx as nx
        # Create a new empty graph
        graph = nx.Graph()
        # Iterate over each block ID and its associated set of read IDs
        for block_id, read_ids in readsInChunks.items():
            attributes = {'reads': list(read_ids), "coverage": len(read_ids)}
            # Add a node for each block ID with the attributes
            graph.add_node(block_id, **attributes)
            # Iterate over other block IDs to find shared read IDs
            for other_block_id, other_read_ids in readsInChunks.items():
                # Skip the same block or blocks that have already been processed
                if other_block_id == block_id or other_block_id in graph.adj[block_id]:
                    continue
                # Check if there are any shared read IDs
                if len(read_ids.intersection(other_read_ids)) > 0:
                    # Add an edge between the blocks if they share at least one read ID
                    attributes = {'reads': list(read_ids.intersection(other_read_ids)), "coverage": len(read_ids.intersection(other_read_ids))}
                    graph.add_edge(block_id, other_block_id, **attributes)
        # Write the graph in GML format to the target path
        nx.write_gml(graph, os.path.join(self.get_output_dir(), "blocks_plot.gml"))
        if not os.path.exists(os.path.join(self.get_output_dir(), "read_clusters")):
            os.mkdir(os.path.join(self.get_output_dir(), "read_clusters"))
        seen_nodes = set()
        chains = []
        for node in tqdm(graph.nodes):
            if not node in seen_nodes:
                component = nx.node_connected_component(graph, node)
                chains.append(list(component))
                for connected_node in component:
                    seen_nodes.add(connected_node)
        for c in range(len(chains)):
            reads = set()
            for unitigId in chains[c]:
                reads.update(readsInChunks[unitigId])
            with open(os.path.join(self.get_output_dir(), "read_clusters", str(c + 1)) + ".txt", "w") as o:
                o.write("\n".join(list(reads)))

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