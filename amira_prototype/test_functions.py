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
                chunks = [tuple(c) for c in chunks] #if all(nodeHash in AMRNodes for nodeHash in c)]
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
    def assign_reads_to_amr_genes(self):
        # isolate nodes containing AMR genes
        AMRNodes = self.get_all_nodes_containing_AMR_genes()
        # get the unitig for each AMR node
        unitigs = set()
        readsPerNode = {}
        for outerNodeHash in tqdm(AMRNodes):
            # get the unitig that this node is on
            u = [n for n in self.get_graph().get_linear_path_for_node(AMRNodes[outerNodeHash], True)]
            for innerNodeHash in u:
                readsPerNode[innerNodeHash] = [r for r in self.get_graph().get_node_by_hash(innerNodeHash).get_reads()]
            if not (tuple(u) in unitigs or tuple(list(reversed(u))) in unitigs):
                unitigs.add(tuple(u))
        unitigId = 0
        readsPerAMRGene = {}
        kmer_size = self.get_graph().get_kmerSize()
        for u in unitigs:
            coverages_for_this_unitig = []
            for nodeHash in list(u):
                node = self.get_graph().get_node_by_hash(nodeHash)
                coverages_for_this_unitig.append(node.get_node_coverage())
            if not len(u) == 1:
                # get the list of genes for this unitig
                genes = self.get_graph().follow_path_to_get_annotations(u)
            else:
                genes = self.get_graph().get_gene_mer_label(self.get_graph().get_node_by_hash(u[0])).split("~~~")
            # match up the nodes with the genes they contain
            nodes = list(u)
            node_indices_for_each_gene_index = {}
            for i in range(len(nodes)):
                gene_indices = [j for j in range(i, i+kmer_size)]
                for g in gene_indices:
                    if genes[g][1:] in self.get_selected_genes():
                        if not g in node_indices_for_each_gene_index:
                            node_indices_for_each_gene_index[g] = []
                        node_indices_for_each_gene_index[g].append(i)
            # get the reads for each gene index
            for i in node_indices_for_each_gene_index:
                gene = genes[i]
                nodes_containing_this_gene = [nodes[j] for j in node_indices_for_each_gene_index[i]]
                reads_per_gene = set()
                for nodeHash in nodes_containing_this_gene:
                    node = self.get_graph().get_node_by_hash(nodeHash)
                    coverages_for_this_unitig.append(node.get_node_coverage())
                    reads_for_this_node = [r for r in node.get_reads()]
                    reads_per_gene.update(reads_for_this_node)
                readsPerAMRGene[f"{gene[1:]}_{str(unitigId)}_{str(i)}"] = list(reads_per_gene)
            unitigId += 1
        import json
        with open(os.path.join(self.get_output_dir(), "reads_per_amr_gene.json"), "w") as o:
            o.write(json.dumps(readsPerAMRGene))

    def make_blocks_plots(self):
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
    if ".gz" in fastq_file:
        with gzip.open(fastq_file, 'rt') as fh:
            # Iterate over the lines in the file
            for identifier, sequence, quality in parse_fastq_lines(fh):
                # Add the identifier and sequence to the results dictionary
                results[identifier] = {"sequence": sequence,
                                    "quality": quality}
    else:
        with open(fastq_file, 'r') as fh:
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


class Unitig:
    def __init__(self,
                listOfNodes: list,
                listOfGenes: list,
                component_ID: int,
                unitig_ID: int):
        self._nodes = listOfNodes
        self._genes = listOfGenes

        def get_reverse_genes(listOfGenes):
            return ["+" + g[1:] if g[0] == "-" else "-" + g[1:] for g in list(reversed(listOfGenes))]

        self._reverse_genes = get_reverse_genes(self.get_genes())
        self._component_ID = component_ID
        self._unitig_ID = unitig_ID

        def get_read_union(listOfNodes):
            read_union = set()
            for node in listOfNodes:
                for read in node.get_reads():
                    read_union.add(read)
            return list(read_union)

        self._reads = get_read_union(self.get_nodes())
        self._coverage = len(self.get_reads())

    def get_nodes(self):
        """ returns an ordered list of nodes in this unitig """
        return self._nodes
    def get_genes(self):
        """ returns an ordered list of genes in this unitig """
        return self._genes
    def get_reverse_genes(self):
        """ returns an ordered list of the reverse complement genes for this unitig """
        return self._reverse_genes
    def get_component_ID(self):
        """ returns an integer of the component containing this unitig """
        return self._component_ID
    def get_unitig_ID(self):
        """ returns the integer identifier of this unitig """
        return self._unitig_ID
    def get_reads(self):
        """ returns a list of the reads in the nodes of this unitig """
        return self._reads
    def get_coverage(self):
        """ returns an integer of the number of reads supporting this unitig """
        return self._coverage