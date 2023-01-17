import os
from tqdm import tqdm

from construct_node import Node
from construct_edge import Edge
from construct_gene_mer import GeneMer
from construct_read import Read
from construct_gene import convert_int_strand_to_string

class GeneMerGraph:

    def __init__(self,
                readDict,
                kmerSize):
        self._reads = readDict
        self._kmerSize = kmerSize
        self._minNodeCoverage = 1
        self._minEdgeCoverage = 1
        self._nodes = {}
        self._edges = {}
        # initialise the graph
        for readId in tqdm(self.get_reads()):
            read = Read(readId,
                        self.get_reads()[readId])
            geneMers = [g for g in read.get_geneMers(self.get_kmerSize())]
            if not len(geneMers) == 0:
                for g in range(len(geneMers) - 1):
                    sourceNode = self.add_node(geneMers[g])
                    sourceNode.increment_node_coverage()
                    sourceNode.add_read(read.get_readId())
                    if not len(geneMers) == 1:
                        targetNode = self.add_node(geneMers[g+1])
                        targetNode.add_read(read.get_readId())
                        sourceToTargetEdge, reverseTargetToSourceEdge = self.add_edge(geneMers[g],
                                                                                    geneMers[g+1])
                        sourceToTargetEdge.increment_edge_coverage()
                        reverseTargetToSourceEdge.increment_edge_coverage()
                if not len(geneMers) == 1:
                    targetNodeHash = geneMers[-1].__hash__()
                    targetNode = self.get_node_by_hash(targetNodeHash)
                    targetNode.increment_node_coverage()
    def get_reads(self):
        """ return a dictionary of all reads and their genes """
        return self._reads
    def get_kmerSize(self):
        """ return an integer of the gene-mer size """
        return self._kmerSize
    def get_minEdgeCoverage(self):
        """ return an integer of the minimum coverage for edges """
        return self._minEdgeCoverage
    def get_nodes(self):
        """ return the node dictionary """
        return self._nodes
    def get_edges(self):
        """ return the edge dictionary """
        return self._edges
    def get_minNodeCoverage(self):
        """ return the minimum node coverage """
        return self._minNodeCoverage
    def get_minEdgeCoverage(self):
        """ return the minimum edge coverage """
        return self._minEdgeCoverage
    def all_nodes(self):
        """ return a generator for all nodes in the graph and their attributes """
        for nodeHash in self.get_nodes():
            yield self.get_nodes()[nodeHash]
    def add_node_to_nodes(self,
                        node,
                        nodeHash):
        """ add a node to the dictionary of nodes """
        self._nodes[nodeHash] = node
    def get_node_by_hash(self,
                        nodeHash):
        """ return a node object corresponding to a node hash """
        return self.get_nodes()[nodeHash]
    def add_node(self,
                geneMer: GeneMer) -> Node:
        """
        Add a gene mer to the graph if it does not exist, else get the node.
        Returns the Node itself
        """
        nodeHash = geneMer.__hash__()
        if not nodeHash in self.get_nodes():
            # convert the gene-mer to a node
            node = Node(geneMer)
            # add the node to the graph
            self.add_node_to_nodes(node,
                                nodeHash)
        return self.get_node_by_hash(nodeHash)
    def get_node(self,
                geneMer: GeneMer) -> Node:
        """ get a node given a GeneMer """
        # the gene-mer and node hashes are the same
        nodeHash = geneMer.__hash__()
        # confirm the node is in the node dictionary
        assert nodeHash in self.get_nodes(), "This gene-mer is not in the graph"
        # return the node
        return self.get_nodes()[nodeHash]
    def get_nodes_containing(self,
                            geneOfInterest: str) -> list:
        """
        Return all nodes that contain a given gene.
        This is useful for later, when we want to get nodes that contain AMR genes.
        """
        # ensure there is no strand information present on the requested genes
        assert not (geneOfInterest[0] == "+" or geneOfInterest[0] == "-"), "Strand information cannot be present for any specified genes"
        selectedNodes = []
        # iterate through the nodes in the graph
        for node in self.all_nodes():
            # get the canonical geneMer for this node
            geneMer = node.get_canonical_geneMer()
            # get the genes in this gene mer
            geneNames = [g.get_name() for g in geneMer]
            # add node to the selectedNodes list if it contains any of the specified genes
            if geneOfInterest in geneNames:
                # add the node to the list of selected nodes
                selectedNodes.append(node)
        return selectedNodes
    def create_edges(self,
                    sourceNode: Node,
                    targetNode: Node,
                    sourceGeneMerDirection: int,
                    targetGeneMerDirection: int):
        """ create and return the forward and reverse edges from source to target and target to source """
        # define the edge that will go from the source gene mer to the target gene mer
        sourceToTargetEdge = Edge(sourceNode,
                                targetNode,
                                sourceGeneMerDirection,
                                targetGeneMerDirection)
        # define the edge that will go from the reverse source gene mer to the reverse target gene mer
        reverseTargetToSourceEdge = Edge(targetNode,
                                        sourceNode,
                                        targetGeneMerDirection * -1,
                                        sourceGeneMerDirection * -1)
        return sourceToTargetEdge, reverseTargetToSourceEdge
    def get_edge_by_hash(self,
                        edgeHash):
        """ return an edge object corresponding to a edge hash """
        return self.get_edges()[edgeHash]
    def add_edge_to_edges(self,
                        edge):
        """Add the edge to the graph edges if it is not present, else increment the edge coverage by 1.
        Returns the edge.
        """
        # get the hash for the edge we want
        edgeHash = edge.__hash__()
        # see if the edge is in the edge dictionary
        if not edgeHash in self.get_edges():
            self._edges[edgeHash] = edge
        return self.get_edge_by_hash(edgeHash)
    def add_edges_to_graph(self,
                        sourceToTargetEdge,
                        reverseTargetToSourceEdge):
        """ add the edges to the dictionary of edges and return the updated edges """
        sourceToTargetEdge = self.add_edge_to_edges(sourceToTargetEdge)
        reverseTargetToSourceEdge = self.add_edge_to_edges(reverseTargetToSourceEdge)
        return sourceToTargetEdge, reverseTargetToSourceEdge
    def add_edge_to_node(self,
                        node,
                        edge):
        """
        Add the edge and reverse edge hashes to the node attributes.
        Return the forward and reverse edge hashes for this node.
        """
        if edge.get_sourceNodeDirection() == 1:
            edgeHash = edge.__hash__()
            node.add_forward_edge_hash(edgeHash)
        if edge.get_sourceNodeDirection() == -1:
            edgeHash = edge.__hash__()
            node.add_backward_edge_hash(edgeHash)
        return node
    def add_edge(self,
                sourceGeneMer: GeneMer,
                targetGeneMer: GeneMer):
        """
        Add a forward and backward edge if they do not exist, else increment the edge coverage by 1.
        Also add the nodes to the graph if they are not present already.
        Returns the updated source and target lists of forward and reverse edges.
        """
        # if the source is not in the graph then add it, else get the node
        sourceNode = self.add_node(sourceGeneMer)
        # if the target is not in the graph then add it, else get the node
        targetNode = self.add_node(targetGeneMer)
        # create the edge objects
        sourceToTargetEdge, reverseTargetToSourceEdge = self.create_edges(sourceNode,
                                                                        targetNode,
                                                                        sourceGeneMer.get_geneMerDirection(),
                                                                        targetGeneMer.get_geneMerDirection())
        # add the edges to the graph
        sourceToTargetEdge, reverseTargetToSourceEdge = self.add_edges_to_graph(sourceToTargetEdge,
                                                                                reverseTargetToSourceEdge)
        # add the hashes of the edges to the source and target nodes
        sourceNode = self.add_edge_to_node(sourceNode,
                                        sourceToTargetEdge)
        targetNode = self.add_edge_to_node(targetNode,
                                        reverseTargetToSourceEdge)
        return sourceToTargetEdge, reverseTargetToSourceEdge
    def get_degree(self,
                node: Node) -> int:
        """ return an integer of the number of neighbours for this node """
        degrees = len(node.get_forward_edge_hashes()) + len(node.get_backward_edge_hashes())
        return degrees
    def get_forward_edges(self,
                        node: Node) -> list:
        """ return a list of integers of node identifiers connected to this node by a forward edge """
        return node.get_forward_edge_hashes()
    def get_backward_edges(self,
                        node: Node) -> list:
        """ return a list of integers of node identifiers connected to this node by a backward edge """
        return node.get_backward_edge_hashes()
    def remove_edge_from_edges(self,
                            edgeHash):
        """ remove an edge object from the dictionary of all edges by hash """
        del self._edges[edgeHash]
        assert not edgeHash in self.get_edges(), "This edge was not removed from the graph successfully"
    def remove_edge(self,
                edgeHash: int):
        """ remove an edge from all graph edges and the node edges """
        # check that the edge exists in the graph
        assert edgeHash in self.get_edges(), "This edge is not in the graph"
        # get the edge object from all graph edges by hash
        edge = self.get_edges()[edgeHash]
        # get the source node object of the edge
        sourceNode = edge.get_sourceNode()
        # get the source node direction of the edge
        sourceNodeDirection = edge.get_sourceNodeDirection()
        # if the direction is forward then remove the edge from the forward edges for the node
        if sourceNodeDirection == 1:
            sourceNode.remove_forward_edge_hash(edgeHash)
        # if the direction is backward then remove the edge from the backward edges for the node
        if sourceNodeDirection == -1:
            sourceNode.remove_backward_edge_hash(edgeHash)
        # remove the edge object from all of the edges in the graph
        del self._edges[edgeHash]
    def get_total_number_of_nodes(self) -> int:
        """ return an integer of the total number of nodes in the filtered graph """
        return len(self.get_nodes())
    def get_total_number_of_edges(self) -> int:
        """ return an integer of the total number of edges in the filtered graph """
        return len(self.get_edges())
    def get_total_number_of_reads(self) -> int:
        """ return an integer of the total number of reads contributing to the filtered graph """
        return len(self.get_reads())
    def remove_node(self,
                node: Node):
        """ remove a node from the graph and all of its edges """
        # get the node hash
        nodeHash = node.__hash__()
        # confirm the node to remove is in the node dictionary
        assert nodeHash in self.get_nodes(), "This node is not in the graph"
        node_to_remove = self.get_nodes()[nodeHash]
        # get the forward edge hashes of the node to remove
        forward_edges_to_remove = node_to_remove.get_forward_edge_hashes()
        # get the backward edge hashes of the node to remove
        backward_edges_to_remove = node_to_remove.get_backward_edge_hashes()
        # remove the forward and backward edges from the edge dictionary
        for edgeHash in forward_edges_to_remove + backward_edges_to_remove:
            # confirm the edge is in the graph
            self.remove_edge(edgeHash)
        # remove the node from the node dictionary
        del self.get_nodes()[nodeHash]
    def set_minNodeCoverage(self,
                            minNodeCoverage):
        """ set the minimum node coverage for the graph and return the new node coverage """
        self._minNodeCoverage = minNodeCoverage
        return self.get_minNodeCoverage()
    def set_minEdgeCoverage(self,
                            minEdgeCoverage):
        """ set the minimum edge coverage for the graph and return the new edge coverage """
        self._minEdgeCoverage = minEdgeCoverage
        return self.get_minEdgeCoverage()
    def filter_graph(self,
                    minNodeCoverage,
                    minEdgeCoverage):
        """ filters the nodes and edges in the graph and then returns the filtered graph """
        # set the new node coverage
        minNodeCoverage = self.set_minNodeCoverage(minNodeCoverage)
        # set the new edge coverage
        minEdgeCoverage = self.set_minEdgeCoverage(minEdgeCoverage)
        # filter the nodes
        nodesToRemove = set()
        for nodeHash in self.get_nodes():
            # mark a node for removal if the coverage is less than the specified coverage
            if not self.get_nodes()[nodeHash].get_node_coverage() > minNodeCoverage - 1:
                nodesToRemove.add(self.get_nodes()[nodeHash])
        # filter the edges
        edgesToRemove = set()
        for edgeHash in self.get_edges():
            # mark an edge for removal if the coverage is less than the specified coverage
            if not self.get_edges()[edgeHash].get_edge_coverage() > minEdgeCoverage - 1:
                edgesToRemove.add(edgeHash)
            # mark an edge for removal if either of the nodes are marked for removal
            if any(n in nodesToRemove for n in [self.get_edges()[edgeHash].get_sourceNode(),
                                                self.get_edges()[edgeHash].get_targetNode()]):
                edgesToRemove.add(edgeHash)
        for e in edgesToRemove:
            self.remove_edge(e)
        for n in nodesToRemove:
            self.remove_node(n)
        return self
    def write_node_entry(self,
                        node_id,
                        node_string,
                        node_coverage,
                        reads):
        """ return a string of a gml node entry """
        if node_coverage > 100:
            node_coverage = 100
        node_entry = "\tnode\t[\n"
        node_entry += "\t\tid\t" + str(node_id) + "\n"
        node_entry += '\t\tlabel\t"' + node_string + '"\n'
        node_entry += "\t\tcoverage\t" + str(node_coverage) + "\n"
        node_entry += '\t\treads\t"' + ",".join(reads) + '"\n'
        node_entry += "\t]"
        return node_entry
    def write_edge_entry(self,
                        source_node,
                        target_node,
                        edge_coverage):
        """ return a string of a gml edge entry """
        edge_entry = "\tedge\t[\n"
        edge_entry += "\t\tsource\t" + str(source_node) + "\n"
        edge_entry += "\t\ttarget\t" + str(target_node) + "\n"
        edge_entry += "\t\tweight\t" + str(edge_coverage) + "\n"
        edge_entry += "\t]"
        return edge_entry
    def assign_Id_to_nodes(self):
        """ assign an integer ID to each node in the graph """
        currentNodeId = 0
        for node in self.all_nodes():
            assignedId = node.assign_node_Id(currentNodeId)
            assert assignedId == currentNodeId, "This node was assigned an incorrect ID"
            currentNodeId += 1
    def write_gml_to_file(self,
                        output_file,
                        gml_content):
        """ Writes the gml to a file titled <output_file>.gml. Returns nothing."""
        # check if there is a directory to make
        if os.path.dirname(output_file) != "":
            # make the directory if it does not exist
            if not os.path.exists(os.path.dirname(output_file)):
                os.mkdir(os.path.dirname(output_file))
        # write the gml to the output file
        with open(output_file + ".gml", "w") as outGml:
            outGml.write("\n".join(gml_content))
    def get_gene_mer_label(self,
                        sourceNode: Node) -> str:
        """ returns a string of the genes in the canonicall gene mer for a node and their strands """
        # get the canonical gene mer for the node
        geneMer = sourceNode.get_canonical_geneMer()
        # get a list of the gene strand + gene name for the canonical gene mer
        geneMerGenes = [convert_int_strand_to_string(g.get_strand()) + g.get_name() for g in geneMer]
        # return a string of the gene mer genes and strands
        return "~~~".join(geneMerGenes)
    def generate_gml(self,
                    output_file: str,
                    geneMerSize: int,
                    min_node_coverage: int,
                    min_edge_coverage: int):
        """ Write a gml of the filtered graph to the output directory. Returns the written content as a list """
        graph_data = ["graph\t["]
        self.assign_Id_to_nodes()
        # iterate through the nodes in the graph
        for sourceNode in self.all_nodes():
            # check that the coverage of the node is greater than the minimum specified node coverage
            # add the node entry
            nodeEntry = self.write_node_entry(sourceNode.get_node_Id(),
                                            self.get_gene_mer_label(sourceNode),
                                            sourceNode.get_node_coverage(),
                                            [read for read in sourceNode.get_reads()])
            graph_data.append(nodeEntry)
            # get the forward edges for this node
            nodeForwardEdgeHashes = sourceNode.get_forward_edge_hashes()
            nodeBackwardEdgeHashes = sourceNode.get_backward_edge_hashes()
            # get the edge objects corresponding to the forward edge hashes
            nodeEdges = [self.get_edge_by_hash(edgeHash) for edgeHash in nodeForwardEdgeHashes + nodeBackwardEdgeHashes]
            for edge in nodeEdges:
                targetNode = edge.get_targetNode()
                edgeEntry = self.write_edge_entry(sourceNode.get_node_Id(),
                                                targetNode.get_node_Id(),
                                                edge.get_edge_coverage())
                graph_data.append(edgeEntry)
        graph_data.append("]")
        output_file = ".".join([output_file,
                                str(geneMerSize),
                                str(min_node_coverage),
                                str(min_edge_coverage)])
        self.write_gml_to_file(output_file,
                            graph_data)
        return graph_data
