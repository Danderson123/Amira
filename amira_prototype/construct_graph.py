import numpy as np

from construct_node import Node
from construct_edge import Edge
from construct_gene_mer import GeneMer

class GeneMerGraph:

    def __init__(self,
                readDict,
                kmerSize,
                minGeneMerCoverage,
                minEdgeCoverage):
        self._reads = readDict
        self._kmerSize = kmerSize
        self._minGeneMerCoverage = minGeneMerCoverage
        self._minEdgeCoverage = minEdgeCoverage
        self._nodes = {}
        self._edges = {}
    def get_reads(self):
        """ return a dictionary of all reads and their genes """
        return self._reads
    def get_kmerSize(self):
        """ return an integer of the gene-mer size """
        return self._kmerSize
    def get_minGeneMerCoverage(self):
        """ return an integer of the minimum coverage for nodes """
        return self._minGeneMerCoverage
    def get_minEdgeCoverage(self):
        """ return an integer of the minimum coverage for edges """
        return self._minEdgeCoverage
    def get_nodes(self):
        """ return the node dictionary """
        return self._nodes
    def get_edges(self):
        """ return the edge dictionary """
        return self._edges
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
        Add a gene mer to the graph if it does not exist, get the node.
        Returns the Node itself
        """
        nodeHash = geneMer.__hash__()
        if not nodeHash in self.get_nodes():
            # convert the gene-mer to a node
            node = Node(geneMer)
            # set the coverage to 1
            node.increment_node_coverage()
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
    def check_no_strand_in_genes(self,
                                listOfGenes: list) -> bool:
        """ returns a bool of whether any genes in a list have a stand prefix """
        return all(g[0] != "+" and g[0] != "-" for g in listOfGenes)
    def get_nodes_containing(self,
                            listOfGenes: list) -> list:
        """
        Return all nodes that contain a given gene.
        This is useful for later, when we want to get nodes that contain AMR genes.
        """
        # ensure there is no strand information present on the requested genes
        assert self.check_no_strand_in_genes(listOfGenes), "Strand information cannot be present for any specified genes"
        selectedNodes = []
        # iterate through the nodes in the graph
        for node in self.all_nodes():
            # get the canonical geneMer for this node
            geneMer = node.get_canonical_geneMer()
            # get the genes in this gene mer
            geneNames = [g.get_name() for g in geneMer]
            # add node to the selectedNodes list if it contains any of the specified genes
            if any(inputGene in geneNames for inputGene in listOfGenes):
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
        """
        Add the edge to the graph edges if it is not present, else increment the edge coverage by 1.
        Returns the edge.
        """
        # get the hash for the edge we want
        edgeHash = edge.__hash__()
        # see if the edge is in the edge dictionary
        if not edgeHash in self.get_edges():
            self._edges[edgeHash] = edge
        self._edges[edgeHash].increment_edge_coverage()
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
        return sourceNode, targetNode
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
        self.remove_edge_from_edges(edgeHash)
    def generate_gml(self,
                    outputDirectory):
        """ write a gml of the filtered graph to the output directory """
