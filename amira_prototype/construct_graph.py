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
        return self._reads
    def get_kmerSize(self):
        return self._kmerSize
    def get_minGeneMerCoverage(self):
        return self._minGeneMerCoverage
    def get_minEdgeCoverage(self):
        return self._minEdgeCoverage
    def get_nodes(self):
        return self._nodes
    def get_edges(self):
        return self._edges
    def add_node_to_nodes(self,
                        node,
                        nodeHash):
        self._nodes[nodeHash] = node
    def add_node(self,
                geneMer: GeneMer,
                readId: str) -> Node:
        """
        Add a gene mer to the graph if it does not exist, else increase the node coverage by 1.
        Returns the Node itself
        """
        nodeHash = geneMer.__hash__()
        if not nodeHash in self.get_nodes():
            # convert the gene-mer to a node
            node = Node(geneMer)
            # add the node to the graph
            self.add_node_to_nodes(node,
                                nodeHash)
        else:
            node = self.get_nodes()[nodeHash]
        # increment the node coverage
        node.increment_node_coverage()
        # add the read ID to the read
        node.add_read(readId)
        return self.get_nodes()[nodeHash]
    def get_node(self,
                geneMer: GeneMer) -> Node:
        """ get a node given a GeneMer """
        nodeHash = geneMer.__hash__()
        assert nodeHash in self.get_nodes(), "This gene-mer is not in the graph"
        return self.get_nodes()[nodeHash]
    def remove_node(self, node: Node):
        """ change the value of a Node to None by index """
    def get_nodes_containing(self,
                            listOfAMRGenes: list) -> list:
        """ return all nodes that contain a given gene. Useful for later, when we want to get nodes that contain AMR genes, and build unitigs from them """
    def add_edge(self,
                sourceGeneMer: GeneMer,
                targetGeneMer: GeneMer,
                read):
        """
        Add a forward and backward edge if they do not exist, else increment the edge coverage by 1.
        Also add the nodes to the graph if they are not present already
        """
    def get_degree(self, node: Node) -> int:
        """ return an integer of the number of neighbours for this node """
    def get_forward_edges(self,
                        node: Node) -> list:
        """ return a list of integers of node identifiers connected to this node by a forward edge """
        return node.get_forward_edges()
    def get_backward_edges(self,
                        node: Node) -> list:
        """ return a list of integers of node identifiers connected to this node by a backward edge """
        return node.get_backward_edges()
    def remove_edge(self, edge: Edge):
        """ remove an edge """
    def get_graph(self) -> list:
        """ return the list of nodes in this graph """
    def get_total_number_of_nodes(self) -> int:
        """ return an integer of the total number of nodes in the filtered graph """
        return len(self.get_graph())
    def get_total_number_of_edges(self) -> int:
        """ return an integer of the total number of edges in the filtered graph """
    def get_total_number_of_reads(self) -> int:
        """ return an integer of the total number of reads contributing to the filtered graph """
        return len(self.get_reads())
    def all_nodes(self):
        """ return an iterator for all nodes in the graph and their attributes """
    def generate_gml(self,
                    outputDirectory):
        """ write a gml of the filtered graph to the output directory """
