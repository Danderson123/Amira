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
        assert self.check_no_strand_in_genes(listOfGenes), "Strand information cannot be present for any specified genes"
        selectedNodes = []
        for node in self.all_nodes():
            geneMer = node.get_canonical_geneMer()
            geneNames = [g.get_name() for g in geneMer]
            if any(inputGene in geneNames for inputGene in listOfGenes):
                selectedNodes.append(node)
        return selectedNodes
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
        # remove the node from the node dictionary
        del self.get_nodes()[nodeHash]
        # remove the forward and backward edges from the edge dictionary
        for edgeHash in forward_edges_to_remove + backward_edges_to_remove:
            # confirm the edge is in the graph
            assert edgeHash in self.get_edges(), "Edge is not in the graph"
            del self.get_edges()[edgeHash]
    def generate_gml(self,
                    outputDirectory):
        """ write a gml of the filtered graph to the output directory """
