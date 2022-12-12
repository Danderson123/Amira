import numpy as np

from construct_node import Node
from construct_edge import Edge
from construct_gene_mer import GeneMer
from construct_gene import Gene

class GeneMerGraph:
    def __init__(self,
                reads,
                kmerSize,
                minGeneMerCoverage,
                minEdgeCoverage):
        self.reads = reads
        self.kmerSize = kmerSize
        self.minGeneMerCoverage = minGeneMerCoverage
        self.minEdgeCoverage = minEdgeCoverage
        self.graph = []
        self.nodeHashes = []
        self.currentNodeId = 0
    def increment_nodeID(self):
        self.currentNodeId += 1
        return self.currentNodeId
    def determine_node_index(self,
                            geneMer):
        geneMerHash = geneMer.__hash__()
        assert geneMerHash in self.nodeHashes, "This gene-mer is not in the graph"
        mask = self.nodeHashes.index(geneMerHash)
        return mask
    def add_node(self,
                geneMer: GeneMer,
                readId: str) -> Node:
        """ add a gene mer to the graph if it does not exist, else increase the node coverage by 1. Returns the Node itself """
        geneMerHash = geneMer.__hash__()
        if not geneMerHash in self.nodeHashes:
            node = Node(geneMer)
            node.assign_nodeId(self.currentNodeId)
            self.graph.append(node)
            self.nodeHashes.append(geneMerHash)
            self.increment_nodeID()
        mask = self.determine_node_index(geneMer)
        self.graph[mask].increment_coverage()
        self.graph[mask].add_read(readId)
        return self.graph[mask]
    def get_node(self,
                geneMer: GeneMer) -> Node:
        """ get a node given a GeneMer """
        mask = self.determine_node_index(geneMer)
        return self.graph[mask]
    def get_nodes_containing_a_gene(self,
                                    gene: Gene) -> list:
        """ return all nodes that contain a given gene. Useful for later, when we want to get nodes that contain AMR genes, and build unitigs from them """
    def remove_node(self, node: Node):
        """ remove a node from the graph"""
        geneMerHash = Node.__hash__()
        assert geneMerHash in self.nodeHashes, "This node is not in the graph"
        mask = self.nodeHashes.index(geneMerHash)
        del self.graph[mask]
        del self.nodeHashes[mask]
    def get_degree(self, node: Node) -> int:
        """ return an integer of the number of neighbours for this node """
        numberOfForwardEdges = len(Node.get_forward_edges())
        numberOfBackwardEdges = len(Node.get_backward_edges())
        return numberOfForwardEdges + numberOfBackwardEdges
    def get_forward_edges(self, node: Node) -> list:
        """ return a list of integers of node identifiers connected to this node by a forward edge """
        return Node.get_forward_edges()
    def get_backward_edges(self, node: Node) -> list:
        """ return a list of integers of node identifiers connected to this node by a backward edge """
        return Node.get_backward_edges()
    def add_edge(self, edge: Edge):
        """ add a forward and backward edge if they do not exist, else increment the edge coverage by 1 """
        sourceNode = edge.get_sourceNode()
        sourceNodeHash = sourceNode.__hash__()
        targetNode = edge.get_targetNode()
        targetNodeHash = targetNode.__hash__()
    def remove_edge(self, edge: Edge):
        """ remove an edge """
    def get_total_number_of_nodes(self) -> int:
        """ return an integer of the total number of nodes in the filtered graph """
        return len(self.graph)
    def get_total_number_of_edges(self) -> int:
        """ return an integer of the total number of edges in the filtered graph """
    def get_total_number_of_reads(self) -> int:
        """ return an integer of the total number of reads contributing to the filtered graph """
        return len(self.reads)
    def all_nodes(self):
        """ return an iterator for all nodes in the graph and their attributes """
        for node in self.graph:
            yield node
    def generate_gml(self,
                    outputDirectory):
        """ write a gml of the filtered graph to the output directory """
