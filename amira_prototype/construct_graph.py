import numpy as np

from construct_node import Node
from construct_edge import Edge
from construct_gene_mer import GeneMer
from construct_gene import Gene
from construct_read import Read
class GeneMerGraph:

    def __init__(self,
                readDict,
                kmerSize,
                minGeneMerCoverage,
                minEdgeCoverage):

        self.readDict = readDict
        self.kmerSize = kmerSize
        self.minGeneMerCoverage = minGeneMerCoverage
        self.minEdgeCoverage = minEdgeCoverage
        self.graph = []
        self.nodeHashes = []
        self.currentNodeId = 0

        for read in readDict:
            this_read = Read(read,
                            readDict[read])
            for geneMer in this_read.get_geneMers(kmerSize):
                self.add_node(geneMer,
                            read)

    def increment_nodeID(self) -> int:
        """ increases the current node ID by 1 and returns an integer of the new value """
        self.currentNodeId += 1
        return self.currentNodeId
    def determine_node_index(self,
                            geneMer) -> int:
        """ returns an integer of the index for this node in the list of nodes in the graph """
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
    def remove_node(self, node: Node):
        """ change the value of a Node to None by index """
        geneMerHash = Node.__hash__()
        assert geneMerHash in self.nodeHashes, "This node is not in the graph"
        mask = self.nodeHashes.index(geneMerHash)
        self.graph[mask] = None
        self.nodeHashes[mask] = None
    def get_nodes_containing(self,
                            listOfAMRGenes: list) -> list:
        """ return all nodes that contain a given gene. Useful for later, when we want to get nodes that contain AMR genes, and build unitigs from them """
        assert all(g[0] != "+" and g[0] != "-" for g in listOfAMRGenes), "Strand information cannot be present for specified genes"
        for node in self.graph:
            geneMer = node.get_geneMer().get_canonical_geneMer()
            geneNames = [g.get_name() for g in geneMer]
            if any(inputGene in geneNames for inputGene in listOfAMRGenes):
                yield node
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
