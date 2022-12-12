from construct_read import Read
from construct_read import GeneMer

class Node:
    def __init__(self,
                geneMer: GeneMer):
        self.geneMer = geneMer
        self.geneMerDirection = geneMer.get_geneMerDirection()
        self.nodeCoverage = 0
        self.setOfReads = set()
        self.forwardEdgeHashes = set()
        self.forwardEdges = []
        self.backwardEdgeHashes = set()
        self.backwardEdges = []
    def get_geneMer(self):
        """ return the GeneMer object represented by this node """
        return self.geneMer
    def assign_nodeId(self,
                    nodeId):
        self.nodeId = nodeId
        return self.nodeId
    def get_nodeId(self):
        """ return a unique integer node ID for this node """
        return self.nodeId
    def get_coverage(self) -> int:
        """ return the number of time we have seen this geneMer """
        return self.nodeCoverage
    def increment_coverage(self) -> int:
        """ increase the coverage of this node by 1 and return the new coverage """
        self.nodeCoverage += 1
        return self.nodeCoverage
    def get_reads(self) -> list:
        """ return a generator for the list of Read objects associated with this node """
        for read in self.setOfReads:
            yield read
    def add_read(self,
                read: Read):
        """ add a read object to a list of reads for this node """
        self.setOfReads.add(read)
    def remove_read(self,
                    read: Read):
        """ remove a read from the list of reads for this node """
        if read in self.setOfReads:
            readList = list(self.setOfReads)
            mask = readList.index(read)
            del readList[mask]
            self.setOfReads = readList
    def add_forward_edge(self,
                        forwardEdge):
        """ add a forward edge if it is not already in the forward edge list """
        edgeHash = forwardEdge.__hash__()
        if not edgeHash in self.forwardEdgeHashes:
            self.forwardEdges.append(forwardEdge)
            self.forwardEdgeHashes.add(edgeHash)
        else:
            mask = list(self.forwardEdgeHashes).index(edgeHash)
            self.forwardEdges[mask].increment_coverage()
    def add_edge(self,
                edge):
        if self.geneMerDirection == -1:
            self.add_forward_edge(edge)
        else:
            self.add_backward_edge(edge)
    def get_forward_edges(self):
        """ return a list of integers of node identifiers connected to this node by a forward edge """
        return self.forwardEdges
    def add_backward_edge(self,
                        backwardEdge):
        """ add a forward edge if it is not already in the forward edge list """
        edgeHash = backwardEdge.__hash__()
        if not edgeHash in self.backwardEdgeHashes:
            self.backwardEdges.append(backwardEdge)
            self.backwardEdgeHashes.add(edgeHash)
        else:
            mask = list(self.backwardEdgeHashes).index(edgeHash)
            self.backwardEdgeHashes[mask].increment_coverage()
    def get_backward_edges(self):
        """ return a list of integers of node identifiers connected to this node by a backward edge """
        return self.backwardEdges
    def __eq__(self,
            otherNode):
        """ check if two nodes are identical """
        return self == otherNode
    def __hash__(self):
        """ return a hash of the canonical gene mer to check if two nodes represent the same gene-mer """
        return self.geneMer.__hash__()
