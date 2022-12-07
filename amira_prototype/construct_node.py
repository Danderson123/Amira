from construct_read import Read
from construct_read import GeneMer

class Node:
    def __init__(self,
                GeneMer: GeneMer,
                nodeId: int):
        self.GeneMer = GeneMer
        self.nodeCoverage = 0
        self.setOfReads = set()
        self.nodeId = nodeId
    def get_GeneMer(self):
        """ return the GeneMer object represented by this node """
        return self.GeneMer
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
    def __eq__(self,
            otherNode):
        """ check if two nodes are identical """
        return self == otherNode
    def __hash__(self):
        """ return a hash of the canonical gene mer to check if two nodes represent the same gene-mer """
        return self.GeneMer.__hash__()
