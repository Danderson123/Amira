from construct_read import Read
from construct_read import GeneMer

class Node:
    def __init__(self,
                geneMer: GeneMer):
        self.geneMer = geneMer
        self.canonicalGeneMer = geneMer.get_canonical_geneMer()
        self.geneMerHash = geneMer.__hash__()
        self.nodeCoverage = 0
        self.listOfReads = []
        self.forwardEdgeHashes = []
        self.backwardEdgeHashes = []
        self._color = None
    def get_geneMer(self):
        """ return the GeneMer object represented by this node """
        return self.geneMer
    def get_canonical_geneMer(self):
        """ return the list of canonical gene objects represented by this node """
        return self.canonicalGeneMer
    def assign_nodeId(self,
                    nodeId):
        self.nodeId = nodeId
        return self.get_nodeId()
    def get_nodeId(self):
        """ return a unique integer node ID for this node """
        return self.nodeId
    def get_node_coverage(self) -> int:
        """ return the number of time we have seen this geneMer """
        return self.nodeCoverage
    def get_color(self):
        """ return the color of this node for debugging purposes """
        return self._color
    def increment_node_coverage(self) -> int:
        """ increase the coverage of this node by 1 and return the new coverage """
        self.nodeCoverage += 1
        return self.get_node_coverage()
    def get_list_of_reads(self) -> list:
        return self.listOfReads
    def get_reads(self) -> list:
        """ return a generator for the list of Read objects associated with this node """
        for read in self.get_list_of_reads():
            yield read
    def add_read(self,
                read: Read):
        """ add a read object to a list of reads for this node """
        if not read in self.get_list_of_reads():
            self.listOfReads.append(read)
    def remove_read(self,
                    read: Read):
        """ remove a read from the list of reads for this node """
        assert read in self.get_list_of_reads(), "This node does not contain the read: " + read.get_readId()
        readList = self.get_list_of_reads()
        mask = readList.index(read)
        del readList[mask]
        self.listOfReads = readList
    def add_forward_edge_hash(self,
                        forwardEdgeHash):
        """ add a forward edge hash if it is not already in the forward edge list, return the edge """
        if not forwardEdgeHash in self.get_forward_edge_hashes():
            self.forwardEdgeHashes.append(forwardEdgeHash)
        return self
    def remove_forward_edge_hash(self,
                                edgeHash):
        """ remove and edge hash from the list of forward edge hashes """
        assert edgeHash in self.get_forward_edge_hashes(), "This edge hash is not in the list of forward edge hashes"
        mask = self.get_forward_edge_hashes().index(edgeHash)
        del self.forwardEdgeHashes[mask]
    def get_forward_edge_hashes(self):
        """ return a list of the hashes of the backward edges """
        return self.forwardEdgeHashes
    def add_backward_edge_hash(self,
                        backwardEdgeHash):
        """ add a backward edge if it is not already in the backward edge list, return the edge """
        if not backwardEdgeHash in self.get_backward_edge_hashes():
            self.backwardEdgeHashes.append(backwardEdgeHash)
        return self
    def remove_backward_edge_hash(self,
                                edgeHash):
        """ remove and edge hash from the list of backward edge hashes """
        assert edgeHash in self.get_backward_edge_hashes(), "This edge hash is not in the list of backward edge hashes"
        mask = self.get_backward_edge_hashes().index(edgeHash)
        del self.backwardEdgeHashes[mask]
    def get_backward_edge_hashes(self):
        """ return a list of the hashes of the backward edges """
        return self.backwardEdgeHashes
    def assign_node_Id(self,
                    nodeId):
        """ assign an integer identifier to this node and return it """
        self.nodeId = nodeId
        return nodeId
    def get_node_Id(self):
        """ return the integer node ID for this node """
        return self.nodeId
    def color_node(self,
                listOfAMRGenes):
        """
        adds a color attribute to this node for debugging purposes
        - 0: node contains no AMR gene
        - 1: node contains an AMR gene and is not at a junction
        - 2: node contains an AMR gene and is at a junction
        - 3: node does not contain an AMR gene but is connected to a node containing an AMR gene
        """
        # get the canonical geneMer for this node
        geneMer = self.get_canonical_geneMer()
        # get the genes in this gene mer
        geneNames = [g.get_name() for g in geneMer]
        # see if there are any AMR genes in this node
        if not any(g in listOfAMRGenes for g in geneNames):
            self._color = 0
        else:
            degree = len(self.get_forward_edge_hashes()) + len(self.get_backward_edge_hashes())
            if not degree > 2:
                self._color = 1
            else:
                self._color = 2
    def __eq__(self,
            otherNode):
        """ check if two nodes are identical """
        return self.__hash__() == otherNode.__hash__() and self.get_node_coverage() == otherNode.get_node_coverage()
    def __hash__(self):
        """ return a hash of the canonical gene mer to check if two nodes represent the same gene-mer """
        return self.geneMerHash
