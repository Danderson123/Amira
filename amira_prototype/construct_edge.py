from construct_node import Node

def extract_node_hashes(firstNode,
                        secondNode):
    """ get the gene mer hashes for the first node and the second node """
    firstNodeHash = firstNode.__hash__()
    secondNodeHash = secondNode.__hash__()
    return firstNodeHash, secondNodeHash

def sort_node_hashes(firstNodeHash,
                    secondNodeHash):
    """ sorted the hashes in ascending order and assign the sourceNode as index 0 and targetNode as index 1 """
    sortedHashes = sorted([firstNodeHash,
                        secondNodeHash])
    # the sourceNode is at index 0
    sourceNode = sortedHashes[0]
    # the targetNode is at index 1
    targetNode = sortedHashes[1]
    return sourceNode, targetNode

def define_source_and_target(firstNode,
                            secondNode):
    """ decides which node will be source and which will be the target based on the smallest of the two hashes """
    # get the hashes of the gene mers for the two nodes
    firstNodeHash, secondNodeHash = extract_node_hashes(firstNode,
                                                        secondNode)
    # decide which gene mer will be the source node and which the target node
    sourceNode, targetNode = sort_node_hashes(firstNodeHash,
                                            secondNodeHash)
    return sourceNode, targetNode

class Edge:
    def __init__(self,
                firstNode,
                secondNode):
        self.sourceNode, self.targetNode = define_source_and_target(firstNode,
                                                                    secondNode)
        self.edgeCoverage = 0
    def get_sourceNode(self) -> Node:
        """ return the assigned source Node object """
        return self.sourceNode
    def get_targetNode(self) -> Node:
        """ return the assigned target Node object """
        return self.targetNode
    def get_coverage(self) -> int:
        """ return the number of time this edge is seen in the data """
        return self.edgeCoverage
    def increment_coverage(self) -> int:
        """ increase the edge coverage by 1 and return the new edge coverage """
        self.edgeCoverage += 1
        return self.edgeCoverage
    def __eq__(self, otherEdge) -> bool:
        """ return a bool of whether the two edge objects are identical """
        return self == otherEdge
    def __hash__(self):
        """ return a hash of a tuple of the source and target gene mer hashes """
        nodeHashes = (self.sourceNode.__hash__(), self.targetNode.__hash__())
        return hash(nodeHashes)