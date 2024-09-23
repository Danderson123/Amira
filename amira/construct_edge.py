from amira.construct_gene import hashlib_hash
from amira.construct_node import Node


def extract_node_hashes(firstNode, secondNode):
    """get the gene mer hashes for the first node and the second node"""
    firstNodeHash = firstNode.__hash__()
    secondNodeHash = secondNode.__hash__()
    return firstNodeHash, secondNodeHash


def sort_node_hashes(firstNodeHash, secondNodeHash):
    """sort the hashes in ascending order"""
    sortedHashes = sorted([firstNodeHash, secondNodeHash])
    # the sourceNode is at index 0
    sourceNode = sortedHashes[0]
    # the targetNode is at index 1
    targetNode = sortedHashes[1]
    return sourceNode, targetNode


def define_source_and_target(firstNode, secondNode):
    """assigns and returns the source node and target node"""
    # get the hashes of the gene mers for the two nodes
    firstNodeHash, secondNodeHash = extract_node_hashes(firstNode, secondNode)
    # decide which gene mer will be the source node and which the target node
    sourceNode, targetNode = sort_node_hashes(firstNodeHash, secondNodeHash)
    return sourceNode, targetNode


class Edge:
    def __init__(self, sourceNode, targetNode, sourceNodeDirection, targetNodeDirection):
        self.sourceNode = sourceNode
        self.targetNode = targetNode
        self.edgeCoverage = 0
        self.sourceNodeDirection = sourceNodeDirection
        self.targetNodeDirection = targetNodeDirection

    def get_sourceNode(self) -> Node:
        """return the assigned source Node object"""
        return self.sourceNode

    def get_targetNode(self) -> Node:
        """return the assigned target Node object"""
        return self.targetNode

    def set_sourceNode(self, new_sourceNode) -> Node:
        """return the assigned source Node object"""
        self.sourceNode = new_sourceNode
        return self.get_sourceNode()

    def set_targetNode(self, new_targetNode) -> Node:
        """return the assigned target Node object"""
        self.targetNode = new_targetNode
        return self.get_targetNode()

    def set_sourceNodeDirection(self, sourceDirection) -> int:
        """return and int of the new direction of the source node"""
        self.sourceNodeDirection = sourceDirection
        return self.get_sourceNodeDirection()

    def get_sourceNodeDirection(self) -> int:
        """modify and return and int of the direction of the source node"""
        return self.sourceNodeDirection

    def set_targetNodeDirection(self, targetDirection) -> int:
        """modify and return and int of the direction of the target node"""
        self.targetNodeDirection = targetDirection
        return self.get_targetNodeDirection()

    def get_targetNodeDirection(self) -> int:
        """return an int of the direction of the target node"""
        return self.targetNodeDirection

    def get_edge_coverage(self) -> int:
        """return the number of time this edge is seen in the data"""
        return self.edgeCoverage

    def increment_edge_coverage(self) -> int:
        """increase the edge coverage by 1 and return the new edge coverage"""
        self.edgeCoverage += 1
        return self.get_edge_coverage()

    def extend_edge_coverage(self, value) -> int:
        """increase the edge coverage by the value and return the new edge coverage"""
        self.edgeCoverage += value
        return self.get_edge_coverage()

    def reduce_edge_coverage(self):
        """reduce the edge coverage by 1 and return the new edge coverage"""
        self.edgeCoverage -= 1
        return self.get_edge_coverage()

    def __eq__(self, otherEdge) -> bool:
        """return a bool of whether the two edges connect the same nodes"""
        sortedEdgeHash = tuple(
            sorted([self.get_sourceNode().__hash__(), self.get_targetNode().__hash__()])
        )
        otherSortedEdgeHash = tuple(
            sorted([otherEdge.get_sourceNode().__hash__(), otherEdge.get_targetNode().__hash__()])
        )
        return sortedEdgeHash == otherSortedEdgeHash

    def __hash__(self):
        """return a hash of a tuple of the source and target gene mer hashes"""
        forwardEdgeHash = hashlib_hash(
            (
                (
                    self.get_sourceNode().__hash__() * self.get_sourceNodeDirection(),
                    self.get_targetNode().__hash__() * self.get_targetNodeDirection(),
                )
            )
        )
        reverseEdgeHash = hashlib_hash(
            (
                (
                    self.get_sourceNode().__hash__() * self.get_sourceNodeDirection() * -1,
                    self.get_targetNode().__hash__() * self.get_targetNodeDirection() * -1,
                )
            )
        )
        sortedEdgeHashes = sorted([forwardEdgeHash, reverseEdgeHash])
        edgeHash = sortedEdgeHashes[0]
        return edgeHash
