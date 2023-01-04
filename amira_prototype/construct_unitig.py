from construct_graph import GeneMerGraph
from construct_node import Node
from construct_gene import Gene

class Unitig:
    def __init__(self,
                graph: GeneMerGraph,
                listOfGenes: list):
            self._graph = graph
            self._listOfGenes = listOfGenes
    def get_graph(self):
        """ returns the geneMerGraph """
        return self._graph
    def get_selected_genes(self):
        """ returns the list of selected genes """
        return self._listOfGenes
    def get_nodes_of_interest(self):
        """ extracts the graph nodes containing the genes of interest and returns them as a list """
        return self.get_graph().get_nodes_containing(self.get_selected_genes())
    def get_forward_node_from_node(self,
                                sourceNode: Node) -> list:
        """ returns a list of nodes in the forward direction from this node until a branch or end of unitig is reached """
        nodeForwardEdges = sourceNode.get_forward_edge_hashes()
        if len(nodeForwardEdges) > 0:
            for edge_hash in nodeForwardEdges:
                edge = self.get_graph().get_edges()[edge_hash]
                targetNode = edge.get_targetNode()
                targetNodeDegree = self.get_graph().get_degree(targetNode)
                if targetNodeDegree == 2 or targetNodeDegree == 1:
                    targetNodeDirection = edge.get_targetNodeDirection()
                    return True, targetNode, targetNodeDirection
                else:
                    return False, None, None
        else:
            return False, None, None
    def get_backward_node_from_node(self,
                                sourceNode: Node) -> list:
        """ returns a list of nodes in the forward direction from this node until a branch or end of unitig is reached """
        nodeBackwardEdges = sourceNode.get_backward_edge_hashes()
        if len(nodeBackwardEdges) > 0:
            for edge_hash in nodeBackwardEdges:
                edge = self.get_graph().get_edges()[edge_hash]
                targetNode = edge.get_targetNode()
                targetNodeDegree = self.get_graph().get_degree(targetNode)
                if targetNodeDegree == 2 or targetNodeDegree == 1:
                    targetNodeDirection = edge.get_targetNodeDirection()
                    return True, targetNode, targetNodeDirection
                else:
                    return False, None, None
        else:
            return False, None, None
    def get_forward_path_from_node(self,
                                node):
        forward_nodes_from_node = []
        forwardExtend, forwardNode, forwardNodeDirection = self.get_forward_node_from_node(node)
        while forwardExtend:
            forward_nodes_from_node.append(forwardNode.get_canonical_geneMer())
            if forwardNodeDirection == 1:
                forwardExtend, forwardNode, forwardNodeDirection = self.get_forward_node_from_node(forwardNode)
            else:
                forwardExtend, forwardNode, forwardNodeDirection = self.get_backward_node_from_node(forwardNode)
        return forward_nodes_from_node
    def get_backward_path_from_node(self,
                                    node):
        backward_nodes_from_node = []
        backwardExtend, backwardNode, backwardNodeDirection = self.get_backward_node_from_node(node)
        while backwardExtend:
            backward_nodes_from_node.append(backwardNode.get_canonical_geneMer())
            if backwardNodeDirection == -1:
                backwardExtend, backwardNode, backwardNodeDirection = self.get_backward_node_from_node(backwardNode)
            else:
                backwardExtend, backwardNode, backwardNodeDirection = self.get_forward_node_from_node(backwardNode)
        return list(reversed(backward_nodes_from_node))
    def get_unitig_for_node(self,
                            node):
        """ builds a unitig starting from the node of interest and expanding in both directions """
        if self.get_graph().get_degree(node) == 2:
            forward_nodes_from_node = self.get_forward_path_from_node(node)
            backward_nodes_from_node = self.get_backward_path_from_node(node)
            unitig = backward_nodes_from_node + [node.get_canonical_geneMer()] + forward_nodes_from_node
            return unitig
    def get_unitigs_of_interest(self):
        nodesOfInterest = self.get_nodes_of_interest()
        all_unitigs = []
        for node in nodesOfInterest:
            node_unitig = self.get_unitig_for_node(node)
            if node_unitig:
                node_unitig = [[g.get_name() for g in geneMer] for geneMer in node_unitig]
                all_unitigs.append(node_unitig)
        return all_unitigs
    def get_path(self) -> list():
        """ returns the linear path of nodes for the genes of interest """
        return
    def contains_node(self,
                    node: Node) -> bool:
        """ returns True if it contains the given node (i.e. it is the unitig for this node) """
        return

class UnitigBuilder:
    def get_unitigs(self,
                amr_genes: list()) -> list():
        """ returns all unitigs containing the specified AMR genes  """
        return
    def visualise_unitigs(self,
                        unitigs: list(),
                        output_directory: str):
        """ generate a figure to visualise the genes, order of genes, direction and lengths of genes on unitigs containing AMR genes """
        return