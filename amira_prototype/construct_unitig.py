from construct_graph import GeneMerGraph
from construct_node import Node

class Unitigs:
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
    def get_nodes_of_interest(self,
                            geneOfInterest):
        """ extracts the graph nodes containing the genes of interest and returns them as a list """
        return self.get_graph().get_nodes_containing(geneOfInterest)
    def get_forward_node_from_node(self,
                                sourceNode: Node) -> list:
        """ returns a list of nodes in the forward direction from this node until a branch or end of unitig is reached """
        # get the list of forward edge hashes for this node
        nodeForwardEdges = sourceNode.get_forward_edge_hashes()
        if len(nodeForwardEdges) > 0:
            # iterate through the edge hashes
            for edge_hash in nodeForwardEdges:
                # get the edge object corresponding to this edge hash
                edge = self.get_graph().get_edges()[edge_hash]
                # get the target node for this edge
                targetNode = edge.get_targetNode()
                # get the degree of the target node
                targetNodeDegree = self.get_graph().get_degree(targetNode)
                # if the degree of the target node is 1 or 2 we can extend the linear path to the next node
                if targetNodeDegree == 2 or targetNodeDegree == 1:
                    # get the direction we are going into the target node
                    targetNodeDirection = edge.get_targetNodeDirection()
                    return True, targetNode, targetNodeDirection
                # else we cannot extend the linear path
                else:
                    return False, None, None
        else:
            return False, None, None
    def get_backward_node_from_node(self,
                                sourceNode: Node) -> list:
        """ returns a list of nodes in the forward direction from this node until a branch or end of unitig is reached """
        # get the list of forward edge hashes for this node
        nodeBackwardEdges = sourceNode.get_backward_edge_hashes()
        if len(nodeBackwardEdges) > 0:
            # iterate through the edge hashes
            for edge_hash in nodeBackwardEdges:
                # get the edge object corresponding to this edge hash
                edge = self.get_graph().get_edges()[edge_hash]
                # get the target node for this edge
                targetNode = edge.get_targetNode()
                # get the degree of the target node
                targetNodeDegree = self.get_graph().get_degree(targetNode)
                # if the degree of the target node is 1 or 2 we can extend the linear path to the next node
                if targetNodeDegree == 2 or targetNodeDegree == 1:
                    # get the direction we are going into the target node
                    targetNodeDirection = edge.get_targetNodeDirection()
                    return True, targetNode, targetNodeDirection
                else:
                    # else we cannot extend the linear path
                    return False, None, None
        else:
            return False, None, None
    def get_forward_path_from_node(self,
                                node):
        forward_nodes_from_node = []
        # get the next node in the forward direction
        forwardExtend, forwardNode, forwardNodeDirection = self.get_forward_node_from_node(node)
        while forwardExtend:
            # if we are extending further in the forward direction, get the next canonical gene mer
            forward_nodes_from_node.append(forwardNode.get_canonical_geneMer())
            # if we enter the next node in the forward direction, we get the next forward node
            if forwardNodeDirection == 1:
                forwardExtend, forwardNode, forwardNodeDirection = self.get_forward_node_from_node(forwardNode)
            # if we enter the next node in the backward direction, we get the next bacckward node
            else:
                forwardExtend, forwardNode, forwardNodeDirection = self.get_backward_node_from_node(forwardNode)
        return forward_nodes_from_node
    def get_backward_path_from_node(self,
                                    node):
        backward_nodes_from_node = []
        # get the next node in the backward direction
        backwardExtend, backwardNode, backwardNodeDirection = self.get_backward_node_from_node(node)
        while backwardExtend:
            backward_nodes_from_node.append(backwardNode.get_canonical_geneMer())
            # if we are extending further in the backward direction, get the next canonical gene mer
            if backwardNodeDirection == -1:
                backwardExtend, backwardNode, backwardNodeDirection = self.get_backward_node_from_node(backwardNode)
            # if we enter the next node in the forward direction, we get the next forward node
            else:
                backwardExtend, backwardNode, backwardNodeDirection = self.get_forward_node_from_node(backwardNode)
        return list(reversed(backward_nodes_from_node))
    def get_unitig_for_node(self,
                            node):
        """ builds a unitig starting from the node of interest and expanding in both directions """
        if self.get_graph().get_degree(node) == 2 or self.get_graph().get_degree(node) == 1:
            # get the forward nodes from this node
            forward_nodes_from_node = self.get_forward_path_from_node(node)
            # get the backward nodes from this node
            backward_nodes_from_node = self.get_backward_path_from_node(node)
            # join the backward nodes, this node, and the forward nodes to get the path of nodes
            unitig = backward_nodes_from_node + [node.get_canonical_geneMer()] + forward_nodes_from_node
            return unitig
    def get_unitigs_of_interest(self):
        """ returns a dictionary of genes of interest and their linear paths in the graph """
        unitigsOfInterest = {}
        # iterate through the list of specified genes
        for geneOfInterest in self.get_selected_genes():
            # get the graph nodes containing this gene
            nodesOfInterest = self.get_nodes_of_interest(geneOfInterest)
            all_unitigs = []
            # iterate through the nodes containing this gene
            for node in nodesOfInterest:
                # get the linear path for this node
                node_unitig = self.get_unitig_for_node(node)
                if node_unitig:
                    if not (node_unitig in all_unitigs or list(reversed(node_unitig)) in all_unitigs):
                        all_unitigs.append(node_unitig)
            # populate a dictionary with the gene name and the corresponding unitigs
            unitigsOfInterest[geneOfInterest] = all_unitigs
        return unitigsOfInterest

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
