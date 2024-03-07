import os
import statistics
import sys
from collections import Counter, deque
from itertools import product
import numpy as np
from tqdm import tqdm

from amira_prototype.construct_edge import Edge
from amira_prototype.construct_gene import convert_int_strand_to_string, Gene
from amira_prototype.construct_gene_mer import GeneMer
from amira_prototype.construct_node import Node
from amira_prototype.construct_read import Read

sys.setrecursionlimit(50000)


def subsample_to_minimizers(kmers, window_size=5):
    minimizers = []
    for i in range(len(kmers) - window_size + 1):
        window = kmers[i : i + window_size]
        minimizer = sorted(window, key=hash)[0]
        if minimizer not in minimizers:
            minimizers.append(minimizer)
    return minimizers


class GeneMerGraph:
    def __init__(self, readDict, kmerSize):
        self._reads = readDict
        self._kmerSize = kmerSize
        self._minNodeCoverage = 1
        self._minEdgeCoverage = 1
        self._nodes = {}
        self._edges = {}
        self._readNodes = {}
        self._readNodeDirections = {}
        self._shortReads = {}
        self._readsToCorrect = set()
        # initialise the graph
        for readId in tqdm(self.get_reads()):
            read = Read(readId, self.get_reads()[readId])
            # get the gene mers on this read
            geneMers = [g for g in read.get_geneMers(self.get_kmerSize())]
            # if there are no gene mers there is nothing to add to the graph
            if len(geneMers) == 0:
                self._shortReads[readId] = read.get_annotatedGenes()
                continue
            if not len(geneMers) == 1:
                # iterate through all the gene-mers except the last
                for g in range(len(geneMers) - 1):
                    # get the read id for this read
                    readId = read.get_readId()
                    # add a source node to the graph
                    sourceNode = self.add_node(geneMers[g], [read])
                    # add the sourceNode to the read
                    self.add_node_to_read(sourceNode, readId, geneMers[g].get_geneMerDirection())
                    # increase the source node coverage by 1
                    sourceNode.increment_node_coverage()
                    # add the target node to the graph
                    targetNode = self.add_node(geneMers[g + 1], [read])
                    # add an edge from source to target and target to source
                    sourceToTargetEdge, reverseTargetToSourceEdge = self.add_edge(
                        geneMers[g], geneMers[g + 1]
                    )
                    # increment the edge coverages by 1
                    sourceToTargetEdge.increment_edge_coverage()
                    reverseTargetToSourceEdge.increment_edge_coverage()
                targetNodeHash = geneMers[-1].__hash__()
                targetNode = self.get_node_by_hash(targetNodeHash)
                # increment the coverage of the target if it is the last gene mer in the read
                targetNode.increment_node_coverage()
                self.add_node_to_read(targetNode, readId, geneMers[-1].get_geneMerDirection())
            else:
                # add a single node to the graph if there is only 1 gene mer
                readId = read.get_readId()
                sourceNode = self.add_node(geneMers[0], [read])
                self.add_node_to_read(sourceNode, readId, geneMers[0].get_geneMerDirection())
                sourceNode.increment_node_coverage()
        # assign component Ids to all nodes in the graph
        self.assign_component_ids()

    def get_reads(self) -> dict[str, list[str]]:
        """return a dictionary of all reads and their genes"""
        return self._reads

    def get_readNodes(self) -> dict[str, list[str]]:
        """return a dictionary of all reads and their node hashes"""
        return self._readNodes

    def get_readNodeDirections(self) -> dict[str, list[int]]:
        """return a dictionary of all reads and their node directions"""
        return self._readNodeDirections

    def get_kmerSize(self) -> int:
        """return an integer of the gene-mer size"""
        return self._kmerSize

    def get_minEdgeCoverage(self) -> int:
        """return an integer of the minimum coverage for edges"""
        return self._minEdgeCoverage

    def get_nodes(self) -> dict[int, Node]:
        """return the node dictionary"""
        return self._nodes

    def get_edges(self) -> dict[int, Edge]:
        """return the edge dictionary"""
        return self._edges

    def get_minNodeCoverage(self) -> int:
        """return the minimum node coverage"""
        return self._minNodeCoverage

    def get_reads_to_correct(self) -> set:
        """return a set of reads that need to be corrected"""
        return self._readsToCorrect

    def all_nodes(self) -> list[Node]:
        """return a generator for all nodes in the graph and their attributes"""
        for nodeHash in self.get_nodes():
            yield self.get_nodes()[nodeHash]

    def get_reads_for_nodes(self, list_of_nodes: list[int]) -> set[str]:
        """return a set of the reads attached to a list of nodes"""
        reads = set()
        for node_hash in list_of_nodes:
            node = self.get_node_by_hash(node_hash)
            reads.update(node.get_list_of_reads())
        return reads

    def add_node_to_read(self, node: Node, readId: str, node_direction: int) -> list[int]:
        # add the read ID to the read dict if it is not present
        if readId not in self.get_readNodes():
            self.get_readNodes()[readId] = []
            self.get_readNodeDirections()[readId] = []
        # add the hash for the node as attributes for this read
        self.get_readNodes()[readId].append(node.__hash__())
        self.get_readNodeDirections()[readId].append(node_direction)
        # each node occurrence will occur in the list (including duplicates)
        return self.get_readNodes()[readId]

    def get_nodes_containing_read(self, readId: str) -> list[Node]:
        """return a list of nodes that contain a read of interest"""
        # get the node hashes that contain this read
        listOfNodeHashes = self.get_readNodes()[readId]
        # this function gets the node for a node hash if it has not been filtered from the graph
        listOfNodes = [self.get_node_by_hash(h) for h in listOfNodeHashes if h in self.get_nodes()]
        return listOfNodes

    def add_node_to_nodes(self, node: Node, nodeHash: int) -> None:
        """add a node to the dictionary of nodes"""
        self._nodes[nodeHash] = node

    def get_node_by_hash(self, nodeHash: int) -> Node:
        """return a node object corresponding to a node hash"""
        return self.get_nodes()[nodeHash]

    def add_node(self, geneMer: GeneMer, reads: list) -> Node:
        """
        Add a gene mer to the graph if it does not exist, else get the node.
        Returns the Node itself
        """
        nodeHash = geneMer.__hash__()
        if nodeHash not in self.get_nodes():
            # convert the gene-mer to a node
            node = Node(geneMer)
            # add the node to the graph
            self.add_node_to_nodes(node, nodeHash)
        else:
            node = self.get_node_by_hash(nodeHash)
            # assert geneMer.get_geneMerDirection() == node.get_geneMer().get_geneMerDirection()
        # add the reads for this node to the dictionary of reads
        for r in reads:
            node.add_read(r.get_readId())
        return node

    def get_node(self, geneMer: GeneMer) -> Node:
        """get a node given a GeneMer"""
        # the gene-mer and node hashes are the same
        nodeHash = geneMer.__hash__()
        # confirm the node is in the node dictionary
        assert nodeHash in self.get_nodes(), "This gene-mer is not in the graph"
        # return the node
        return self.get_nodes()[nodeHash]

    def get_nodes_containing(self, geneOfInterest: str) -> list[Node]:
        """
        Return all nodes that contain a given gene.
        This is useful for later, when we want to get nodes that contain AMR genes.
        """
        # ensure there is no strand information present on the requested genes
        assert not (
            geneOfInterest[0] == "+" or geneOfInterest[0] == "-"
        ), "Strand information cannot be present for any specified genes"
        assert isinstance(geneOfInterest, str), "Gene of interest is the wrong type"
        selectedNodes = []
        # iterate through the nodes in the graph
        for node in self.all_nodes():
            # get the canonical geneMer for this node
            geneMer = node.get_canonical_geneMer()
            # get the genes in this gene mer
            geneNames = [g.get_name() for g in geneMer]
            # add node to the selectedNodes list if it contains any of the specified genes
            if geneOfInterest in geneNames:
                # add the node to the list of selected nodes
                selectedNodes.append(node)
        return selectedNodes

    def create_edges(
        self,
        sourceNode: Node,
        targetNode: Node,
        sourceGeneMerDirection: int,
        targetGeneMerDirection: int,
    ):
        """Create and return the forward and reverse edges."""
        # the edge that will go from the source to target
        sourceToTargetEdge = Edge(
            sourceNode, targetNode, sourceGeneMerDirection, targetGeneMerDirection
        )
        # the edge that will go from the reverse source to reverse target
        reverseTargetToSourceEdge = Edge(
            targetNode, sourceNode, targetGeneMerDirection * -1, sourceGeneMerDirection * -1
        )
        return sourceToTargetEdge, reverseTargetToSourceEdge

    def get_edge_by_hash(self, edgeHash: int) -> Edge:
        """return an edge object corresponding to a edge hash"""
        return self.get_edges()[edgeHash]

    def add_edge_to_edges(self, edge: Edge) -> Edge:
        """Add the edge to the graph if not present, else increment the coverage by 1.
        Returns the edge.
        """
        # get the hash for the edge we want
        edgeHash = edge.__hash__()
        # see if the edge is in the edge dictionary
        if edgeHash not in self.get_edges():
            self._edges[edgeHash] = edge
        return self.get_edge_by_hash(edgeHash)

    def add_edges_to_graph(
        self, sourceToTargetEdge: Edge, reverseTargetToSourceEdge: Edge
    ) -> tuple:
        """add the edges to the dictionary of edges and return the updated edges"""
        sourceToTargetEdge = self.add_edge_to_edges(sourceToTargetEdge)
        reverseTargetToSourceEdge = self.add_edge_to_edges(reverseTargetToSourceEdge)
        return sourceToTargetEdge, reverseTargetToSourceEdge

    def add_edge_to_node(self, node: Node, edge: Edge) -> Node:
        """
        Add the edge and reverse edge hashes to the node attributes.
        Return the forward and reverse edge hashes for this node.
        """
        if edge.get_sourceNodeDirection() == 1:
            edgeHash = edge.__hash__()
            node.add_forward_edge_hash(edgeHash)
        if edge.get_sourceNodeDirection() == -1:
            edgeHash = edge.__hash__()
            node.add_backward_edge_hash(edgeHash)
        return node

    def add_edge(self, sourceGeneMer: GeneMer, targetGeneMer: GeneMer) -> tuple:
        """
        Add a forward and backward edge if they do not exist, else increment the edge coverage by 1.
        Also add the nodes to the graph if they are not present already.
        Returns the updated source and target lists of forward and reverse edges.
        """
        # if the source is not in the graph then add it, else get the node
        sourceNode = self.add_node(sourceGeneMer, [])
        # if the target is not in the graph then add it, else get the node
        targetNode = self.add_node(targetGeneMer, [])
        # create the edge objects
        sourceToTargetEdge, reverseTargetToSourceEdge = self.create_edges(
            sourceNode,
            targetNode,
            sourceGeneMer.get_geneMerDirection(),
            targetGeneMer.get_geneMerDirection(),
        )
        # add the edges to the graph
        sourceToTargetEdge, reverseTargetToSourceEdge = self.add_edges_to_graph(
            sourceToTargetEdge, reverseTargetToSourceEdge
        )
        # add the hashes of the edges to the source and target nodes
        sourceNode = self.add_edge_to_node(sourceNode, sourceToTargetEdge)
        targetNode = self.add_edge_to_node(targetNode, reverseTargetToSourceEdge)
        return sourceToTargetEdge, reverseTargetToSourceEdge

    def get_degree(self, node: Node) -> int:
        """return an integer of the number of neighbours for this node"""
        degrees = len(node.get_forward_edge_hashes()) + len(node.get_backward_edge_hashes())
        return degrees

    def get_forward_neighbors(self, node: Node) -> list:
        """return a list of nodes corresponding to the forward neighbors for this node"""
        return [edge.get_targetNode() for edge in self.get_forward_edges(node)]

    def get_backward_neighbors(self, node: Node) -> list:
        """return a list of nodes corresponding to the backward neighbors for this node"""
        return [edge.get_targetNode() for edge in self.get_backward_edges(node)]

    def get_all_neighbors(self, node: Node) -> list[Node]:
        """return a list of combined forward and reverse neighbor nodes for this node"""
        return [
            neighbor
            for neighbor in self.get_forward_neighbors(node) + self.get_backward_neighbors(node)
        ]

    def get_all_neighbor_hashes(self, node: Node) -> set[int]:
        """return a set of combined forward and reverse node hashes for this node"""
        return set([neighbor.__hash__() for neighbor in self.get_all_neighbors(node)])

    def get_forward_edges(self, node: Node) -> list[Edge]:
        """return a list of integers of node identifiers connected to this node by a forward edge"""
        return [self.get_edge_by_hash(edgeHash) for edgeHash in node.get_forward_edge_hashes()]

    def get_backward_edges(self, node: Node) -> list[Edge]:
        """return a list of backward node identifiers"""
        return [self.get_edge_by_hash(edgeHash) for edgeHash in node.get_backward_edge_hashes()]

    def check_if_nodes_are_adjacent(self, sourceNode: Node, targetNode: Node) -> bool:
        """returns a bool of if the sourceNode and targetNode are neighbors"""
        return targetNode.__hash__() in self.get_all_neighbor_hashes(
            sourceNode
        ) and sourceNode.__hash__() in self.get_all_neighbor_hashes(targetNode)

    def get_edge_hashes_between_nodes(self, sourceNode: Node, targetNode: Node) -> tuple:
        """return a tuple of the source to target edge hashes"""
        # check that the two nodes are adjacent
        assert self.check_if_nodes_are_adjacent(sourceNode, targetNode)
        # get the edge hash from source to target
        sourceNodeEdges = []
        for edge in self.get_forward_edges(sourceNode) + self.get_backward_edges(sourceNode):
            if edge.get_targetNode() == targetNode:
                sourceNodeEdgeHash = edge.__hash__()
                sourceNodeEdges.append(sourceNodeEdgeHash)
        # get the edge hash from target to source
        targetNodeEdges = []
        for edge in self.get_forward_edges(targetNode) + self.get_backward_edges(targetNode):
            if edge.get_targetNode() == sourceNode:
                targetNodeEdgeHash = edge.__hash__()
                targetNodeEdges.append(targetNodeEdgeHash)
        assert not (
            sourceNodeEdgeHash == [] or targetNodeEdgeHash == []
        ), "There are edges missing from the source and target nodes"
        if not (len(sourceNodeEdges) > 1 or len(targetNodeEdges) > 1):
            return (sourceNodeEdges[0], targetNodeEdges[0])
        else:
            return (sourceNodeEdges, targetNodeEdges)

    def get_edges_between_nodes(self, sourceNode: Node, targetNode: Node) -> list[Edge]:
        # get the edge hashes between these nodes
        sourceToTargetHash, targetToSourceHash = self.get_edge_hashes_between_nodes(
            sourceNode, targetNode
        )
        if not (isinstance(sourceToTargetHash, list) or isinstance(targetToSourceHash, list)):
            return self.get_edge_by_hash(sourceToTargetHash), self.get_edge_by_hash(
                targetToSourceHash
            )
        else:
            return [self.get_edge_by_hash(h) for h in sourceToTargetHash], [
                self.get_edge_by_hash(h) for h in targetToSourceHash
            ]

    def remove_edge_from_edges(self, edgeHash: int) -> None:
        """remove an edge object from the dictionary of all edges by hash"""
        del self._edges[edgeHash]
        assert (
            edgeHash not in self.get_edges()
        ), "This edge was not removed from the graph successfully"

    def remove_edge(self, edgeHash: int) -> None:
        """remove an edge from all graph edges and the node edges"""
        # check that the edge exists in the graph
        # assert edgeHash in self.get_edges(), "This edge is not in the graph"
        if not edgeHash in self.get_edges():
            return
        # get the edge object from all graph edges by hash
        edge = self.get_edges()[edgeHash]
        # get the source node object of the edge
        sourceNode = edge.get_sourceNode()
        # get the source node direction of the edge
        sourceNodeDirection = edge.get_sourceNodeDirection()
        # if the direction is forward then remove the edge from the forward edges for the node
        if sourceNodeDirection == 1:
            sourceNode.remove_forward_edge_hash(edgeHash)
        # if the direction is backward then remove the edge from the backward edges for the node
        if sourceNodeDirection == -1:
            sourceNode.remove_backward_edge_hash(edgeHash)
        # remove the edge object from all of the edges in the graph
        del self._edges[edgeHash]

    def get_total_number_of_nodes(self) -> int:
        """return an integer of the total number of nodes in the filtered graph"""
        return len(self.get_nodes())

    def get_total_number_of_edges(self) -> int:
        """return an integer of the total number of edges in the filtered graph"""
        return len(self.get_edges())

    def get_total_number_of_reads(self) -> int:
        """return an integer of the total number of reads contributing to the filtered graph"""
        return len(self.get_reads())

    def remove_node_from_reads(self, node_to_remove: Node) -> None:
        """remove a node from the list of nodes and return the modified node list"""
        for readId in node_to_remove.get_reads():
            new_nodes = []
            new_directions = []
            for i in range(len(self.get_readNodes()[readId])):
                if self.get_readNodes()[readId][i] != node_to_remove.__hash__():
                    new_nodes.append(self.get_readNodes()[readId][i])
                    new_directions.append(self.get_readNodeDirections()[readId][i])
                else:
                    new_nodes.append(None)
                    new_directions.append(None)
            # keep track of all of the reads we have modified
            self.get_reads_to_correct().add(readId)
            self.get_readNodes()[readId] = new_nodes
            self.get_readNodeDirections()[readId] = new_directions

    def remove_node(self, node: Node):
        """remove a node from the graph and all of its edges"""
        # get the node hash
        nodeHash = node.__hash__()
        # confirm the node to remove is in the node dictionary
        assert nodeHash in self.get_nodes(), "This node is not in the graph"
        assert node == self.get_node_by_hash(nodeHash)
        # remove the nodeHash from the readNodes
        self.remove_node_from_reads(node)
        # get the forward edge hashes of the node to remove
        forward_edges_to_remove = node.get_forward_edge_hashes()
        # get the backward edge hashes of the node to remove
        backward_edges_to_remove = node.get_backward_edge_hashes()
        # remove the forward and backward edges from the edge dictionary
        for edgeHash in set(forward_edges_to_remove + backward_edges_to_remove):
            targetNode = self.get_edge_by_hash(edgeHash).get_targetNode()
            edgesToTarget = self.get_edge_hashes_between_nodes(node, targetNode)
            for e in edgesToTarget:
                # confirm the edge is in the graph
                self.remove_edge(e)
        # remove the node from the node dictionary
        del self.get_nodes()[nodeHash]

    def set_minNodeCoverage(self, minNodeCoverage: int):
        """set the minimum node coverage for the graph and return the new node coverage"""
        self._minNodeCoverage = minNodeCoverage
        return self.get_minNodeCoverage()

    def set_minEdgeCoverage(self, minEdgeCoverage: int):
        """set the minimum edge coverage for the graph and return the new edge coverage"""
        self._minEdgeCoverage = minEdgeCoverage
        return self.get_minEdgeCoverage()

    def list_nodes_to_remove(self, minNodeCoverage: int) -> list:
        """returns a list of node objects to be filtered from the graph"""
        nodesToRemove = set()
        for nodeHash in self.get_nodes():
            # mark a node for removal if the coverage is less than the specified coverage
            if not self.get_nodes()[nodeHash].get_node_coverage() > minNodeCoverage - 1:
                nodesToRemove.add(self.get_nodes()[nodeHash])
        return nodesToRemove

    def list_edges_to_remove(self, minEdgeCoverage: int, nodesToRemove: list) -> list:
        """returns a list of edge hashes to be filtered from the graph"""
        edgesToRemove = set()
        for edgeHash in self.get_edges():
            # mark an edge for removal if the coverage is less than the specified coverage
            if not self.get_edges()[edgeHash].get_edge_coverage() > minEdgeCoverage - 1:
                edgesToRemove.add(edgeHash)
            # mark an edge for removal if either of the nodes are marked for removal
            if any(
                n in nodesToRemove
                for n in [
                    self.get_edges()[edgeHash].get_sourceNode(),
                    self.get_edges()[edgeHash].get_targetNode(),
                ]
            ):
                edgesToRemove.add(edgeHash)
        return edgesToRemove

    def filter_graph(self, minNodeCoverage: int, minEdgeCoverage: int):
        """filters the nodes and edges in the graph and then returns the filtered graph"""
        # set the new node coverage
        minNodeCoverage = self.set_minNodeCoverage(minNodeCoverage)
        # set the new edge coverage
        minEdgeCoverage = self.set_minEdgeCoverage(minEdgeCoverage)
        # filter the nodes
        nodesToRemove = self.list_nodes_to_remove(minNodeCoverage)
        # filter the edges
        edgesToRemove = self.list_edges_to_remove(minEdgeCoverage, nodesToRemove)
        # remove the edges
        for e in edgesToRemove:
            self.remove_edge(e)
        # remove the nodes
        for n in nodesToRemove:
            # add an edge between the neighbors of the removed node
            self.remove_node(n)
        return self

    def write_node_entry(self, node_id, node_string, node_coverage, reads, component_ID, nodeColor):
        """return a string of a gml node entry"""
        node_entry = "\tnode\t[\n"
        node_entry += "\t\tid\t" + str(node_id) + "\n"
        node_entry += '\t\tlabel\t"' + node_string + '"\n'
        node_entry += "\t\tcoverage\t" + str(node_coverage) + "\n"
        if component_ID:
            node_entry += "\t\tcomponent\t" + str(component_ID) + "\n"
        node_entry += '\t\treads\t"' + ",".join(reads) + '"\n'
        if nodeColor:
            node_entry += '\t\tcolor\t"' + str(nodeColor) + '"\n'
        node_entry += "\t]"
        return node_entry

    def write_edge_entry(self, source_node, target_node, edge_direction, edge_coverage):
        """return a string of a gml edge entry"""
        edge_entry = "\tedge\t[\n"
        edge_entry += "\t\tsource\t" + str(source_node) + "\n"
        edge_entry += "\t\ttarget\t" + str(target_node) + "\n"
        edge_entry += "\t\tdirection\t" + str(edge_direction) + "\n"
        edge_entry += "\t\tweight\t" + str(edge_coverage) + "\n"
        edge_entry += "\t]"
        return edge_entry

    def assign_Id_to_nodes(self):
        """assign an integer ID to each node in the graph"""
        currentNodeId = 0
        for node in self.all_nodes():
            assignedId = node.assign_node_Id(currentNodeId)
            assert assignedId == currentNodeId, "This node was assigned an incorrect ID"
            currentNodeId += 1

    def write_gml_to_file(self, output_file, gml_content):
        """Writes the gml to a file titled <output_file>.gml. Returns nothing."""
        # check if there is a directory to make
        if os.path.dirname(output_file) != "":
            # make the directory if it does not exist
            if not os.path.exists(os.path.dirname(output_file)):
                os.mkdir(os.path.dirname(output_file))
        # write the gml to the output file
        with open(output_file + ".gml", "w") as outGml:
            outGml.write("\n".join(gml_content))

    def get_gene_mer_genes(self, sourceNode: Node) -> list:
        """return a list of genes in the canonical gene mer and their strands"""
        # get the canonical gene mer for the node
        geneMer = sourceNode.get_canonical_geneMer()
        return [convert_int_strand_to_string(g.get_strand()) + g.get_name() for g in geneMer]

    def get_reverse_gene_mer_genes(self, sourceNode: Node) -> list:
        """return a list of genes in the canonical gene mer and their strands"""
        # get the reversed gene mer for the node
        geneMer = sourceNode.get_reverse_geneMer()
        return [convert_int_strand_to_string(g.get_strand()) + g.get_name() for g in geneMer]

    def get_gene_mer_label(self, sourceNode: Node) -> str:
        """returns a string of the genes in the canonical gene mer for a node and their strands"""
        # get a list of the gene strand + gene name for the canonical gene mer
        geneMerGenes = self.get_gene_mer_genes(sourceNode)
        # return a string of the gene mer genes and strands
        return "~~~".join(geneMerGenes)

    def get_nodes_with_degree(self, degree: int):
        """return a list of node objects with the specified degree"""
        assert isinstance(degree, int), "The input degree must be an integer."
        nodesWithDegree = []
        for node in self.all_nodes():
            nodeDegree = self.get_degree(node)
            if nodeDegree == degree:
                nodesWithDegree.append(node)
        return nodesWithDegree

    def insert_between_nodes_on_read(self, readNodes, preInsertion, postInsertion, node_to_insert):
        """returns a modified list with "*" marking site to insert node"""
        if not len(readNodes) == 1:
            newReadNodes = []
            i = 0
            while i < len(readNodes):
                newReadNodes.append(readNodes[i])
                if i + 1 < len(readNodes) and (
                    (readNodes[i] == preInsertion and readNodes[i + 1] == postInsertion)
                    or (readNodes[i] == postInsertion and readNodes[i + 1] == preInsertion)
                ):
                    newReadNodes.append(node_to_insert)
                i += 1
            # check if either of the nodes are at the ends but the node adjacent is missing
            if (readNodes[0] == preInsertion and not readNodes[1] == postInsertion) or (
                readNodes[0] == postInsertion and not readNodes[1] == preInsertion
            ):
                newReadNodes.insert(0, node_to_insert)
            if (readNodes[-1] == preInsertion and not readNodes[-2] == postInsertion) or (
                readNodes[-1] == postInsertion and not readNodes[-2] == preInsertion
            ):
                newReadNodes.append(node_to_insert)
        else:
            newReadNodes = readNodes
        return newReadNodes

    def make_replacement_dict(self, old_path, new_path):
        # insert a star to indicate where the new node will go
        assert not len(new_path) % 2 == 0, "The gene-mer size must be an odd number"
        if len(old_path) == 0:
            old_path.append("*")
        else:
            assert old_path[0] == new_path[0] and old_path[-1] == new_path[-1], [
                old_path[0],
                new_path[0],
                old_path[-1],
                new_path[-1],
            ]
        # If the list length is even, insert a star at the midpoint
        if len(old_path) % 2 == 0:
            midpoint = len(old_path) // 2
            old_path.insert(midpoint, "*")
        # make a dictionary that will be used to replace nodes later
        replacementDict = {}
        for n in range(len(old_path)):
            replacementDict[old_path[n]] = new_path[n]
        return replacementDict

    def modify_readNode_list(self, readNodes, replacementDict):
        modifiedReadNodes = readNodes[:]
        to_replace = list(replacementDict.keys())
        for i in range(len(to_replace)):
            if to_replace[i] == "*":
                preInsertion = to_replace[i - 1]
                postInsertion = to_replace[i + 1]
                # get the indices of node hashes immediately adjacent to the insertion site
                modifiedReadNodes = self.insert_between_nodes_on_read(
                    modifiedReadNodes, preInsertion, postInsertion, "*"
                )
        return modifiedReadNodes

    def get_genes_in_unitig(self, listOfNodes):
        if len(listOfNodes) == 1:
            return self.get_gene_mer_genes(self.get_node_by_hash(listOfNodes[0]))
        # store the new annotations in a list
        newAnnotations = []
        # iterate through the new readNodes
        for n in range(len(listOfNodes) - 1):
            sourceNode = self.get_node_by_hash(listOfNodes[n])
            targetNode = self.get_node_by_hash(listOfNodes[n + 1])
            # get the edge between the source and target node
            edgeHashes = self.get_edge_hashes_between_nodes(sourceNode, targetNode)
            edge = self.get_edge_by_hash(edgeHashes[0])
            # add either the forward or reverse gene mer depending on the node direction
            if n == 0:
                if edge.get_sourceNodeDirection() == 1:
                    newAnnotations += self.get_gene_mer_genes(sourceNode)
                else:
                    newAnnotations += self.get_reverse_gene_mer_genes(sourceNode)
            if edge.get_targetNodeDirection() == 1:
                newAnnotations.append(self.get_gene_mer_genes(targetNode)[-1])
            else:
                newAnnotations.append(self.get_reverse_gene_mer_genes(targetNode)[-1])
        return newAnnotations

    def old_follow_path_to_get_annotations(self, listOfNodes, read_ID, read_node_directions=None):
        # get the read node directions
        if not read_node_directions:
            read_node_directions = self.get_readNodeDirections()[read_ID]
        # return the gene-mer or rc gene-mer if there is one node
        if len(listOfNodes) == 1:
            geneMer_direction = read_node_directions[0]
            if geneMer_direction == 1:
                return self.get_gene_mer_genes(self.get_node_by_hash(listOfNodes[0]))
            else:
                return self.get_reverse_gene_mer_genes(self.get_node_by_hash(listOfNodes[0]))
        # store the new annotations in a list
        newAnnotations = []
        # iterate through the new readNodes
        assert listOfNodes[0], f"The first node in the list of nodes for read {read_ID} is None."
        assert len(listOfNodes) == len(read_node_directions)
        for n in range(len(listOfNodes)):
            # get the node we are adding
            node = self.get_node_by_hash(listOfNodes[n])
            # get the node direction
            node_direction = read_node_directions[n]
            # add all the gene-mer genes if this is the first node
            if n == 0:
                if node_direction == 1:
                    newAnnotations += self.get_gene_mer_genes(node)
                else:
                    newAnnotations += self.get_reverse_gene_mer_genes(node)
            else:
                if node_direction:
                    if node_direction == 1:
                        newAnnotations.append(self.get_gene_mer_genes(node)[-1])
                    else:
                        newAnnotations.append(self.get_reverse_gene_mer_genes(node)[-1])
                else:
                    newAnnotations.append(None)
        for i in range(len(newAnnotations)):
            if newAnnotations[i] is None:
                # Adjust the node_index calculation to avoid negative indices
                node_index = i - (self.get_kmerSize() - 1)
                # Ensure node_index is within bounds of listOfNodes
                node = self.get_node_by_hash(listOfNodes[node_index])
                last_k_minus_one_genes = newAnnotations[i - (self.get_kmerSize() - 1) : i]

                forward_gene_mer = self.get_gene_mer_genes(node)
                reverse_gene_mer = self.get_reverse_gene_mer_genes(node)

                if last_k_minus_one_genes == forward_gene_mer[: len(last_k_minus_one_genes)]:
                    newAnnotations[i] = forward_gene_mer[len(last_k_minus_one_genes)]
                elif last_k_minus_one_genes == reverse_gene_mer[: len(last_k_minus_one_genes)]:
                    newAnnotations[i] = reverse_gene_mer[len(last_k_minus_one_genes)]
                else:
                    # Consider handling this situation differently if it's not necessarily an error
                    print(
                        [self.get_gene_mer_label(self.get_node_by_hash(n)) for n in listOfNodes],
                        len(newAnnotations),
                    )
                    print(
                        f"Annotations: {newAnnotations}\nIndex: {i}\nLast k-1 Genes: {last_k_minus_one_genes}"
                    )
                    print(
                        f"Forward Gene Mer: {forward_gene_mer}\nReverse Gene Mer: {reverse_gene_mer}"
                    )
                    raise ValueError(
                        "Neither the forward nor reverse gene-mer overlaps with the existing gene calls."
                    )
        return newAnnotations

    def remove_short_linear_paths(self, min_length):
        """remove nodeHashes on reads if on a linear path of length < min_length. \
            Returns a list of nodeHashes that have bene removed"""
        paths_to_remove = []
        for node in self.all_nodes():
            if self.get_degree(node) == 1:
                path = self.get_linear_path_for_node(node)
                if len(path) > 0 and (
                    len(path) < min_length
                ):  # or all(self.get_node_by_hash(n).get_node_coverage() < 2 for n in path)):
                    paths_to_remove.append(path)
        removed = set()
        for path in paths_to_remove:
            for nodeHash in path:
                if nodeHash not in removed:
                    self.remove_node(self.get_node_by_hash(nodeHash))
                    removed.add(nodeHash)
        return list(removed)

    def get_forward_converging_paths(self, startNode):
        pair_paths = []
        # get the forward neighbors
        fowardNodes = self.get_forward_neighbors(startNode)
        # if there are 2 forward neighbors
        if len(fowardNodes) == 2:
            paths = []
            for fwNode in fowardNodes:
                if not self.get_degree(fwNode) > 2:
                    # get the edge between the forward neighbour and the start node
                    sourceToTargetEdge, targetToSourceEdge = self.get_edges_between_nodes(
                        startNode, fwNode
                    )
                    # get the linear path for this forward node
                    path = [startNode.__hash__()] + self.get_forward_path_from_node(
                        fwNode, sourceToTargetEdge.get_targetNodeDirection(), True
                    )
                else:
                    path = [startNode.__hash__(), fwNode.__hash__()]
                # if the path does not start with the start node we need to reverse it
                if not path[0] == startNode.__hash__():
                    path = list(reversed(path))
                assert path[0] == startNode.__hash__()
                # add the path and its mean coverage to a list
                path_edge_coverages = [
                    self.get_edges_between_nodes(
                        self.get_node_by_hash(path[i]), self.get_node_by_hash(path[i + 1])
                    )[0].get_edge_coverage()
                    for i in range(len(path) - 1)
                ]
                paths.append((path, statistics.mean(path_edge_coverages)))
            # sort the paths based on their mean coverage
            paths.sort(key=lambda x: x[1], reverse=False)
            # see if the paths converge to the same node, if so this is a bubble
            if len(paths) > 1:
                first_elements = [sublist[0][0] for sublist in paths]
                last_elements = [sublist[0][-1] for sublist in paths]
                if len(set(first_elements)) == 1 and len(set(last_elements)) == 1:
                    pair_paths.append(paths)
        return pair_paths

    def get_backward_converging_paths(self, startNode):
        pair_paths = []
        # get the backward neighbors
        backwardNodes = self.get_backward_neighbors(startNode)
        # if there are 2 backward neighbors
        if len(backwardNodes) == 2:
            paths = []
            for bwNode in backwardNodes:
                if not self.get_degree(bwNode) > 2:
                    # get the edge between the forward neighbour and the start node
                    sourceToTargetEdge, targetToSourceEdge = self.get_edges_between_nodes(
                        startNode, bwNode
                    )
                    # get the linear path for this forward node
                    path = self.get_backward_path_from_node(
                        bwNode, sourceToTargetEdge.get_targetNodeDirection(), True
                    ) + [startNode.__hash__()]
                else:
                    path = [bwNode.__hash__(), startNode.__hash__()]
                # if the path does not end with the end node we need to reverse it
                if not path[-1] == startNode.__hash__():
                    path = list(reversed(path))
                assert path[-1] == startNode.__hash__()
                # add the path and its mean coverage to a list
                path_edge_coverages = [
                    self.get_edges_between_nodes(
                        self.get_node_by_hash(path[i]), self.get_node_by_hash(path[i + 1])
                    )[0].get_edge_coverage()
                    for i in range(len(path) - 1)
                ]
                paths.append((path, statistics.mean(path_edge_coverages)))
            # sort the paths based on their mean coverage
            paths.sort(key=lambda x: x[1], reverse=False)
            # see if the paths converge to the same node, if so this is a bubble
            if len(paths) > 1:
                first_elements = [sublist[0][0] for sublist in paths]
                last_elements = [sublist[0][-1] for sublist in paths]
                if len(set(first_elements)) == 1 and len(set(last_elements)) == 1:
                    pair_paths.append(paths)
        return pair_paths

    def replace_nodes_on_read(self, modifiedReadNodes, replacementDict):
        corrected_list = []
        for nodeHash in modifiedReadNodes:
            if nodeHash not in replacementDict:
                corrected_list.append(nodeHash)
            else:
                corrected_list.append(replacementDict[nodeHash])
        return corrected_list

    def correct_incorrectly_found_gene(self, lowest_coverage_path, highest_coverage_path):
        assert not len(lowest_coverage_path) % 2 == 0, "The gene-mer size must be an odd number"
        # we have to remove the middle node from the lowest coverage path
        modified_lowest_coverage_path = lowest_coverage_path[:]
        node_index_to_delete = int((len(lowest_coverage_path) + 1) / 2) - 1
        del modified_lowest_coverage_path[node_index_to_delete]
        removed_node = lowest_coverage_path[node_index_to_delete]
        # make the replacement dictionary
        replacementDict = {}
        for i in range(len(modified_lowest_coverage_path)):
            replacementDict[modified_lowest_coverage_path[i]] = highest_coverage_path[i]
        # get the reads we are going to correct
        pathReads = set()
        for nodeHash in lowest_coverage_path:
            for readId in self.get_node_by_hash(nodeHash).get_reads():
                pathReads.add(readId)
        # correct the reads
        for readId in list(pathReads):
            modified_readNodes = [n for n in self.get_readNodes()[readId] if not n == removed_node]
            # replace the nodes on the read
            corrected_list = self.replace_nodes_on_read(modified_readNodes, replacementDict)

    def get_converging_paths(self, startNode):
        # get the forward neighbors
        fowardNodes = self.get_forward_neighbors(startNode)
        # get the backward neighbors
        backwardNodes = self.get_backward_neighbors(startNode)
        # if there are 2 forward neighbors
        paths = []
        for node in fowardNodes + backwardNodes:
            # get the linear path for this forward node
            path = self.get_linear_path_for_node(node, True)
            # if the path does not start with the start node we need to reverse it
            if not path[0] == startNode.__hash__():
                path = list(reversed(path))
            if path[0] == startNode.__hash__():
                # add the path and its mean coverage to a list
                paths.append(
                    (
                        path,
                        statistics.mean(
                            [self.get_node_by_hash(h).get_node_coverage() for h in path]
                        ),
                    )
                )
        # group the paths if they share a start and end node
        groups = {}
        for lst in paths:
            key = (lst[0][0], lst[0][-1])
            if key not in groups:
                groups[key] = []  # Initialize a new group
            groups[key].append(lst)
        to_delete = []
        for terminals in groups:
            # sort the paths by coverage if there are more than 1 of them
            if len(groups[terminals]) > 1:
                groups[terminals].sort(key=lambda x: x[1], reverse=False)
            # remove the paths if there is only one of them
            else:
                to_delete.append(terminals)
        for t in to_delete:
            del groups[t]
        return groups

    def get_triangles(self):
        more_than_two_neighbors = [
            n for n in self.get_nodes() if self.get_degree(self.get_node_by_hash(n)) > 2
        ]
        triangles = []
        seen_triangles = set()
        for node_hash in more_than_two_neighbors:
            start_node = self.get_node_by_hash(node_hash)
            for start_direction, get_edges in [
                (1, self.get_forward_edges),
                (-1, self.get_backward_edges),
            ]:
                for start_edge in get_edges(start_node):
                    next_node = start_edge.get_targetNode()
                    target_node_direction = start_edge.get_targetNodeDirection()
                    if target_node_direction == start_direction:
                        next_edges = get_edges(next_node)
                    else:
                        if get_edges == self.get_forward_edges:
                            next_edges = self.get_backward_edges(next_node)
                        else:
                            next_edges = self.get_forward_edges(next_node)
                    for next_edge in next_edges:
                        final_node_hash = next_edge.get_targetNode().__hash__()
                        node_hashes = set(
                            e.get_targetNode().__hash__() for e in get_edges(start_node)
                        )
                        if final_node_hash in node_hashes:
                            self.process_triangle(
                                start_node, next_node, final_node_hash, triangles, seen_triangles
                            )
        return triangles

    def process_triangle(self, start_node, next_node, final_node_hash, triangles, seen_triangles):
        longest_path = [start_node.__hash__(), next_node.__hash__(), final_node_hash]
        longest_path_edges = [
            self.get_edges_between_nodes(
                self.get_node_by_hash(longest_path[i]), self.get_node_by_hash(longest_path[i + 1])
            )[0]
            for i in range(len(longest_path) - 1)
        ]
        shortest_path = [start_node.__hash__(), final_node_hash]
        shortest_path_edges = [
            self.get_edges_between_nodes(
                self.get_node_by_hash(shortest_path[i]), self.get_node_by_hash(shortest_path[i + 1])
            )[0]
            for i in range(len(shortest_path) - 1)
        ]
        # try:
        longest_path_coverage = [e.get_edge_coverage() for e in longest_path_edges]
        shortest_path_coverage = [e.get_edge_coverage() for e in shortest_path_edges]
        # except:
        #     [print(e.get_targetNode().__hash__()) for e in longest_path_edges]
        this_triangle = tuple(sorted(longest_path))
        if this_triangle not in seen_triangles:
            triangles.append(
                (
                    (longest_path, statistics.mean(longest_path_coverage)),
                    (shortest_path, shortest_path_coverage[0]),
                )
            )
            seen_triangles.add(this_triangle)

    def get_forward_node_from_node(self, sourceNode) -> list:
        # get the list of forward edge hashes for this node
        nodeForwardEdges = sourceNode.get_forward_edge_hashes()
        if len(nodeForwardEdges) == 1:
            # iterate through the edge hashes
            for edge_hash in nodeForwardEdges:
                # get the edge object corresponding to this edge hash
                edge = self.get_edges()[edge_hash]
                # get the target node for this edge
                targetNode = edge.get_targetNode()
                # get the degree of the target node
                targetNodeDegree = self.get_degree(targetNode)
                # get the direction we are going into the target node
                targetNodeDirection = edge.get_targetNodeDirection()
                # see if we can extend the linear path to the next node
                if (targetNodeDegree == 2 or targetNodeDegree == 1) and targetNode != sourceNode:
                    return True, targetNode, targetNodeDirection
                else:
                    return False, targetNode, targetNodeDirection
        return False, None, None

    def get_forward_path_from_node(
        self, node: Node, startDirection: int, wantBranchedNode=False
    ) -> list:
        """return a list of node hashes in the forward direction from this node"""
        forward_nodes_from_node = [node.__hash__()]
        if startDirection == 1:
            # get the next node in the forward direction
            forwardExtend, forwardNode, forwardNodeDirection = self.get_forward_node_from_node(node)
        else:
            # get the next node in the backward direction
            forwardExtend, forwardNode, forwardNodeDirection = self.get_backward_node_from_node(
                node
            )
        # if we are extending further in the forward direction, get the next canonical gene mer
        while forwardExtend:
            # if we enter the next node in the forward direction, we get the next forward node
            if forwardNodeDirection == 1:
                if not forward_nodes_from_node[0] == forwardNode.__hash__():
                    forward_nodes_from_node.append(forwardNode.__hash__())
                    forwardExtend, forwardNode, forwardNodeDirection = (
                        self.get_forward_node_from_node(forwardNode)
                    )
                else:
                    forwardExtend = False
            # if we enter the next node in the backward direction, we get the next backward node
            else:
                if not forward_nodes_from_node[0] == forwardNode.__hash__():
                    forward_nodes_from_node.append(forwardNode.__hash__())
                    forwardExtend, forwardNode, forwardNodeDirection = (
                        self.get_backward_node_from_node(forwardNode)
                    )
                else:
                    forwardExtend = False
        # if we want the final branched node too, then we add this to the path
        if wantBranchedNode and forwardNode:
            forward_nodes_from_node.append(forwardNode.__hash__())
        return forward_nodes_from_node

    def get_backward_node_from_node(self, sourceNode) -> list:
        # get the list of forward edge hashes for this node
        nodeBackwardEdges = sourceNode.get_backward_edge_hashes()
        if len(nodeBackwardEdges) > 0:
            # iterate through the edge hashes
            for edge_hash in nodeBackwardEdges:
                # get the edge object corresponding to this edge hash
                edge = self.get_edges()[edge_hash]
                # get the target node for this edge
                targetNode = edge.get_targetNode()
                # get the degree of the target node
                targetNodeDegree = self.get_degree(targetNode)
                # get the direction we are going into the target node
                targetNodeDirection = edge.get_targetNodeDirection()
                # see if we can extend the linear path to the next node
                if (targetNodeDegree == 2 or targetNodeDegree == 1) and targetNode != sourceNode:
                    # get the direction we are going into the target node
                    return True, targetNode, targetNodeDirection
                else:
                    # else we cannot extend the linear path
                    return False, targetNode, targetNodeDirection
        return False, None, None

    def get_backward_path_from_node(
        self, node: Node, startDirection, wantBranchedNode=False
    ) -> list:
        """return a list of node hashes in the backward direction from this node"""
        backward_nodes_from_node = [node.__hash__()]
        if startDirection == -1:
            # get the next node in the backward direction
            backwardExtend, backwardNode, backwardNodeDirection = self.get_backward_node_from_node(
                node
            )
        else:
            # get the next node in the forward direction
            backwardExtend, backwardNode, backwardNodeDirection = self.get_forward_node_from_node(
                node
            )
        # if we are extending further in the backward direction, get the next canonical gene mer
        count = 0
        while backwardExtend:
            if backwardNodeDirection == -1:
                if not backward_nodes_from_node[-1] == backwardNode.__hash__():
                    backward_nodes_from_node.insert(0, backwardNode.__hash__())
                    (
                        backwardExtend,
                        backwardNode,
                        backwardNodeDirection,
                    ) = self.get_backward_node_from_node(backwardNode)
                else:
                    backwardExtend = False
            # if we enter the next node in the forward direction, we get the next forward node
            else:
                if not backward_nodes_from_node[-1] == backwardNode.__hash__():
                    backward_nodes_from_node.insert(0, backwardNode.__hash__())
                    (
                        backwardExtend,
                        backwardNode,
                        backwardNodeDirection,
                    ) = self.get_forward_node_from_node(backwardNode)
                else:
                    backwardExtend = False
            count += 1
        # if we want the final branched node too, then we add this to the path
        if wantBranchedNode and backwardNode:
            backward_nodes_from_node.insert(0, backwardNode.__hash__())
        return backward_nodes_from_node

    def get_linear_path_for_node(self, node: Node, wantBranchedNode=False) -> list:
        """return a list of nodes in the linear path that contains the specified node. \
            Does not include the terminal nodes with a degree of more than 2"""
        backward_path = self.get_backward_path_from_node(
            node, -1 * node.get_geneMer().get_geneMerDirection(), wantBranchedNode
        )
        assert backward_path[-1] == node.__hash__()
        forward_path = self.get_forward_path_from_node(
            node, node.get_geneMer().get_geneMerDirection(), wantBranchedNode
        )
        assert forward_path[0] == node.__hash__()
        linear_path = backward_path[:-1] + [node.__hash__()] + forward_path[1:]
        return linear_path

    def get_all_node_coverages(self):
        """return an unordered list of all node coverages in the graph"""
        node_coverages = [node.get_node_coverage() for node in self.all_nodes()]
        return node_coverages

    def get_mean_node_coverage(self):
        """return an integer of the mean node coverage in the graph"""
        node_coverages = self.get_all_node_coverages()
        mean_node_coverage = sum(node_coverages) / len(node_coverages)
        return mean_node_coverage

    def generate_gml(
        self, output_file: str, geneMerSize: int, min_node_coverage: int, min_edge_coverage: int
    ):
        """Write a gml of the filtered graph to the output directory. \
            Returns the written content as a list"""
        graph_data = ["graph\t[", "multigraph 1"]
        self.assign_Id_to_nodes()
        # iterate through the nodes in the graph
        for sourceNode in self.all_nodes():
            # add the node entry
            nodeEntry = self.write_node_entry(
                sourceNode.get_node_Id(),
                self.get_gene_mer_label(sourceNode),
                sourceNode.get_node_coverage(),
                [read for read in sourceNode.get_reads()],
                sourceNode.get_component(),
                sourceNode.get_color(),
            )
            graph_data.append(nodeEntry)
            for edge in self.get_forward_edges(sourceNode) + self.get_backward_edges(sourceNode):
                if edge.get_edge_coverage() == 0:
                    continue
                targetNode = edge.get_targetNode()
                edgeEntry = self.write_edge_entry(
                    sourceNode.get_node_Id(),
                    targetNode.get_node_Id(),
                    edge.get_targetNodeDirection(),
                    edge.get_edge_coverage(),
                )
                graph_data.append(edgeEntry)
        graph_data.append("]")
        output_file = ".".join(
            [output_file, str(geneMerSize), str(min_node_coverage), str(min_edge_coverage)]
        )
        self.write_gml_to_file(output_file, graph_data)
        return graph_data

    def dfs_component(self, start, component_id, visited=None):
        if visited is None:
            visited = set()
        visited.add(start.__hash__())
        start.set_component(component_id)
        for neighbor in self.get_all_neighbors(start):
            if neighbor.__hash__() not in visited:
                self.dfs_component(neighbor, component_id, visited)

    def assign_component_ids(self):
        """use a depth first search algorithm to assign component IDs to nodes"""
        visited = set()
        component_id = 1
        for nodeHash in self.get_nodes():
            if nodeHash not in visited:
                self.dfs_component(self.get_node_by_hash(nodeHash), component_id, visited)
                component_id += 1

    def get_nodes_in_component(self, component):
        """return a list of nodes in the specified component"""
        nodes_in_component = []
        for nodeHash in self.get_nodes():
            node = self.get_node_by_hash(nodeHash)
            if node.get_component() == int(component):
                nodes_in_component.append(node)
        return nodes_in_component

    def components(self) -> list:
        """return a sorted list of all the possible component IDs"""
        components = set()
        for nodeHash in self.get_nodes():
            node = self.get_node_by_hash(nodeHash)
            components.add(node.get_component())
        return sorted(list(components))

    def get_number_of_component(self) -> int:
        """return an integer of the number of connected components in the graph"""
        return len(self.components())

    def get_highest_coverage_node_at_end_of_component(self, component_ID) -> Node:
        """return the node that has the highest coverage in the specified component"""
        # a place holder for the highest coverage node
        highest_coverage = 0
        # iterate through the nodes in this component
        for node in self.get_nodes_in_component(component_ID):
            if (
                len(self.get_forward_neighbors(node)) == 0
                or len(self.get_backward_neighbors(node)) == 0
            ):
                # get the coverage of the node
                coverage = node.get_node_coverage()
                # see if the coverage of this node is greater than the previous highest
                if coverage > highest_coverage:
                    highest_coverage = coverage
                    start_node = node
        return start_node

    def follow_path_by_coverage(self, startNode, traversal_direction, coverage):
        # ensure that the variables are in the correct form
        assert traversal_direction == 1 or traversal_direction == -1
        assert coverage == "highest" or coverage == "lowest"
        # use the coverage variable to decide how we are sorting the nodes
        if coverage == "highest":
            reverse_sort = True
        else:
            reverse_sort = False
        # start the chain at the start node
        path = [startNode.__hash__()]
        seen_nodes = {startNode.__hash__(): 1}
        current_node = startNode
        path_coverages = [startNode.get_node_coverage()]
        # use the direction to decide if we get the forward or reverse neighbors
        if traversal_direction == 1:
            next_neighbors = self.get_forward_neighbors(startNode)
        else:
            next_neighbors = self.get_backward_neighbors(startNode)
        # continue the chain until we reach the end of the component
        while len(next_neighbors) > 0:
            next_nodes = sorted(
                next_neighbors, key=lambda x: x.get_node_coverage(), reverse=reverse_sort
            )
            if not all(n.__hash__() in seen_nodes for n in next_nodes):
                for next_node in next_nodes:
                    next_node_hash = next_node.__hash__()
                    if next_node_hash not in seen_nodes:
                        path.append(next_node_hash)
                        if next_node_hash not in seen_nodes:
                            seen_nodes[next_node_hash] = 0
                        seen_nodes[next_node_hash] += 1
                        path_coverages.append(next_node.get_node_coverage())
                        sourceToTargetEdge, targetToSourceEdge = self.get_edges_between_nodes(
                            current_node, next_node
                        )
                        current_node = next_node
                        if sourceToTargetEdge.get_targetNodeDirection() == 1:
                            next_neighbors = self.get_forward_neighbors(current_node)
                        else:
                            next_neighbors = self.get_backward_neighbors(current_node)
                        break
            else:
                break
        return path, path_coverages

    def remove_low_coverage_components(self, min_component_coverage):
        # remove components that have a coverage less than min_component_coverage
        for component_ID in self.components():
            nodes_in_component = self.get_nodes_in_component(component_ID)
            if all(
                node.get_node_coverage() < min_component_coverage for node in nodes_in_component
            ):
                for node in nodes_in_component:
                    self.remove_node(node)

    def new_get_heaviest_path_through_component(self, component_ID) -> list:
        """return the heaviest path through a component"""
        # get the highest coverage node in this component
        start_nodes = []
        for node in self.get_nodes_in_component(component_ID):
            if len(self.get_forward_neighbors(node)) == 0:
                start_nodes.append((node, -1))
            elif len(self.get_backward_neighbors(node)) == 0:
                start_nodes.append((node, 1))
            else:
                pass
        # start_node = self.get_highest_coverage_node_at_end_of_component(component_ID)
        paths = []
        for start_node in start_nodes:
            # get the heaviest path in the forward direction
            heaviest_path, heaviest_coverages = self.follow_path_by_coverage(
                start_node[0], start_node[1], "highest"
            )
            paths.append(heaviest_path)
        sorted_paths = sorted(paths, key=len, reverse=True)
        heaviest_path = sorted_paths[0]
        # combine the foward and backward paths to get the heaviest path through the component
        return heaviest_path

    def reverse_list_of_genes(self, list_of_genes: list[str]) -> list[str]:
        return [("-" if g[0] == "+" else "+") + g[1:] for g in reversed(list_of_genes)]

    def make_distance_matrix(self, correct, incorrect):
        # get the forward and reverse list of genes for each correct node
        from collections import Counter

        correct_genes = []
        for node in correct:
            fw = self.get_gene_mer_label(node).split("~~~")
            rv = self.reverse_list_of_genes(fw)
            correct_genes.append((fw, rv))
        replacement = {}
        for node in tqdm(incorrect):
            list_of_genes = self.get_gene_mer_label(node).split("~~~")
            counter1 = Counter(list_of_genes)
            distances = []
            for i in range(len(correct_genes)):
                closest = sorted(
                    [correct_genes[i][0][:], correct_genes[i][1][:]],
                    key=lambda x: len(set(list_of_genes).intersection(set(x))),
                    reverse=True,
                )[0]
                counter2 = Counter(closest)
                shared_count = len(list((counter1 & counter2).elements()))
                if not shared_count == 0:
                    union_count = len(list((counter1 | counter2).elements()))
                    distances.append(shared_count / len(list_of_genes))
            if not distances == []:
                max_index, max_value = max(enumerate(distances), key=lambda x: x[1])
                assert not distances.count(max_value) > 1, distances
                print(max_value)

    def replace_low_coverage_nodes(self, node_min_coverage):
        """separates nodes into correct and incorrect nodes based on coverage and corrects all incorrect nodes to the closest correct node"""
        self.set_minNodeCoverage(node_min_coverage)
        correct = []
        incorrect = []
        for node_hash in self.get_nodes():
            node = self.get_node_by_hash(node_hash)
            if node.get_node_coverage() >= self.get_minNodeCoverage():
                correct.append(node)
            else:
                incorrect.append(node)
        self.make_distance_matrix(correct, incorrect)
        return self

    def get_high_coverage_junction_nodes_in_component(self, component):
        """return a list of nodes with >2 high coverage neighbors"""
        # get all the nodes in the component
        nodes_in_component = self.get_nodes_in_component(component)
        # get a set of the high coverage nodes in this component
        high_coverage_nodes = set()
        for node in nodes_in_component:
            if node.get_node_coverage() >= 10:
                high_coverage_nodes.add(node.__hash__())
        # get a list of the high coverage nodes with more than 2 high coverage neighbors
        high_coverage_junctions = []
        for node_hash in high_coverage_nodes:
            node = self.get_node_by_hash(node_hash)
            neighbor_hashes = self.get_all_neighbor_hashes(node)
            high_coverage_count = 0
            for nh in neighbor_hashes:
                if nh in high_coverage_nodes:
                    high_coverage_count += 1
            if high_coverage_count > 2:
                high_coverage_junctions.append(node)
        return high_coverage_junctions

    def check_if_component_is_resolved(self, component):
        """returns a bool of whether a component is resolved or not"""
        # get all the nodes in the component
        high_coverage_junctions = self.get_high_coverage_junction_nodes_in_component(component)
        if not len(high_coverage_junctions) == 0:
            return False, high_coverage_junctions
        else:
            return True, high_coverage_junctions

    def traverse_graph_from_high_coverage_junctions(self, high_coverage_junctions):
        high_coverage_junction_hashes = set([n.__hash__() for n in high_coverage_junctions])
        segments_to_correct_to = set()
        fw_segments = set()
        bw_segments = set()
        for node in high_coverage_junctions:
            # segments_to_correct_to[node.__hash__()] = set()
            # get the high coverage neighbors
            fw_nodes = [n for n in self.get_forward_neighbors(node) if n.get_node_coverage() >= 10]
            bw_nodes = [n for n in self.get_backward_neighbors(node) if n.get_node_coverage() >= 10]
            # get the heaviest path from each neighbor to the end of the component or another high coverage neighbor
            for neighbor in tqdm(fw_nodes):
                # get the start direction
                sourceToTargetEdge, targetToSourceEdge = self.get_edges_between_nodes(
                    node, neighbor
                )
                start_direction = sourceToTargetEdge.get_targetNodeDirection()
                # get the heaviest path in the start direction
                heaviest_path, heaviest_coverages = self.follow_path_by_coverage(
                    neighbor, start_direction, "highest"
                )
                # slice the path if there are any other junction nodes
                for index, item in enumerate(heaviest_path[1:]):
                    if item in high_coverage_junction_hashes:
                        heaviest_path = heaviest_path[: index + 2]
                        assert heaviest_path[-1] in high_coverage_junction_hashes
                        break
                heaviest_path = [node.__hash__()] + heaviest_path
                assert heaviest_path[0] == node.__hash__()
                assert heaviest_path[1] == neighbor.__hash__()
                segments_to_correct_to.add(
                    tuple(sorted([heaviest_path, list(reversed(heaviest_path))])[0])
                )
                fw_segments.add(tuple(heaviest_path))
                # segments_to_correct_to[node.__hash__()].add(tuple(sorted([heaviest_path, list(reversed(heaviest_path))])[0]))
            for neighbor in tqdm(bw_nodes):
                # get the start direction
                sourceToTargetEdge, targetToSourceEdge = self.get_edges_between_nodes(
                    node, neighbor
                )
                start_direction = sourceToTargetEdge.get_targetNodeDirection()
                # get the heaviest path in the start direction
                heaviest_path, heaviest_coverages = self.follow_path_by_coverage(
                    neighbor, start_direction, "highest"
                )
                # slice the path if there are any other junction nodes
                for index, item in enumerate(heaviest_path[1:]):
                    if item in high_coverage_junction_hashes:
                        heaviest_path = heaviest_path[: index + 2]
                        assert heaviest_path[-1] in high_coverage_junction_hashes
                        break
                heaviest_path = [node.__hash__()] + heaviest_path
                assert heaviest_path[0] == node.__hash__()
                assert heaviest_path[1] == neighbor.__hash__()
                segments_to_correct_to.add(
                    tuple(sorted([heaviest_path, list(reversed(heaviest_path))])[0])
                )
                bw_segments.add(tuple(list(reversed(heaviest_path))))
        #                segments_to_correct_to[node.__hash__()].add(tuple(sorted([heaviest_path, list(reversed(heaviest_path))])[0]))
        return segments_to_correct_to, list(fw_segments), list(bw_segments)

    def enumerate_paths(self, graph):

        def find_paths_local(node, path):
            path.append(node)
            # If the node has no adjacent nodes or all paths from it are explored, add the current path to the set of paths
            if node not in graph or all(neighbor in path for neighbor in graph[node]):
                node_path = list(path[0]) + [p[-1] for p in path[1:]]
                all_paths.add(tuple(node_path))
                return
            # Explore adjacent nodes
            for neighbor in graph[node]:
                if neighbor not in path:  # Avoid loops
                    find_paths_local(neighbor, path.copy())

        all_paths = set()
        for start_node in graph:
            find_paths_local(start_node, [])
        return all_paths

    def join_segments(self, fw_segments, bw_segments):
        new_segments = {}
        for b in bw_segments:
            new_segments[b] = set()
            for f in fw_segments:
                if b[-1] == f[0]:
                    new_segments[b].add(f)
        return

    def get_all_reads_not_containing_nodes(self, setOfNodeHashes: set, component: int) -> set:
        reads_to_correct = set()
        for node in self.get_nodes_in_component(component):
            node_hash = node.__hash__()
            if not node_hash in setOfNodeHashes:
                for readId in node.get_reads():
                    reads_to_correct.add(readId)
        return reads_to_correct

    def get_AMR_nodes(self, listOfGenes):
        AMR_nodes = {}
        # iterate through the list of specified genes
        for geneOfInterest in tqdm(listOfGenes):
            # get the graph nodes containing this gene
            nodesOfInterest = self.get_nodes_containing(geneOfInterest)
            # add nodes of interest to the AMR node set
            for node in nodesOfInterest:
                AMR_nodes[node.__hash__()] = node
        return AMR_nodes

    def get_AMR_anchors_and_junctions(self, AMRNodes):
        # initialise the node anchor and junction sets
        nodeAnchors = set()
        nodeJunctions = set()
        # get nodes that are anchors for traversing in the forward direction of the graph
        for nodeHash in AMRNodes:
            node = AMRNodes[nodeHash]
            if not self.get_degree(node) > 2:
                # get the forward neighbors for this node
                forward_neighbors = self.get_forward_neighbors(node)
                # get the backward neighbors for this node
                backward_neighbors = self.get_backward_neighbors(node)
                # get the forward neighbors that contain an AMR gene
                forwardAMRNodes = [n for n in forward_neighbors if n.__hash__() in AMRNodes]
                # get the backward neighbors that contain an AMR gene
                backwardAMRNodes = [n for n in backward_neighbors if n.__hash__() in AMRNodes]
                if len(backwardAMRNodes) == 0 or len(forwardAMRNodes) == 0:
                    nodeAnchors.add(nodeHash)
            # if the number of backward neighbors is not 1 or 0 then this is a forward anchor
            else:
                nodeJunctions.add(nodeHash)
        return nodeAnchors, nodeJunctions

    def create_adjacency_matrix(self, nodeHashesOfInterest):
        size = len(nodeHashesOfInterest)
        matrix = np.zeros((size, size), dtype=int)
        node_index = {n: idx for idx, n in enumerate(nodeHashesOfInterest)}
        for node_hash in nodeHashesOfInterest:
            node = self.get_node_by_hash(node_hash)
            for neighbor_hash in self.get_all_neighbor_hashes(node):
                if neighbor_hash in node_index:
                    matrix[node_index[node_hash], node_index[neighbor_hash]] = 1
        return matrix

    def find_paths(self, matrix, start, end, path=[]):
        path = path + [start]
        if start == end:
            return [path]
        paths = []
        for neighbor, connected in enumerate(matrix[start]):
            if connected and neighbor not in path:
                newpaths = self.find_paths(matrix, neighbor, end, path)
                for newpath in newpaths:
                    paths.append(newpath)
        return paths

    def all_paths_for_subgraph(self, nodeHashesOfInterest, anchor_nodes):
        matrix = self.create_adjacency_matrix(nodeHashesOfInterest)
        paths = {}
        for i in range(len(nodeHashesOfInterest)):
            for j in range(len(nodeHashesOfInterest)):
                sorted_indices = tuple(sorted([i, j]))
                sorted_pairs = (
                    nodeHashesOfInterest[sorted_indices[0]],
                    nodeHashesOfInterest[sorted_indices[1]],
                )
                if (
                    i != j
                    and sorted_pairs not in paths
                    and (
                        nodeHashesOfInterest[i] in anchor_nodes
                        and nodeHashesOfInterest[j] in anchor_nodes
                    )
                ):
                    returned_paths = [
                        [nodeHashesOfInterest[x] for x in p]
                        for p in self.find_paths(
                            matrix, sorted_indices[0], sorted_indices[1])
                    ]
                    if not returned_paths == []:
                        paths[sorted_pairs] = returned_paths
        return paths

    def get_anchors_of_interest(self, nodeHashesOfInterest):
        # initialise the node anchor and junction sets
        nodeAnchors = set()
        nodeJunctions = set()
        # get nodes that are anchors for traversing in the forward direction of the graph
        for node_hash in nodeHashesOfInterest:
            node = self.get_node_by_hash(node_hash)
            forwardAMRNodes = [
                    n
                    for n in self.get_forward_neighbors(node)
                    if n.__hash__() in nodeHashesOfInterest
                ]
            backwardAMRNodes = [
                    n
                    for n in self.get_backward_neighbors(node)
                    if n.__hash__() in nodeHashesOfInterest
                ]
            if len(backwardAMRNodes) == 0 or len(forwardAMRNodes) == 0:
                nodeAnchors.add(node_hash)
            if len(backwardAMRNodes) > 1 or len(forwardAMRNodes) > 1:
                nodeJunctions.add(node_hash)
        return nodeAnchors, nodeJunctions

    def is_sublist(self, long_list, sub_list):
        """Check if list is a sublist of long_list."""
        len_sub = len(sub_list)
        return any(
            sub_list == long_list[i : i + len_sub] for i in range(len(long_list) - len_sub + 1)
        )

    def assign_reads_to_genes(self, listOfGenes):
        # initialise a dictionary to store the paths
        clusters_of_interest = {}
        # iterate through the genes we are interested in
        for geneOfInterest in tqdm(listOfGenes):
            # get the graph nodes containing this gene
            nodeOfInterest = self.get_nodes_containing(geneOfInterest)
            # get the node hashes containing this gene
            nodeHashesOfInterest = [n.__hash__() for n in nodeOfInterest]
            # get the nodes that contain this gene but are at the end of a path
            anchor_nodes, junction_nodes = self.get_anchors_of_interest(nodeHashesOfInterest)
            # get all of the possible paths for these nodes
            paths = self.all_paths_for_subgraph(nodeHashesOfInterest, anchor_nodes)
            # get the reads covering the node hashes of interest
            reads_of_interest = set()
            for node in nodeOfInterest:
                for read_id in node.get_reads():
                    reads_of_interest.add(read_id)
            # check if any of the paths are longer than k, this means there are duplicates
            for pair in paths:
                for p in range(len(paths[pair])):
                    path = paths[pair][p]
                    if len(path) > self.get_kmerSize():
                        # get the list of genes for this unitig
                        genes = self.get_genes_in_unitig(path)
                        nodes = list(path)
                        node_indices_for_each_gene_index = {}
                        for i in range(len(nodes)):
                            gene_indices = [j for j in range(i, i + self.get_kmerSize())]
                            for g in gene_indices:
                                if genes[g][1:] == geneOfInterest:
                                    if g not in node_indices_for_each_gene_index:
                                        node_indices_for_each_gene_index[g] = []
                                    node_indices_for_each_gene_index[g].append(i)
                        new_paths = []
                        for i in node_indices_for_each_gene_index:
                            nodes_containing_this_gene = [
                                nodes[j] for j in node_indices_for_each_gene_index[i]
                            ]
                            new_paths.append(nodes_containing_this_gene)
                        # remove the old path
                        del paths[pair][p]
                        # add the new paths
                        paths[pair] += new_paths
            # assign reads to paths if the path or reversed path is in the list of nodes
            clustered_reads = {}
            for read_id in reads_of_interest:
                # get the nodes on the read
                nodes_on_read = self.get_readNodes()[read_id]
                # reversed the nodes on the read
                reversed_node_on_read = list(reversed(nodes_on_read))
                # iterate through the paths
                for pair in paths:
                    for p in paths[pair]:
                        if any(n in junction_nodes for n in p):
                            if self.is_sublist(nodes_on_read, p) or self.is_sublist(
                                reversed_node_on_read, p
                            ):
                                path_tuple = tuple(p)
                                if not path_tuple in clustered_reads:
                                    clustered_reads[path_tuple] = set()
                                clustered_reads[path_tuple].add(read_id)
                        else:
                            if len(set(p).intersection(set(nodes_on_read))) > 0:
                                path_tuple = tuple(p)
                                if not path_tuple in clustered_reads:
                                    clustered_reads[path_tuple] = set()
                                clustered_reads[path_tuple].add(read_id)
            # store the paths
            copy_count = 0
            for path in clustered_reads:
                clusters_of_interest[
                    f"{geneOfInterest}_{copy_count}_{len(clustered_reads[path])}"
                ] = list(clustered_reads[path])
                copy_count += 1
        return clusters_of_interest

    def extract_elements(self, lst):
        result = []
        for i in range(len(lst)):
            if lst[i] != 0:
                result.append(lst[i])
            elif i < len(lst) - 1 and lst[i + 1] != 0:
                result.append(lst[i])
        return result

    def find_paths_between_nodes(
        self, start_hash, end_hash, distance, current_direction, path=None, seen_nodes=None
    ):
        if path is None:
            path = []
        if seen_nodes is None:
            seen_nodes = set()

        path = path + [(start_hash, current_direction)]
        seen_nodes.add(start_hash)  # Add the current node to the seen nodes to avoid revisiting

        if end_hash:
            if start_hash == end_hash and len(path) - 1 <= distance:
                return [path]
        else:
            if len(path) == distance:
                return [path]

        if len(path) > distance:
            return []

        paths = []
        next_nodes = []
        if current_direction == 1:
            next_nodes = self.get_node_by_hash(start_hash).get_forward_edge_hashes()
        elif current_direction == -1:
            next_nodes = self.get_node_by_hash(start_hash).get_backward_edge_hashes()

        if not all(
            self.get_edge_by_hash(edge_hash).get_targetNode().__hash__() in seen_nodes
            for edge_hash in next_nodes
        ):
            for edge_hash in next_nodes:
                edge = self.get_edge_by_hash(edge_hash)
                node_hash = edge.get_targetNode().__hash__()
                if node_hash not in seen_nodes:
                    newpaths = self.find_paths_between_nodes(
                        node_hash,
                        end_hash,
                        distance,
                        edge.get_targetNodeDirection(),
                        path,
                        seen_nodes.copy(),  # Pass a copy of seen_nodes
                    )
                    for newpath in newpaths:
                        paths.append(newpath)

        return paths

    def insert_valid_paths(self, replacements, node_list, node_directions):
        offset = 0  # Keep track of the change in length of the list
        modified_node_directions = node_directions[:]
        assert len(node_list) == len(modified_node_directions)
        for (start, end), values in sorted(replacements.items(), key=lambda x: x[0][0]):
            # Adjust indices based on the current offset
            adjusted_start = start + offset
            adjusted_end = end + offset + 1  # +1 to make 'end' inclusive
            # Calculate the difference in length between the new values and the original range
            length_difference = len(values) - (adjusted_end - adjusted_start)
            # Replace the elements from 'adjusted_start' to 'adjusted_end' with 'values'
            node_list[adjusted_start:adjusted_end] = values
            modified_node_directions[adjusted_start + 1 : adjusted_end - 1] = [None] * (
                len(values) - 2
            )
            # Update the offset for the next replacement
            offset += length_difference
        assert len(node_list) == len(modified_node_directions)
        return node_list, modified_node_directions, offset

    def correct_reads(self):
        """Main function to correct reads in the dataset."""
        readNodes = self.get_readNodes()  # Retrieve read nodes
        corrected_genes = {}
        for read_id in tqdm(readNodes):
            list_of_genes = self.correct_single_read(read_id, readNodes)
            if len(list_of_genes) > 0:
                corrected_genes[read_id] = list_of_genes
        return corrected_genes

    def correct_single_read(self, read_id, readNodes):
        """Correct a single read based on its ID and the readNodes data."""
        if read_id not in self.get_reads_to_correct():
            return self.get_reads()[read_id]

        if not all(n is None for n in readNodes[read_id]):
            start, end = self.find_read_boundaries(readNodes[read_id])
            return self.process_read_correction(read_id, readNodes, start, end)
        return []
        # return self.get_reads()[read_id]

    def find_read_boundaries(self, readNode):
        """Find the first and last non-None positions in a read."""
        start, end = 0, len(readNode) - 1
        for i, node in enumerate(readNode):
            if node:
                start = i
                break
        for i, node in enumerate(reversed(readNode)):
            if node:
                end = len(readNode) - 1 - i
                break
        return start, end

    def insert_elements(self, base_list, insert_dict):
        from itertools import product

        if len(insert_dict) == 0:
            return [base_list]

        def get_all_combinations(dict_of_lists):
            # Transform each list into a list of (key, element) tuples
            lists_with_keys = [
                [(key, element) for element in lst] for key, lst in dict_of_lists.items()
            ]

            # Generate all combinations using the transformed lists
            return list(product(*lists_with_keys))

        # Create a list of all combinations of insertions
        all_combinations = get_all_combinations(insert_dict)

        # Function to insert a single combination of elements
        def insert_single_combination(base_list, combination):
            inserted_list = base_list[:]
            offset = 0
            for (start, end), path in list(combination):
                insertion_point = start + offset
                # Remove the elements that will be replaced by the new path
                del inserted_list[insertion_point : end + offset + 1]
                # Insert the new path
                inserted_list[insertion_point:insertion_point] = path
                insertion_point += len(path)
                offset += len(path) - (end - start + 1)
            return inserted_list

        # Generate all possible lists with inserted elements
        result_lists = []
        for combination in all_combinations:
            result_lists.append(insert_single_combination(base_list, combination))

        return result_lists

    def get_possible_paths(self, nodes_on_read, replacementDict, start, end):
        # get all the possible middle paths we could correct to
        possible_middle_paths = self.insert_elements(nodes_on_read, replacementDict)
        # Initialize lists for upstream and downstream paths
        upstream_paths, downstream_paths = [], []
        # get all the possible upstream paths
        if start != 0:
            upstream_paths = self.find_paths_between_nodes(
                nodes_on_read[start][0], None, start * 2, -nodes_on_read[start][1]
            )
        # get all the possible downstream paths
        if end != len(nodes_on_read) - 1:
            downstream_paths = self.find_paths_between_nodes(
                nodes_on_read[end][0],
                None,
                (len(nodes_on_read) - end - 1) * 2,
                nodes_on_read[end][1],
            )
        possible_paths = []
        # for now we are not going to add the upstream and downstream paths
        #upstream_paths, downstream_paths = [], []
        # Iterate over each possible middle path
        for corrected in possible_middle_paths:
            # Extract nodes and directions from the middle path, ignoring entries without a node
            corrected_path = [n[0] for n in corrected if n[0]]
            corrected_directions = [n[1] for n in corrected if n[0]]

            # Ensure there's at least an empty list to iterate over if upstream_paths is empty
            upstream_paths = upstream_paths or [[]]
            # Ensure there's at least an empty list to iterate over if downstream_paths is empty
            downstream_paths = downstream_paths or [[]]
            # Iterate over each upstream path
            for u_path in upstream_paths:
                # Reverse and negate directions for upstream path, ignore if u_path is empty
                u_nodes = [n[0] for n in reversed(u_path)] if u_path else []
                u_directions = [-n[1] for n in reversed(u_path)] if u_path else []

                # Iterate over each downstream path
                for d_path in downstream_paths:
                    # Extract nodes and directions from downstream path, ignore if d_path is empty
                    d_nodes = [n[0] for n in d_path] if d_path else []
                    d_directions = [n[1] for n in d_path] if d_path else []

                    # Combine the paths, excluding the last node of the upstream and the first node of the downstream
                    combined_path = (
                        (u_nodes[:-1] if u_nodes else [])
                        + corrected_path
                        + (d_nodes[1:] if d_nodes else [])
                    )
                    combined_directions = (
                        (u_directions[:-1] if u_directions else [])
                        + corrected_directions
                        + (d_directions[1:] if d_directions else [])
                    )

                    # Add the combined path and its directions to the list of possible paths
                    possible_paths.append((combined_path, combined_directions))
        return possible_paths

    def get_coverage_of_path(self, path):
        coverages = [self.get_node_by_hash(n).get_node_coverage() for n in path]
        return statistics.mean(coverages)

    def process_read_correction(self, read_id, readNodes, start, end):
        nodes_on_read = [
            (readNodes[read_id][i], self.get_readNodeDirections()[read_id][i])
            for i in range(len(readNodes[read_id]))
        ]
        # correct filtered nodes in the middle of a read
        path_terminals = self.identify_path_terminals(readNodes[read_id], start, end)
        # chop the list of genes if the ends have been cut off
        if len(path_terminals) == 0:
            new_nodes = [n[0] for n in nodes_on_read[start: end + 1]]
            new_directions = [n[1] for n in nodes_on_read[start: end + 1]]
            return self.get_annotation_for_read(new_nodes, new_directions, read_id)
        # find alternative paths where nodes have been deleted
        replacementDict = {}
        for pair in path_terminals:
            replacementDict.update(self.generate_replacement_dict(nodes_on_read, pair))
        # recursively get all possible corrected path options
        possible_paths = self.get_possible_paths(nodes_on_read, replacementDict, start, end)
        if possible_paths == []:
            return self.get_reads()[read_id]
        # get the genes for each potential path
        distance = 0
        coverage = 0
        for path in possible_paths:
            path_mean_coverage = self.get_coverage_of_path(path[0])
            genes = self.get_annotation_for_read(path[0], path[1], read_id)
            this_distance = len(set(genes).intersection(self.get_reads()[read_id]))
            if this_distance > distance:
                closest = genes
                distance = this_distance
                coverage = path_mean_coverage
            elif this_distance == distance and path_mean_coverage > coverage:
                closest = genes
                distance = this_distance
                coverage = path_mean_coverage
            else:
                pass
        return closest

    def get_annotation_for_read(self, listOfNodes, listOfNodeDirections, read_id):
        # check that the length of the nodes and node directions is equal
        assert len(listOfNodes) == len(
            listOfNodeDirections
        ), f"The number of nodes and node directions for read {read_id} are not the same"
        # get the read node directions
        if not listOfNodeDirections:
            listOfNodeDirections = self.get_readNodeDirections()[read_id]
        # return the gene-mer or rc gene-mer if there is one node
        if len(listOfNodes) == 1:
            geneMer_direction = listOfNodeDirections[0]
            if geneMer_direction == 1:
                return self.get_gene_mer_genes(self.get_node_by_hash(listOfNodes[0]))
            elif geneMer_direction == -1:
                return self.get_reverse_gene_mer_genes(self.get_node_by_hash(listOfNodes[0]))
            else:
                raise ValueError(
                    f"Gene-mer direction for a node with 1 read cannot be {geneMer_direction}"
                )
        # initialise a list to get the new annotations
        new_annotations = []
        for n in range(len(listOfNodes)):
            # get the node we are adding
            node = self.get_node_by_hash(listOfNodes[n])
            # get the node direction
            node_direction = listOfNodeDirections[n]
            # add all the gene-mer genes if this is the first node
            if n == 0:
                if node_direction == 1:
                    new_annotations += self.get_gene_mer_genes(node)[:-1]
                else:
                    new_annotations += self.get_reverse_gene_mer_genes(node)[:-1]
            # if the node direction is not None
            if node_direction:
                if node_direction:
                    if node_direction == 1:
                        new_annotations.append(self.get_gene_mer_genes(node)[-1])
                    else:
                        new_annotations.append(self.get_reverse_gene_mer_genes(node)[-1])
                else:
                    new_annotations.append(None)
        assert None not in new_annotations
        return new_annotations

    def identify_path_terminals(self, corrected, start, end):
        """Identify terminals for paths that need correction within a read."""
        path_terminals = []
        for i in range(len(corrected)):
            if i >= start and i <= end:
                if not corrected[i]:
                    if corrected[i - 1]:
                        path_start = i - 1
                    if corrected[i + 1]:
                        path_end = i + 1
                        path_terminals.append((path_start, path_end))
        return path_terminals

    def generate_replacement_dict(self, corrected, pair):
        """Generate a dictionary with replacements for invalid nodes."""
        paths = self.find_paths_between_nodes(
            corrected[pair[0]][0],
            corrected[pair[1]][0],
            self.get_kmerSize() * 2,
            corrected[pair[0]][1],
        )
        return {pair: paths}

    def remove_junk_reads(self, error_rate):
        new_reads = {}
        for read_id in self.get_readNodes():
            number_of_nodes = len(self.get_readNodes()[read_id])
            expected_filtered_nodes = round(number_of_nodes * (1 - error_rate))
            number_of_filtered_nodes = self.get_readNodes()[read_id].count(None)
            if not number_of_filtered_nodes > expected_filtered_nodes:
                new_reads[read_id] = self.get_reads()[read_id]
        return new_reads

    def get_valid_reads_only(self):
        valid_reads = {}
        for read_id in self.get_reads():
            if not read_id in self.get_reads_to_correct():
                valid_reads[read_id] = self.get_reads()[read_id]
        return valid_reads

    def old_pop_bubbles(self):
        # get nodes with >1 neighbour in the forward or backward direction
        potential_bubble_starts = []
        all_bubble_nodes = set()
        for node in self.all_nodes():
            # get the forward edges
            fw_edges = node.get_forward_edge_hashes()
            # get the backward edges
            bw_edges = node.get_backward_edge_hashes()
            if len(fw_edges) > 1:
                potential_bubble_starts.append((node, 1))
                all_bubble_nodes.add(node.__hash__())
            if len(bw_edges) > 1:
                potential_bubble_starts.append((node, -1))
                all_bubble_nodes.add(node.__hash__())
        # define the maximum distance
        max_distance = self.get_kmerSize() + 2
        # iterate through the potential bubble starts
        seen_bubbles = set()
        for start in potential_bubble_starts:
            node = start[0]
            direction = start[1]
            if direction == 1:
                next_nodes = [self.get_edge_by_hash(edge_hash).get_targetNode() for edge_hash in node.get_forward_edge_hashes()]
            if direction == -1:
                next_nodes = [self.get_edge_by_hash(edge_hash).get_targetNode() for edge_hash in node.get_backward_edge_hashes()]
            potential_paths = {}
            for neighbor in next_nodes:
                path = self.get_linear_path_for_node(neighbor, True)
                terminal_nodes = tuple(sorted([path[0], path[-1]]))
                if not terminal_nodes in seen_bubbles:
                    if 2 < len(path) <= max_distance:
                        if not terminal_nodes in potential_paths:
                            potential_paths[terminal_nodes] = []
                        path_coverage = statistics.mean([self.get_node_by_hash(n).get_node_coverage() for n in path[1:-1]])
                        potential_paths[terminal_nodes].append((path, path_coverage))
            for pair in potential_paths:
                if len(potential_paths[pair]) > 1:
                    sorted_by_length = sorted(potential_paths[pair], key=lambda x: x[1], reverse=True)
                    if sorted_by_length[-1][1] < 10:
                        for node_hash in sorted_by_length[-1][0][1:-1]:
                            node = self.get_node_by_hash(node_hash)
                            self.remove_node(node)
                    seen_bubbles.add(pair)

    def pop_bubbles(self):
        # get nodes with >1 neighbour in the forward or backward direction
        potential_bubble_starts = []
        all_bubble_nodes = set()
        for node in self.all_nodes():
            # get the forward edges
            fw_edges = node.get_forward_edge_hashes()
            # get the backward edges
            bw_edges = node.get_backward_edge_hashes()
            if len(fw_edges) > 1:
                potential_bubble_starts.append((node, 1))
                all_bubble_nodes.add(node.__hash__())
            if len(bw_edges) > 1:
                potential_bubble_starts.append((node, -1))
                all_bubble_nodes.add(node.__hash__())
        # define the maximum distance
        max_distance = self.get_kmerSize() + 2
        # iterate through the potential bubble starts
        seen_bubbles = set()
        for start in potential_bubble_starts:
            start_node = start[0]
            start_node_hash = start_node.__hash__()
            direction = start[1]
            potential_paths = {}
            for end_hash in all_bubble_nodes:
                bubble_terminals = tuple(sorted([start_node_hash, end_hash]))
                if not bubble_terminals in seen_bubbles:
                    paths = self.find_paths_between_nodes(start_node_hash, end_hash, self.get_kmerSize() + 2, direction)
                    for p in paths:
                        if 2 < len(p) <= max_distance:
                            if not bubble_terminals in potential_paths:
                                potential_paths[bubble_terminals] = set()
                            path_coverage = statistics.mean([self.get_node_by_hash(n[0]).get_node_coverage() for n in p[1:-1]])
                            potential_paths[bubble_terminals].add((tuple([n[0] for n in p]), path_coverage))
            for pair in potential_paths:
                if len(potential_paths[pair]) > 1:
                    sorted_by_length = sorted(list(potential_paths[pair]), key=lambda x: x[1], reverse=True)
                    if sorted_by_length[-1][1] < 10:
                        for node_hash in list(sorted_by_length[-1][0])[1:-1]:
                            if not node_hash in all_bubble_nodes:
                                try:
                                    node = self.get_node_by_hash(node_hash)
                                    self.remove_node(node)
                                except:
                                    pass
                    seen_bubbles.add(pair)

    def needleman_wunsch(self, x, y):
        N, M = len(x), len(y)
        # Scoring function: returns 1 if elements are equal, 0 otherwise
        s = lambda a, b: int(a == b)
        # Direction constants for traceback
        DIAG, LEFT, UP = (-1, -1), (-1, 0), (0, -1)
        # Initialize score (F) and pointer (Ptr) matrices
        F, Ptr = {}, {}
        F[-1, -1] = 0
        ##### I HAVE GOT RID OF GAP PENALTIES AT THE START OF THE ALIGNMENT #####
        # Initial scoring for gaps along x
        for i in range(N):
            F[i, -1] = -i
        # Initial scoring for gaps along y
        for j in range(M):
            F[-1, j] = -j
        #########################################################################
        # Option for Ptr to trace back alignment
        option_Ptr = DIAG, LEFT, UP
        # Fill F and Ptr tables
        for i, j in product(range(N), range(M)):
            # Score options: match/mismatch, gap in x, gap in y
            option_F = (
                F[i - 1, j - 1] + s(x[i], y[j]),  # Match/mismatch
                F[i - 1, j] - 1,  # Gap in x
                F[i, j - 1] - 1,  # Gap in y
            )
            # Choose best option for F and Ptr
            F[i, j], Ptr[i, j] = max(zip(option_F, option_Ptr))
        # Trace back to get the alignment
        alignment = deque()
        i, j = N - 1, M - 1
        while i >= 0 and j >= 0:
            direction = Ptr[i, j]
            # Add aligned elements or gaps based on direction
            if direction == DIAG:
                element = x[i], y[j]
            elif direction == LEFT:
                element = x[i], "*"  # Insert gap in y
            elif direction == UP:
                element = "*", y[j]  # Insert gap in x
            alignment.appendleft(element)
            di, dj = direction
            i, j = i + di, j + dj
        # Add remaining gaps if any
        while i >= 0:
            alignment.appendleft((x[i], "*"))  # Gap in y
            i -= 1
        while j >= 0:
            alignment.appendleft(("*", y[j]))  # Gap in x
            j -= 1
        return list(alignment)

#     def identify_potential_bubble_starts(self, node_min_coverage):
#         potential_bubble_starts = []
#         all_bubble_nodes = set()
#         for node in self.all_nodes():
#             #if node.get_node_coverage() >= node_min_coverage:
# #            if node.get_component() in {2, 3}:
#             fw_edges = node.get_forward_edge_hashes()  # Forward edges
#             bw_edges = node.get_backward_edge_hashes()  # Backward edges
#             if len(fw_edges) > 1 and len(bw_edges) == 1:
#                 potential_bubble_starts.append((node, 1))
#                 all_bubble_nodes.add(node.__hash__())
#             if len(bw_edges) > 1 and len(fw_edges) == 1:
#                 potential_bubble_starts.append((node, -1))
#                 all_bubble_nodes.add(node.__hash__())
#         return potential_bubble_starts, all_bubble_nodes

    def old_find_potential_paths(self, start, all_bubble_nodes, seen_bubbles, max_distance):
        potential_paths = {}
        start_node, direction = start
        start_node_hash = start_node.__hash__()
        for end_hash in all_bubble_nodes:
            bubble_terminals = tuple(sorted([start_node_hash, end_hash]))
            if bubble_terminals not in seen_bubbles:
                paths = self.find_paths_between_nodes(start_node_hash, end_hash, max_distance, direction)
                for p in paths:
                    if 2 < len(p) <= max_distance:
                        if bubble_terminals not in potential_paths:
                            potential_paths[bubble_terminals] = set()
                        path_coverage = self.calculate_path_coverage(p)
                        potential_paths[bubble_terminals].add((tuple([n[0] for n in p]), path_coverage))
        return potential_paths

    def calculate_path_coverage(self, path):
        #return statistics.mean([self.get_node_by_hash(n[0]).get_node_coverage() for n in path[1:-1]])
        return statistics.mean([self.get_node_by_hash(n[0]).get_node_coverage() for n in path[1:-1]])

    def process_path_pair(self, paths, all_bubble_nodes):
        sorted_by_coverage = sorted(list(paths), key=lambda x: x[1], reverse=True)
        highest_coverage_path = sorted_by_coverage[0][0]
        forward_high_coverage_path_genes = self.get_genes_in_unitig(highest_coverage_path)
        backward_high_coverage_path_genes = self.reverse_list_of_genes(forward_high_coverage_path_genes)
        paths_to_correct = self.identify_paths_to_correct(sorted_by_coverage, forward_high_coverage_path_genes, all_bubble_nodes)
        #self.correct_paths(paths_to_correct, forward_high_coverage_path_genes, backward_high_coverage_path_genes)

    def split_into_chunks(self, path, all_bubble_nodes):
        chunks = []
        this_chunk = []
        for n in path:
            if n not in all_bubble_nodes:
                this_chunk.append(n)
            else:
                if this_chunk:  # Only add non-empty chunks
                    chunks.append(this_chunk)
                    this_chunk = []
        if this_chunk:  # Add the last chunk if not empty
            chunks.append(this_chunk)
        return chunks

    def get_position_of_genes_in_alignment(self, alignment):
        false_index = 0
        true_index = 0
        false_path_alignment_to_gene, true_path_alignment_to_gene = {}, {}
        for c in range(len(alignment)):
            col = alignment[c]
            if not col[1] == "*":
                false_path_alignment_to_gene[false_index] = c
                false_index += 1
            if not col[0] == "*":
                true_path_alignment_to_gene[true_index] = c
                true_index += 1
        return false_path_alignment_to_gene, true_path_alignment_to_gene

    def get_gene_to_node_mapping(self, path):
        gene_to_node_mapping = {}
        for i in range(len(path)):
            gene_indices = [j for j in range(i, i + self.get_kmerSize())]
            for g in gene_indices:
                if g not in gene_to_node_mapping:
                    gene_to_node_mapping[g] = []
                gene_to_node_mapping[g].append(i)
        return gene_to_node_mapping

    def make_node_alignment(self,
                        true_path_alignment_to_gene,
                        false_path_alignment_to_gene,
                        true_gene_to_node_mapping,
                        false_gene_to_node_mapping,
                        true_path,
                        false_path,
                        alignment_length):
        node_alignment = []
        for i in range(alignment_length-self.get_kmerSize() + 1):
            if i in true_path_alignment_to_gene and i in false_path_alignment_to_gene:
                assert len(true_gene_to_node_mapping[true_path_alignment_to_gene[i]]) == len(false_gene_to_node_mapping[false_path_alignment_to_gene[i]])
                node_alignment.append((true_path[true_gene_to_node_mapping[true_path_alignment_to_gene[i]][-1]], false_path[false_gene_to_node_mapping[false_path_alignment_to_gene[i]][-1]]))
            if i in true_path_alignment_to_gene and i not in false_path_alignment_to_gene:
                #for j in range(len(true_gene_to_node_mapping[true_path_alignment_to_gene[i]])):
                node_alignment.append((true_path[true_gene_to_node_mapping[true_path_alignment_to_gene[i]][-1]], "*"))
            if i not in true_path_alignment_to_gene and i in false_path_alignment_to_gene:
                #for j in range(len(false_path_alignment_to_gene[false_path_alignment_to_gene[i]])):
                node_alignment.append(("*", false_path[false_gene_to_node_mapping[false_path_alignment_to_gene[i]][-1]]))
        return node_alignment

    def convert_gene_alignment_to_node_alignment(self, alignment, heaviest_path, path_to_correct):
        # get a mapping of the indices of genes in the path we are correcting to the alignment column
        false_path_alignment_to_gene, true_path_alignment_to_gene = self.get_position_of_genes_in_alignment(alignment)
        # get a mapping of the indices of genes in the path we are correction to the nodes containing them
        false_gene_to_node_mapping = self.get_gene_to_node_mapping(path_to_correct)
        true_gene_to_node_mapping= self.get_gene_to_node_mapping(heaviest_path)
        # make the node alignment
        node_alignment = self.make_node_alignment(true_path_alignment_to_gene,
                                                false_path_alignment_to_gene,
                                                true_gene_to_node_mapping,
                                                false_gene_to_node_mapping,
                                                heaviest_path,
                                                path_to_correct,
                                                len(alignment))
        return node_alignment

    def identify_paths_to_correct(self, sorted_paths, forward_genes, all_bubble_nodes):
        reads_to_correct = set()
        correction_dict = {}
        for other_path in sorted_paths[1:]:
            reads = self.collect_reads_in_path(other_path[0])
            if not all(n in self.get_nodes() for n in other_path[0]):
                continue
            genes_in_path = self.get_genes_in_unitig(other_path[0])
            alignment = self.needleman_wunsch(forward_genes, genes_in_path)
            difference_count = len([col for col in alignment if col[0] != col[1]])
            if difference_count <= 2:
                # convert the alignment into node space
                node_alignment = self.convert_gene_alignment_to_node_alignment(alignment, sorted_paths[0][0], other_path[0])
                # replace the nodes on the read in the false path
                for col in node_alignment:
                    if col[0] != col[1] and col[1] != "*":
                        if col[1] in self.get_nodes():
                            self.remove_node(self.get_node_by_hash(col[1]))
        return reads_to_correct, correction_dict

    def correct_paths(self, paths_to_correct, forward_genes, backward_genes):
        fw_counter = Counter(forward_genes)
        bw_counter = Counter(backward_genes)
        for path in paths_to_correct:
            reads_in_path = self.collect_reads_in_path(path)
            for read_id in reads_in_path:
                genes_on_read = self.get_reads()[read_id]
                shared_fw = len(list((fw_counter & Counter(genes_on_read)).elements()))
                shared_bw = len(list((bw_counter & Counter(genes_on_read)).elements()))
                assert shared_fw != 0 or shared_bw != 0, forward_genes
                if shared_fw > shared_bw:
                    alignment = self.needleman_wunsch(genes_on_read, forward_genes)
                elif shared_bw > shared_fw:
                    alignment = self.needleman_wunsch(genes_on_read, backward_genes)
                else:
                    print(genes_on_read)
                    print(forward_genes)
                    raise ValueError("Forward and backward lists of genes share the same number of common elements with the highest coverage path")
                if read_id == "SRR23044226.121397_0":
                    print([self.get_gene_mer_label(self.get_node_by_hash(h)) for h in path])
                    print(alignment)
                new_genes = self.get_new_genes_from_alignment(alignment)
                self.get_reads()[read_id] = new_genes

    def collect_reads_in_path(self, path):
        reads_in_path = set()
        for node_hash in list(path):
            if not node_hash in self.get_nodes():
                continue
            for read in self.get_node_by_hash(node_hash).get_reads():
                reads_in_path.add(read)
        return reads_in_path

    def get_new_genes_from_alignment(self, alignment):
        new_genes = []
        for col in alignment:
            if not col[1] == "*":
                new_genes.append(col[1])
            else:
                new_genes.append(col[0])
        return new_genes

    def get_direction_between_two_nodes(self, source_node_hash, target_node_hash):
        # get the edges between the target and the source
        source_to_target_edge, target_to_source_edge = self.get_edges_between_nodes(self.get_node_by_hash(source_node_hash), self.get_node_by_hash(target_node_hash))
        # get the direction from source to target
        source_to_target_direction = source_to_target_edge.get_targetNodeDirection()
        return source_to_target_direction * -1

    def newest_correct_bubbles(self):
        # get the potential start nodes
        potential_bubble_starts, all_bubble_nodes = self.identify_potential_bubble_starts()
        # define a maximum path length
        for component in self.components():
            # convert the bubbles nodes to a list
            potential_bubble_starts_component = potential_bubble_starts[component]
            # collect the bubbles
            bubbles = {}
            # iterate through the bubble nodes
            for start_hash, start_direction in potential_bubble_starts_component:
                # get the start node
                start_node = self.get_node_by_hash(start_hash)
                # get the neighbour nodes
                if start_direction == 1:
                    neighbors = [self.get_edge_by_hash(e).get_targetNode() for e in start_node.get_forward_edge_hashes()]
                if start_direction == -1:
                    neighbors = [self.get_edge_by_hash(e).get_targetNode() for e in start_node.get_backward_edge_hashes()]
                # get the unitigs for the neighbors
                for neighbor in neighbors:
                    unitig = self.get_linear_path_for_node(neighbor, True)
                    if len(unitig) > 2:
                        terminals = tuple(sorted([unitig[0], unitig[-1]]))
                        if not terminals in bubbles:
                            bubbles[terminals] = set()
                        path_coverage = self.calculate_path_coverage(unitig)
                        fw_path = (tuple(unitig), path_coverage)
                        rv_path = (tuple(list(reversed(unitig))), path_coverage)
                        bubbles[terminals].add(tuple(sorted([fw_path, rv_path])[0]))
            self.new_correct_bubble_paths(bubbles)
        return self.get_reads()

    def split_paths_by_coverage(self, paths, node_min_coverage):
        valid = []
        invalid = []
        for p in paths:
            if p[1] >= node_min_coverage:
                valid.append(p)
            else:
                invalid.append(p)
        return valid, invalid

    def new_correct_bubble_paths(self, bubbles):
        for pair in bubbles:
            if len(bubbles[pair]) > 1:
                # sort the paths from highest to lower coverage
                paths = sorted(list(bubbles[pair]), key=lambda x: x[1], reverse=True)
                # split the paths into valid and invalid paths
                valid_paths, invalid_paths = self.split_paths_by_coverage(paths, 10)
                if valid_paths == []:
                    valid_paths = [invalid_paths[0]]
                    invalid_paths = invalid_paths[1:]
                # store the corrected paths
                corrected_paths = set()
                # iterate through the valid paths
                for path, path_coverage in valid_paths:
                    path = [n[0] for n in path]
                    # get the genes in the path
                    genes_in_path = self.get_genes_in_unitig(list(path))
                    # make a counter for the genes in the path
                    gene_counts = Counter(genes_in_path)
                    # iterate through the invalid paths
                    for other_path, other_path_coverage in invalid_paths:
                        other_path = tuple([n[0] for n in other_path])
                        if not other_path in corrected_paths:
                            # get the genes in the path
                            fw_genes_in_path = self.get_genes_in_unitig(list(other_path))
                            # reverse the genes in the path
                            bw_genes_in_path = self.reverse_list_of_genes(fw_genes_in_path)
                            # decide whether the forward or backward path is closest and align the genes in the paths
                            if len(gene_counts & Counter(fw_genes_in_path)) > len(gene_counts & Counter(bw_genes_in_path)):
                                alignment = self.needleman_wunsch(genes_in_path, fw_genes_in_path)
                                other_genes_in_path = fw_genes_in_path
                            else:
                                alignment = self.needleman_wunsch(genes_in_path, bw_genes_in_path)
                                other_genes_in_path = bw_genes_in_path
                            # get the SNP differences
                            SNP_count = len([col for col in alignment if col[0] != col[1] and col[0] != "*" and col[1] != "*"])
                            # get the INDEL differences
                            INDEL_count = len([col for col in alignment if col[0] != col[1] and (col[0] == "*" or col[1] == "*")])
                            if SNP_count <= 1 and INDEL_count <= 2:
                                # get the indices mapping the original lists of genes to alignment columns
                                false_path_alignment_to_gene, true_path_alignment_to_gene = self.get_position_of_genes_in_alignment(alignment)
                                # get the sublists we are going to replace
                                replacement_dict = {}
                                reversed_replacement_dict = {}
                                for i in range(len(other_genes_in_path) - self.get_kmerSize() + 1):
                                    alignment_indices = [false_path_alignment_to_gene[j] for j in range(i, i+self.get_kmerSize())]
                                    alignment_columns = alignment[min(alignment_indices): max(alignment_indices) + 1]
                                    replacement = [col[0] for col in alignment_columns if not col[0] == "*"]
                                    to_replace = other_genes_in_path[i: i + self.get_kmerSize()]
                                    replacement_dict["~~~".join(to_replace)] = "~~~".join(replacement)
                                    reversed_replacement_dict["~~~".join(self.reverse_list_of_genes(to_replace))] = "~~~".join(self.reverse_list_of_genes(replacement))
                                # get the reads in the path
                                nodes_to_correct = [n for n in list(other_path)[1:-1]]
                                reads_in_path = self.collect_reads_in_path(nodes_to_correct)
                                for read_id in reads_in_path:
                                    genes_on_read = self.get_reads()[read_id][:]
                                    # convert the list of genes to a string
                                    genes_on_read = "~~~".join(genes_on_read)
                                    fw_count = 0
                                    rv_count = 0
                                    # count the number of keys in the list of genes
                                    for k in replacement_dict:
                                        if k in genes_on_read:
                                            fw_count += genes_on_read.count(k)
                                    for k in reversed_replacement_dict:
                                        if k in genes_on_read:
                                            rv_count += genes_on_read.count(k)
                                    # Determine which replacement dict to use based on shared counts
                                    replacement = replacement_dict if fw_count > rv_count else reversed_replacement_dict
                                    # replace the genes on the read
                                    for to_replace in replacement:
                                        genes_on_read = genes_on_read.replace(to_replace, replacement[to_replace])
                                    if read_id == "SRR23044210.27420_0":
                                        print(self.get_reads()[read_id][:])
                                        print(genes_on_read.split("~~~"))
                                        print(replacement,"\n")
                                    #     djddjdjk
                                    # Correct the read using the appropriate replacement dictionary
                                    self.get_reads()[read_id] = genes_on_read.split("~~~")
                                    corrected_paths.add(other_path)

    def correct_bubbles(self, gap_distances, node_min_coverage):
        # get the potential start nodes
        potential_bubble_starts, all_bubble_nodes = self.identify_potential_bubble_starts()
        # define a maximum distance
        max_distance = self.get_kmerSize() * 3
        # iterate through the components
        for component in self.components():
            #if not component == 2:
            #    continue
            # convert the bubbles nodes to a list
            if not component in potential_bubble_starts:
                continue
            potential_bubble_starts_component = potential_bubble_starts[component]
            # store the paths
            unique_paths = set()
            seen_pairs = set()
            # iterate through the bubble nodes
            for start_hash, start_direction in potential_bubble_starts_component:
                for stop_hash, stop_direction in potential_bubble_starts_component:
                    if start_hash == stop_hash:
                        continue
                    if tuple(sorted([start_hash, stop_hash])) in seen_pairs:
                        continue
                    seen_pairs.add(tuple(sorted([start_hash, stop_hash])))
                    paths_between_nodes = self.new_find_paths_between_nodes(start_hash, stop_hash, max_distance, start_direction)
                    # remove paths that do not start and stop in the same direction
                    valid_paths_between_nodes = [p for p in paths_between_nodes if p[0] == (start_hash, start_direction) and (p[-1][0], self.get_direction_between_two_nodes(p[-2][0], p[-1][0])) == (stop_hash, stop_direction)]
                    # check if there is more than one option
                    if len(valid_paths_between_nodes) > 1:
                        for p in valid_paths_between_nodes:
                            unique_paths.add(tuple(sorted([p, list(reversed([(t[0], t[1]*-1) for t in p]))])[0]))
            # convert to list
            unique_paths = list(unique_paths)
            # Assuming unique_paths is a list of lists or a similar iterable of iterables.
            filtered_paths = []
            for p in unique_paths:
                p_list = list(p)
                # Skip self-comparison and check for sublists in both forward and reversed directions
                if not any(self.is_sublist(list(q), p) or self.is_sublist(list(q), list(reversed(p))) for q in unique_paths if q != p):
                    filtered_paths.append((p_list, self.calculate_path_coverage(p_list)))
            # sort by coverage
            sorted_filtered_paths = sorted(filtered_paths, key=lambda x:x[1], reverse=True)
            paired_sorted_filtered_paths = {}
            junctions_crossed = set()
            for p in sorted_filtered_paths:
                sorted_terminals = tuple(sorted([p[0][0][0], p[0][-1][0]]))
                if not sorted_terminals in paired_sorted_filtered_paths:
                    paired_sorted_filtered_paths[sorted_terminals] = []
                paired_sorted_filtered_paths[sorted_terminals].append(p)
                # see which junctions this path crosses
                juncs = [n[0] for n in p[0][1:-1] if n[0] in all_bubble_nodes[component]]
                junctions_crossed.update(juncs)
            # refine further
            refined_pairs = {}
            for pair in paired_sorted_filtered_paths:
                if pair[0] not in junctions_crossed or pair[1] not in junctions_crossed:
                    refined_pairs[pair] = paired_sorted_filtered_paths[pair]
            # clean the paths
            self.new_correct_bubble_paths(refined_pairs)
        return self.get_reads()

    def correct_read(self, alignment_map, read):
        corrected_read = []
        for gene in read:
            corrected_gene = alignment_map.get(gene)
            if corrected_gene:
                corrected_read.append(corrected_gene)
            else:
                continue
        return corrected_read

    def get_node_for_genes(self, listOfGenes):
        gene_objects = [Gene(g) for g in listOfGenes]
        return Node(GeneMer(gene_objects)).__hash__()

    def new_correct_bubble_paths(self, bubbles):
        for pair in bubbles:
            if len(bubbles[pair]) > 1:
                # sort the paths from highest to lower coverage
                paths = sorted(list(bubbles[pair]), key=lambda x: x[1], reverse=True)
                # split the paths into valid and invalid paths
                valid_paths, invalid_paths = self.split_paths_by_coverage(paths, 10)
                if valid_paths == []:
                    valid_paths = [invalid_paths[0]]
                    invalid_paths = invalid_paths[1:]
                # store the corrected paths
                corrected_paths = set()
                # iterate through the valid paths
                for path, path_coverage in valid_paths:
                    path = [n[0] for n in path]
                    # get the genes in the path
                    genes_in_path = self.get_genes_in_unitig(list(path))
                    # make a counter for the genes in the path
                    gene_counts = Counter(genes_in_path)
                    # iterate through the invalid paths
                    for other_path, other_path_coverage in invalid_paths:
                        other_path = tuple([n[0] for n in other_path])
                        if not other_path in corrected_paths:
                            # get the genes in the path
                            fw_genes_in_path = self.get_genes_in_unitig(list(other_path))
                            # reverse the genes in the path
                            bw_genes_in_path = self.reverse_list_of_genes(fw_genes_in_path)
                            # decide whether the forward or backward path is closest and align the genes in the paths
                            if len(gene_counts & Counter(fw_genes_in_path)) > len(gene_counts & Counter(bw_genes_in_path)):
                                alignment = self.needleman_wunsch(genes_in_path, fw_genes_in_path)
                                other_genes_in_path = fw_genes_in_path
                            else:
                                alignment = self.needleman_wunsch(genes_in_path, bw_genes_in_path)
                                other_genes_in_path = bw_genes_in_path
                            # get the SNP differences
                            SNP_count = len([col for col in alignment if col[0] != col[1] and col[0] != "*" and col[1] != "*"])
                            # get the INDEL differences
                            INDEL_count = len([col for col in alignment if col[0] != col[1] and (col[0] == "*" or col[1] == "*")])
                            if SNP_count <= 1 and INDEL_count <= 2:
                                # get the indices mapping the original lists of genes to alignment columns
                                false_path_alignment_to_gene, true_path_alignment_to_gene = self.get_position_of_genes_in_alignment(alignment)
                                # get the sublists we are going to replace
                                replacement_dict = {}
                                for i in range(len(other_genes_in_path) - self.get_kmerSize() + 1):
                                    alignment_indices = [false_path_alignment_to_gene[j] for j in range(i, i+self.get_kmerSize())]
                                    alignment_columns = alignment[min(alignment_indices): max(alignment_indices) + 1]
                                    replacement = [col[0] for col in alignment_columns if not col[0] == "*"]
                                    to_replace = other_genes_in_path[i: i + self.get_kmerSize()]
                                    replacement_dict[self.get_node_for_genes(to_replace)] = self.get_node_for_genes(replacement)
                                print(replacement_dict)
                                # get the reads in the path
                                nodes_to_correct = [n for n in list(other_path)[1:-1]]
                                reads_in_path = self.collect_reads_in_path(nodes_to_correct)
                                for read_id in reads_in_path:
                                    genes_on_read = self.get_reads()[read_id][:]
                                    # convert the list of genes to a string
                                    genes_on_read = "~~~".join(genes_on_read)
                                    # Determine which replacement dict to use based on shared counts
                                    # replace the genes on the read
                                    for to_replace in replacement:
                                        genes_on_read = genes_on_read.replace(to_replace, replacement_dict[to_replace])
                                    # Correct the read using the appropriate replacement dictionary
                                    self.get_reads()[read_id] = genes_on_read.split("~~~")
                                    corrected_paths.add(other_path)

    def identify_potential_bubble_starts(self):
        potential_bubble_starts = {}
        all_bubble_nodes = {}
        for node in self.all_nodes():
            fw_edges = node.get_forward_edge_hashes()  # Forward edges
            bw_edges = node.get_backward_edge_hashes()  # Backward edges
            if len(fw_edges) > 1 and len(bw_edges) < 2:
                if not node.get_component() in potential_bubble_starts:
                    potential_bubble_starts[node.get_component()] = []
                    all_bubble_nodes[node.get_component()] = set()
                potential_bubble_starts[node.get_component()].append((node.__hash__(), 1))
                all_bubble_nodes[node.get_component()].add(node.__hash__())
            if len(bw_edges) > 1 and len(fw_edges) < 2:
                if not node.get_component() in potential_bubble_starts:
                    potential_bubble_starts[node.get_component()] = []
                    all_bubble_nodes[node.get_component()] = set()
                potential_bubble_starts[node.get_component()].append((node.__hash__(), -1))
                all_bubble_nodes[node.get_component()].add(node.__hash__())
        return potential_bubble_starts, all_bubble_nodes

    def find_potential_paths(self, start, all_bubble_nodes, max_distance):
        start_node, direction = start
        start_node_hash = start_node.__hash__()
        # get all paths of length max distance from the start node
        paths = self.new_find_paths_between_nodes(start_node_hash, None, max_distance, direction)
        # slice the paths to the last element that is the start of a bubble
        valid_paths = set()
        for p in paths:
            junction_nodes_in_path = [(i, v) for i, v in enumerate(p) if v[0] in all_bubble_nodes]
            index = max([t[0] for t in junction_nodes_in_path], default=-1)
            assert not index == -1
            sliced = p[:index + 1]
            if not len(sliced) == 0:
                valid_paths.add(tuple(sliced))
        paths_from_start = {}
        for p in valid_paths:
            p = list(p)
            if 2 < len(p):
                terminals = (p[0][0], p[-1][0])
                if not terminals in paths_from_start:
                    paths_from_start[terminals] = []
                path_coverage = self.calculate_path_coverage(p)
                paths_from_start[terminals].append((([n[0] for n in p]), path_coverage))
        return paths_from_start

    def new_find_paths_between_nodes(self, start_hash, end_hash, distance, current_direction, path=None, seen_nodes=None):
        if path is None:
            path = []
        if seen_nodes is None:
            seen_nodes = set()

        path.append((start_hash, current_direction))  # Append the current node and direction to the path
        seen_nodes.add(start_hash)  # Add the current node to the seen nodes to avoid revisiting

        # Check if we've reached the end node or if end_hash is None and the path length equals the specified distance
        if (end_hash and start_hash == end_hash and len(path) - 1 <= distance) or (end_hash is None and len(path) - 1 == distance):
            path_copy = path.copy()  # Make a copy of the path to return
            path.pop()  # Remove the last element before returning to ensure correct backtracking
            return [path_copy]  # Return the path copy

        if len(path) - 1 > distance:  # Check if the current path length exceeds the allowed distance
            path.pop()  # Remove the last element added to path before returning
            return []

        paths = []  # Initialize the list to collect paths
        next_nodes = []
        if current_direction == 1:
            next_nodes = self.get_node_by_hash(start_hash).get_forward_edge_hashes()
        elif current_direction == -1:
            next_nodes = self.get_node_by_hash(start_hash).get_backward_edge_hashes()
        for edge_hash in next_nodes:
            edge = self.get_edge_by_hash(edge_hash)
            node_hash = edge.get_targetNode().__hash__()
            if node_hash not in seen_nodes:  # Check if the node has not been visited
                # Make a copy of seen_nodes for each recursive call to ensure correct tracking of visited nodes
                new_seen_nodes = seen_nodes.copy()
                new_seen_nodes.add(node_hash)
                newpaths = self.new_find_paths_between_nodes(
                    node_hash,
                    end_hash,
                    distance,
                    edge.get_targetNodeDirection(),
                    path.copy(),  # Pass a copy of the path to avoid modification by reference
                    new_seen_nodes,  # Pass the updated copy of seen_nodes
                )
                paths.extend(newpaths)  # Add the new paths found to the paths list

        path.pop()  # Remove the last element added to path before returning to the caller
        return paths
