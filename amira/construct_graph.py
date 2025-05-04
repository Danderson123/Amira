import os
import statistics
import sys
from collections import Counter, defaultdict, deque
from itertools import product
from multiprocessing import Pool

import numpy as np
import sourmash
from joblib import Parallel, delayed
from suffix_tree import Tree
from tqdm import tqdm

from amira.construct_edge import Edge
from amira.construct_gene import convert_int_strand_to_string
from amira.construct_gene_mer import GeneMer
from amira.construct_node import Node
from amira.construct_read import Read
from amira.path_finding_utils import (
    construct_suffix_tree,
    filter_blocks,
    get_suffixes_from_initial_tree,
    process_anchors,
    process_combinations_for_i,
)

sys.setrecursionlimit(50000)


class GeneMerGraph:
    def __init__(self, readDict, kmerSize, gene_positions=None):
        self._reads = readDict
        self._kmerSize = kmerSize
        self._minNodeCoverage = 1
        self._minEdgeCoverage = 1
        self._genePositions = gene_positions
        self._nodes = {}
        self._edges = {}
        self._readNodes = {}
        self._readNodeDirections = {}
        self._readNodePositions = {}
        self._shortReads = {}
        self._readsToCorrect = set()
        # initialise the graph
        for readId in self.get_reads():
            if self.get_gene_positions():
                read = Read(readId, self.get_reads()[readId], self.get_gene_positions()[readId])
            else:
                read = Read(readId, self.get_reads()[readId], None)
            # get the gene mers on this read
            geneMers, geneMerPositions = read.get_geneMers(self.get_kmerSize())
            # if there are no gene mers there is nothing to add to the graph
            if len(geneMers) == 0:
                self._shortReads[readId] = read.get_annotatedGenes()
                continue
            if not len(geneMers) == 1:
                # iterate through all the gene-mers except the last
                for g in range(len(geneMers) - 1):
                    # get the read id for this read
                    this_readId = read.get_readId()
                    # add a source node to the graph
                    sourceNode = self.add_node(geneMers[g], [this_readId])
                    # add the sourceNode to the read
                    self.add_node_to_read(
                        sourceNode,
                        this_readId,
                        geneMers[g].get_geneMerDirection(),
                        geneMerPositions[g],
                    )
                    # increase the source node coverage by 1
                    sourceNode.increment_node_coverage()
                    # get the target node readId
                    targetReadId = read.get_readId()
                    # add the target node to the graph
                    targetNode = self.add_node(geneMers[g + 1], [targetReadId])
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
                self.add_node_to_read(
                    targetNode,
                    this_readId,
                    geneMers[-1].get_geneMerDirection(),
                    geneMerPositions[-1],
                )
            else:
                # add a single node to the graph if there is only 1 gene mer
                this_readId = read.get_readId()
                sourceNode = self.add_node(geneMers[0], [this_readId])
                self.add_node_to_read(
                    sourceNode, this_readId, geneMers[0].get_geneMerDirection(), geneMerPositions[0]
                )
                sourceNode.increment_node_coverage()
        # assign component Ids to all nodes in the graph
        self.assign_component_ids()

    def get_reads(self) -> dict[str, list[str]]:
        """return a dictionary of all reads and their genes"""
        return self._reads

    def get_short_read_annotations(self):
        return self._shortReads

    def get_gene_positions(self):
        return self._genePositions

    def get_short_read_gene_positions(self):
        return {r: self.get_gene_positions()[r] for r in self.get_short_read_annotations()}

    def get_readNodes(self) -> dict[str, list[str]]:
        """return a dictionary of all reads and their node hashes"""
        return self._readNodes

    def get_readNodeDirections(self) -> dict[str, list[int]]:
        """return a dictionary of all reads and their node directions"""
        return self._readNodeDirections

    def get_readNodePositions(self) -> dict[str, list[tuple[int, int]]]:
        return self._readNodePositions

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

    def add_node_to_read(
        self, node: Node, readId: str, node_direction: int, node_position=None
    ) -> list[int]:
        # add the read ID to the read dict if it is not present
        if readId not in self.get_readNodes():
            self.get_readNodes()[readId] = []
            self.get_readNodeDirections()[readId] = []
            self.get_readNodePositions()[readId] = []
        # add the hash for the node as attributes for this read
        self.get_readNodes()[readId].append(node.__hash__())
        self.get_readNodeDirections()[readId].append(node_direction)
        self.get_readNodePositions()[readId].append(node_position)
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
        # add the reads for this node to the dictionary of reads
        for r in reads:
            node.add_read(r)
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
        if edgeHash not in self.get_edges():
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
            new_positions = []
            for i in range(len(self.get_readNodes()[readId])):
                if self.get_readNodes()[readId][i] != node_to_remove.__hash__():
                    new_nodes.append(self.get_readNodes()[readId][i])
                    new_directions.append(self.get_readNodeDirections()[readId][i])
                    new_positions.append(self.get_readNodePositions()[readId][i])
                else:
                    new_nodes.append(None)
                    new_directions.append(None)
                    new_positions.append(None)
            # keep track of all of the reads we have modified
            self.get_reads_to_correct().add(readId)
            self.get_readNodes()[readId] = new_nodes
            self.get_readNodeDirections()[readId] = new_directions
            self.get_readNodePositions()[readId] = new_positions

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

    def write_edge_entry(
        self, source_node, target_node, source_edge_direction, target_edge_direction, edge_coverage
    ):
        """return a string of a gml edge entry"""
        edge_entry = "\tedge\t[\n"
        edge_entry += "\t\tsource\t" + str(source_node) + "\n"
        edge_entry += "\t\ttarget\t" + str(target_node) + "\n"
        edge_entry += "\t\tsource_direction\t" + str(source_edge_direction) + "\n"
        edge_entry += "\t\ttarget_direction\t" + str(target_edge_direction) + "\n"
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

    def get_genes_in_unitig(self, listOfNodes):
        if len(listOfNodes) == 1:
            return self.get_gene_mer_genes(self.get_node_by_hash(listOfNodes[0]))
        # store the new annotations in a list
        newAnnotations = []
        # iterate through the new readNodes
        errored = False
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
            fw_genes = self.get_gene_mer_genes(targetNode)
            bw_genes = self.get_reverse_gene_mer_genes(targetNode)
            # Check if forward genes match, if not try reverse
            try:
                if fw_genes[:-1] == newAnnotations[-self.get_kmerSize() + 1 :]:
                    newAnnotations.append(fw_genes[-1])
                else:
                    assert bw_genes[:-1] == newAnnotations[-self.get_kmerSize() + 1 :]
                    newAnnotations.append(bw_genes[-1])
            except AssertionError:
                errored = True
                break

        # Handle errors and try the alternative path
        if errored:
            newAnnotations = []
            for n in range(len(listOfNodes) - 1):
                sourceNode = self.get_node_by_hash(listOfNodes[n])
                targetNode = self.get_node_by_hash(listOfNodes[n + 1])
                # Get the edge between the source and target node
                edgeHashes = self.get_edge_hashes_between_nodes(sourceNode, targetNode)
                edge = self.get_edge_by_hash(edgeHashes[0])
                # Add either forward or reverse gene mer depending on first node direction
                if n == 0:
                    newAnnotations += (
                        self.get_gene_mer_genes(sourceNode)
                        if edge.get_sourceNodeDirection() == 1
                        else self.get_reverse_gene_mer_genes(sourceNode)
                    )
                # Get forward and backward genes for target node
                fw_genes = self.get_gene_mer_genes(targetNode)
                bw_genes = self.get_reverse_gene_mer_genes(targetNode)
                # Check if forward or backward genes match
                try:
                    if fw_genes[1:] == newAnnotations[: self.get_kmerSize() - 1]:
                        newAnnotations.insert(0, fw_genes[0])
                    else:
                        assert bw_genes[1:] == newAnnotations[: self.get_kmerSize() - 1]
                        newAnnotations.insert(0, bw_genes[0])
                except AssertionError:
                    raise ValueError("Gene sequences do not match in alternative path.")
        return newAnnotations

    def remove_short_linear_paths(self, min_length, sample_genesOfInterest={}):
        """remove nodeHashes on reads if on a linear path of length < min_length. \
            Returns a list of nodeHashes that have bene removed"""
        paths_to_remove = {}
        for node in self.all_nodes():
            if self.get_degree(node) == 1:
                path = self.get_linear_path_for_node(node)
                if len(path) > 0 and (
                    len(path) < min_length
                ):  # or all(self.get_node_by_hash(n).get_node_coverage() < 2 for n in path)):
                    if node.get_component() not in paths_to_remove:
                        paths_to_remove[node.get_component()] = []
                    paths_to_remove[node.get_component()].append(path)
        # get the nodes in the graph containing AMR genes
        AMR_nodes = self.get_AMR_nodes(sample_genesOfInterest)
        # remove paths that do not contain AMR genes
        removed = set()
        for component in paths_to_remove:
            if component is not None:
                nodes_in_component = set(
                    [n.__hash__() for n in self.get_nodes_in_component(component)]
                )
            else:
                nodes_in_component = []
            for path in paths_to_remove[component]:
                if component is not None and len(nodes_in_component.intersection(path)) == len(
                    nodes_in_component
                ):
                    continue
                for nodeHash in path:
                    if nodeHash in AMR_nodes:
                        continue
                    if nodeHash not in removed:
                        self.remove_node(self.get_node_by_hash(nodeHash))
                        removed.add(nodeHash)
        return list(removed)

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
        return statistics.mean(node_coverages)

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
                    edge.get_sourceNodeDirection(),
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

    def remove_low_coverage_components(self, min_component_coverage):
        # remove components that have a coverage less than min_component_coverage
        for component_ID in self.components():
            nodes_in_component = self.get_nodes_in_component(component_ID)
            if all(
                node.get_node_coverage() < min_component_coverage for node in nodes_in_component
            ):
                for node in nodes_in_component:
                    self.remove_node(node)

    def reverse_list_of_genes(self, list_of_genes: list[str]) -> list[str]:
        return [("-" if g[0] == "+" else "+") + g[1:] for g in reversed(list_of_genes)]

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
                        for p in self.find_paths(matrix, sorted_indices[0], sorted_indices[1])
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
                n for n in self.get_forward_neighbors(node) if n.__hash__() in nodeHashesOfInterest
            ]
            backwardAMRNodes = [
                n for n in self.get_backward_neighbors(node) if n.__hash__() in nodeHashesOfInterest
            ]
            if len(backwardAMRNodes) == 0 or len(forwardAMRNodes) == 0:
                nodeAnchors.add(node_hash)
            if (
                len(self.get_backward_neighbors(node)) > 1
                or len(self.get_forward_neighbors(node)) > 1
            ):
                nodeJunctions.add(node_hash)
        return nodeAnchors, nodeJunctions

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

    def correct_reads(self, fastq_data):
        """Main function to correct reads in the dataset."""
        readNodes = self.get_readNodes()  # Retrieve read nodes
        corrected_genes = {}
        corrected_gene_positions = {}
        for read_id in readNodes:
            list_of_genes = self.correct_single_read(read_id, readNodes, fastq_data)
            if len(list_of_genes) > 0:
                corrected_genes[read_id] = list_of_genes
                if self.get_gene_positions():
                    corrected_gene_positions[read_id] = self.get_gene_positions()[read_id]
        return corrected_genes, corrected_gene_positions

    def correct_single_read(self, read_id, readNodes, fastq_data):
        """Correct a single read based on its ID and the readNodes data."""
        if read_id not in self.get_reads_to_correct():
            return self.get_reads()[read_id]

        if not all(n is None for n in readNodes[read_id]):
            start, end = self.find_read_boundaries(readNodes[read_id])
            # get the new list of genes on each read
            new_genes_on_read = self.process_read_correction(
                read_id, readNodes, start, end, fastq_data
            )
            if self.get_gene_positions():
                assert len(new_genes_on_read) == len(self.get_gene_positions()[read_id])
            return new_genes_on_read
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
        # # get all the possible upstream paths
        # if start != 0:
        #     upstream_paths = self.find_paths_between_nodes(
        #         nodes_on_read[start][0], None, start * 2, -nodes_on_read[start][1]
        #     )
        # # get all the possible downstream paths
        # if end != len(nodes_on_read) - 1:
        #     downstream_paths = self.find_paths_between_nodes(
        #         nodes_on_read[end][0],
        #         None,
        #         (len(nodes_on_read) - end - 1) * 2,
        #         nodes_on_read[end][1],
        #     )
        possible_paths = []
        # for now we are not going to add the upstream and downstream paths
        upstream_paths, downstream_paths = [], []
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

                    # Combine the paths
                    # excludes the last node of the upstream and the first node of the downstream
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

    def process_read_correction(self, read_id, readNodes, start, end, fastq_data):
        nodes_on_read = [
            (readNodes[read_id][i], self.get_readNodeDirections()[read_id][i])
            for i in range(len(readNodes[read_id]))
        ]
        # correct filtered nodes in the middle of a read
        path_terminals = self.identify_path_terminals(readNodes[read_id], start, end)
        # chop the list of genes if the ends have been cut off
        if len(path_terminals) == 0:
            new_nodes = [n[0] for n in nodes_on_read[start : end + 1]]
            new_directions = [n[1] for n in nodes_on_read[start : end + 1]]
            if self.get_gene_positions():
                # modify the gene positions dictionary to account for removed nodes
                self.get_gene_positions()[read_id] = self.get_gene_positions()[read_id][
                    start : end + self.get_kmerSize()
                ]
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
        # align the new and old lists of genes
        alignment = self.needleman_wunsch(closest, self.get_reads()[read_id])
        # modify the gene positions
        current_index = 0
        new_positions = []
        for i in range(len(alignment)):
            col = alignment[i]
            if col[0] != "*":
                if col[1] != col[0]:
                    new_positions.append((None, None))
                else:
                    new_positions.append(self.get_gene_positions()[read_id][current_index])
                    current_index += 1
            else:
                current_index += 1
        # infer the missing gene positions
        new_positions = self.replace_invalid_gene_positions(new_positions, fastq_data, read_id)
        self.get_gene_positions()[read_id] = new_positions
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
        paths = self.new_find_paths_between_nodes(
            corrected[pair[0]][0],
            corrected[pair[1]][0],
            self.get_kmerSize() * 2,
            corrected[pair[0]][1],
        )
        return {pair: paths}

    def remove_junk_reads(self, error_rate):
        new_reads = {}
        new_positions = {}
        rejected_reads = {}
        rejected_read_positions = {}
        read_nodes = self.get_readNodes()
        reads = self.get_reads()
        gene_positions = self.get_gene_positions()
        for read_id in read_nodes:
            nodes = read_nodes[read_id]
            number_of_nodes = len(nodes)
            expected_filtered_nodes = round(number_of_nodes * (1 - error_rate))
            number_of_filtered_nodes = 0
            for n in nodes:
                if n is None:
                    number_of_filtered_nodes += 1
            if number_of_filtered_nodes <= expected_filtered_nodes:
                new_reads[read_id] = reads[read_id]
                new_positions[read_id] = gene_positions[read_id]
            else:
                rejected_reads[read_id] = reads[read_id]
                rejected_read_positions[read_id] = gene_positions[read_id]
        return new_reads, new_positions, rejected_reads, rejected_read_positions

    def get_valid_reads_only(self):
        valid_reads = {}
        for read_id in self.get_reads():
            if read_id not in self.get_reads_to_correct():
                valid_reads[read_id] = self.get_reads()[read_id]
        return valid_reads

    def score(self, a, b):
        # Scoring function: returns 1 if elements are equal, 0 otherwise
        return int(a == b)

    def needleman_wunsch(self, x, y):
        N, M = len(x), len(y)
        # Direction constants for traceback
        DIAG, LEFT, UP = (-1, -1), (-1, 0), (0, -1)
        # Initialize score (F) and pointer (Ptr) matrices
        F, Ptr = {}, {}
        F[-1, -1] = 0
        # Initial scoring for gaps along x
        for i in range(N):
            F[i, -1] = -i
        # Initial scoring for gaps along y
        for j in range(M):
            F[-1, j] = -j
        # Option for Ptr to trace back alignment
        option_Ptr = DIAG, LEFT, UP
        # Fill F and Ptr tables
        for i, j in product(range(N), range(M)):
            # Score options: match/mismatch, gap in x, gap in y
            option_F = (
                F[i - 1, j - 1] + self.score(x[i], y[j]),  # Match/mismatch
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

    def calculate_path_coverage(self, path):
        return statistics.mean(
            [self.get_node_by_hash(n[0]).get_node_coverage() for n in path[1:-1]]
        )

    def get_gene_to_node_mapping(self, path):
        gene_to_node_mapping = {}
        for i in range(len(path)):
            gene_indices = [j for j in range(i, i + self.get_kmerSize())]
            for g in gene_indices:
                if g not in gene_to_node_mapping:
                    gene_to_node_mapping[g] = []
                gene_to_node_mapping[g].append(i)
        return gene_to_node_mapping

    def collect_reads_in_path(self, path):
        reads_in_path = set()
        for node_hash in list(path):
            if node_hash not in self.get_nodes():
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
        source_to_target_edge, target_to_source_edge = self.get_edges_between_nodes(
            self.get_node_by_hash(source_node_hash), self.get_node_by_hash(target_node_hash)
        )
        # get the direction from source to target
        source_to_target_direction = source_to_target_edge.get_targetNodeDirection()
        return source_to_target_direction * -1

    def reverse_gene(self, gene):
        if gene[0] == "+":
            return "-" + gene[1:]
        if gene[0] == "-":
            return "+" + gene[1:]
        if gene[0] == "*":
            return "*"

    def reverse_gene_alignment(self, alignment):
        # reverse the alignment
        reversed_alignment = list(reversed(alignment))
        # reverse the gene directions
        reversed_alignment = [
            (self.reverse_gene(col[0]), self.reverse_gene(col[1])) for col in reversed_alignment
        ]
        return reversed_alignment

    def get_read_sequence_for_path(self, read_id, path, fastq_data):
        # check that the number of nodes and number of node positions is the same
        assert len(self.get_readNodes()[read_id]) == len(self.get_readNodePositions()[read_id])
        nodes_on_read = self.get_readNodes()[read_id]
        # Find the range of nodes on the read that are shared with the path
        shared_node_indices = [i for i, node_hash in enumerate(nodes_on_read) if node_hash in path]
        if not shared_node_indices:
            return None
        # get the sequence of the read
        entire_read_sequence = fastq_data[read_id]["sequence"]
        first_index, last_index = shared_node_indices[0], shared_node_indices[-1]
        start_pos = self.get_readNodePositions()[read_id][first_index][0]
        if start_pos is None:
            assert first_index == 0, self.get_readNodePositions()[read_id]
            start_pos = 0
        end_pos = self.get_readNodePositions()[read_id][last_index][1]
        if end_pos is None:
            assert (
                last_index == len(self.get_readNodePositions()[read_id]) - 1
            ), self.get_readNodePositions()[read_id]
            end_pos = len(entire_read_sequence) - 1
        assert start_pos >= 0
        assert end_pos < len(entire_read_sequence)
        read_sequence = entire_read_sequence[start_pos : end_pos + 1]
        return read_sequence

    def get_minhash_for_path(self, path, reads_in_path, fastq_data):
        minhash = sourmash.MinHash(n=0, ksize=9, scaled=1)
        reads_with_positions = set()
        for read_id in reads_in_path:
            read_sequence = self.get_read_sequence_for_path(read_id, path, fastq_data)
            if read_sequence is not None:
                minhash.add_sequence(read_sequence, force=True)
                reads_with_positions.add(f"{read_id}")
        return minhash, reads_with_positions

    def count_snps_in_alignment(self, aln):
        return len([col for col in aln if col[0] != col[1] and col[0] != "*" and col[1] != "*"])

    def count_indels_in_alignment(self, aln):
        return len([col for col in aln if col[0] != col[1] and (col[0] == "*" or col[1] == "*")])

    def get_gene_mer_strings(self, genes_on_read):
        gene_mers_on_read = []
        for i in range(len(genes_on_read) - (self.get_kmerSize() - 1)):
            # take a slice of the list of Genes from index i to i + kmerSize
            gene_mer = genes_on_read[i : i + self.get_kmerSize()]
            gene_mers_on_read.append(tuple(gene_mer))
        return gene_mers_on_read

    def reorient_alignment(
        self,
        gene_mers_on_read,
        fw_genes_in_path_counter,
        bw_genes_in_path_counter,
        fw_alignment,
        rv_alignment,
    ):
        # make a counter for the genes on the read
        genes_on_read_counter = Counter(gene_mers_on_read)
        # choose the alignment or rv alignment
        fw_count = len(genes_on_read_counter & fw_genes_in_path_counter)
        rv_count = len(genes_on_read_counter & bw_genes_in_path_counter)
        # Determine which replacement dict to use based on shared counts
        if fw_count > rv_count:
            return fw_alignment
        elif rv_count > fw_count:
            return rv_alignment
        elif fw_count == 0 and rv_count == 0:
            return None
        else:
            return None

    #            raise ValueError("Forward and reverse path alignments are equally distant.")

    def correct_genes_on_read(
        self,
        genes_on_read,
        first_shared_read_index,
        last_shared_read_index,
        alignment_subset,
        read_id,
    ):
        prefix = genes_on_read[:first_shared_read_index]
        suffix = genes_on_read[last_shared_read_index + 1 :]
        core = [c[0] for c in alignment_subset if c[0] != "*"]
        self.get_reads()[read_id] = prefix + core + suffix
        return self.get_reads()[read_id]

    def get_gene_position_prefix(self, gene_positions, first_shared_read_index):
        return gene_positions[:first_shared_read_index]

    def get_gene_position_suffix(self, gene_positions, last_shared_read_index):
        return gene_positions[last_shared_read_index + 1 :]

    def get_gene_position_core(
        self, gene_positions, first_shared_read_index, last_shared_read_index
    ):
        return gene_positions[first_shared_read_index : last_shared_read_index + 1]

    def get_new_gene_position_core(self, alignment_subset, core_gene_positions):
        core_position_index = 0
        new_core_gene_positions = []
        for i in range(len(alignment_subset)):
            col = alignment_subset[i]
            if col[0] != "*":
                if col[1] != col[0]:
                    new_core_gene_positions.append((None, None))
                else:
                    new_core_gene_positions.append(core_gene_positions[core_position_index])
                    core_position_index += 1
            else:
                core_position_index += 1
        return new_core_gene_positions

    def join_gene_position_ends_with_core(
        self, position_prefix, position_suffix, new_core_gene_positions
    ):
        if len(position_prefix) != 0 and len(position_suffix) != 0:
            new_positions = position_prefix + new_core_gene_positions + position_suffix
        elif len(position_prefix) != 0 and len(position_suffix) == 0:
            new_positions = position_prefix + new_core_gene_positions
        elif len(position_prefix) == 0 and len(position_suffix) != 0:
            new_positions = new_core_gene_positions + position_suffix
        else:
            new_positions = new_core_gene_positions
        return new_positions

    def replace_invalid_gene_positions(self, new_positions, fastq_data, read_id):
        prev_end = 0
        for i, (start, end) in enumerate(new_positions):
            if end is not None:
                prev_end = end
            if start is None and end is None:
                next_start = None
                for j in range(i + 1, len(new_positions)):
                    if new_positions[j][0] is not None:
                        next_start = new_positions[j][0]
                        break
                if prev_end is not None and next_start is not None:
                    new_positions[i] = (prev_end, next_start)
                elif next_start is None and prev_end is not None:
                    new_positions[i] = (
                        prev_end,
                        len(fastq_data[read_id]["sequence"]) - 1,
                    )
                else:
                    print(new_positions)
                    raise AttributeError("Could not find a valid gene start or end position.")
                assert None not in list(new_positions[i]), new_positions
        return new_positions

    def correct_gene_positions_on_read(
        self,
        first_shared_read_index,
        last_shared_read_index,
        alignment_subset,
        read_id,
        fastq_data,
    ):
        # get the gene positions
        gene_positions = self.get_gene_positions()[read_id][:]
        # get the new gene position prefix
        position_prefix = self.get_gene_position_prefix(gene_positions, first_shared_read_index)
        # get the new gene position suffix
        position_suffix = self.get_gene_position_suffix(gene_positions, last_shared_read_index)
        # get the existing gene position core
        core_gene_positions = self.get_gene_position_core(
            gene_positions, first_shared_read_index, last_shared_read_index
        )
        # get the new gene position core
        new_core_gene_positions = self.get_new_gene_position_core(
            alignment_subset, core_gene_positions
        )
        # join the prefix, core and suffix to get the new gene positions
        new_positions = self.join_gene_position_ends_with_core(
            position_prefix, position_suffix, new_core_gene_positions
        )
        # genes we have had to add have invalid start and stop positions
        corrected_new_positions = self.replace_invalid_gene_positions(
            new_positions, fastq_data, read_id
        )
        # replace the gene positions with the corrected new positions
        self.get_gene_positions()[read_id] = corrected_new_positions
        # double check that the number of genes is equal to the number of gene positions
        assert len(self.get_reads()[read_id]) == len(self.get_gene_positions()[read_id]), (
            str(len(self.get_reads()[read_id])) + "/" + str(len(self.get_gene_positions()[read_id]))
        )
        return self.get_gene_positions()[read_id]

    def modify_alignment_subset(self, alignment_subset, genes_on_read):
        true_path = [c[0] for c in alignment_subset if c[0] != "*"]
        if true_path == genes_on_read:
            return alignment_subset
        return self.needleman_wunsch(true_path, genes_on_read)

    def compare_paths(self, lower_coverage_genes, fw_higher_coverage_genes):
        # get the alignment of this path
        fw_alignment = self.needleman_wunsch(fw_higher_coverage_genes, lower_coverage_genes)
        rv_alignment = self.reverse_gene_alignment(fw_alignment)
        # get the SNP differences
        SNP_count = self.count_snps_in_alignment(fw_alignment)
        # get the INDEL differences
        INDEL_count = self.count_indels_in_alignment(fw_alignment)
        return fw_alignment, rv_alignment, SNP_count, INDEL_count

    def get_minimizers_from_minhashes(self, path, path_minimizers):
        path_minimizer_set = set()
        for minHash in path_minimizers[tuple(path)]:
            path_minimizer_set.update(minHash.hashes)
        return path_minimizer_set

    def define_correction_operations(
        self,
        paths,
        path_coverages,
        reads_to_correct,
        correction_operations,
        path_minimizers,
        seen_nodes,
        threshold,
    ):
        # initialise a set to keep track of the paths we have already corrected
        corrected_paths = set()
        # iterate through the paths
        for p in paths:
            path_coverages.append(p[1])
        for i in range(len(paths)):
            # split the path tuple
            higher_coverage_path, higher_coverage = paths[i]
            # get the nodes in the path
            higher_coverage_path = [n[0] for n in higher_coverage_path]
            higher_coverage_path_set = set(higher_coverage_path)
            higher_coverage_path_tuple = tuple(higher_coverage_path)
            # check if the current high coverage path has been corrected
            if higher_coverage_path_tuple not in corrected_paths:
                if any(n in seen_nodes for n in higher_coverage_path):
                    continue
                # get the genes in the higher coverage path
                high_coverage_minimizers = self.get_minimizers_from_minhashes(
                    higher_coverage_path, path_minimizers
                )
                # iterate through the remaining paths
                for lower_coverage_path, lower_coverage in paths[i + 1 :]:
                    # get the nodes in the path
                    lower_coverage_path = [n[0] for n in lower_coverage_path]
                    lower_coverage_path_tuple = tuple(lower_coverage_path)
                    # skip the correction if we have already corrected this path
                    if lower_coverage_path_tuple in corrected_paths:
                        continue
                    if any(n in seen_nodes for n in lower_coverage_path):
                        continue
                    low_coverage_minimizers = self.get_minimizers_from_minhashes(
                        lower_coverage_path, path_minimizers
                    )
                    # see how similar the paths are based on their minimizers
                    containment = max(
                        [
                            len(high_coverage_minimizers & low_coverage_minimizers)
                            / len(low_coverage_minimizers),
                            len(high_coverage_minimizers & low_coverage_minimizers)
                            / len(high_coverage_minimizers),
                        ]
                    )
                    # if the containment is greater than the threshold
                    correct_path = False
                    if containment > threshold:
                        correct_path = True
                    if correct_path is True:
                        # define the operation
                        operation = (
                            lower_coverage_path_tuple,
                            higher_coverage_path_tuple,
                        )
                        # store the correction operation
                        correction_operations.add(operation)
                        corrected_paths.add(lower_coverage_path_tuple)
                        # store the nodes we haave seen
                        for n in lower_coverage_path:
                            if n not in higher_coverage_path_set:
                                seen_nodes[n] = operation
        return path_coverages

    def get_path_reads_to_correct(self, reads_to_correct, seen_nodes):
        for n in seen_nodes:
            operation = seen_nodes[n]
            for read in self.get_node_by_hash(n).get_reads():
                if read not in reads_to_correct:
                    reads_to_correct[read] = operation

    def correct_bubble_paths(
        self,
        bubbles,
        fastq_data,
        path_minimizers,
        genesOfInterest,
        min_path_coverage,
        threshold=0.80,
    ):
        seen_nodes = {}
        correction_operations = set()
        reads_to_correct = {}
        path_coverages = []
        for pair in bubbles:
            if len(bubbles[pair]) > 1:
                # sort the paths from highest to lower coverage
                paths = sorted(list(bubbles[pair]), key=lambda x: x[1], reverse=True)
                # look through the lower coverage paths in this cluster
                path_coverages = self.define_correction_operations(
                    paths,
                    path_coverages,
                    reads_to_correct,
                    correction_operations,
                    path_minimizers,
                    seen_nodes,
                    threshold,
                )
        # get the reads to correct and the operations
        self.get_path_reads_to_correct(reads_to_correct, seen_nodes)
        # align the paths
        fw_alignments = {}
        bw_alignments = {}
        fw_counters = {}
        bw_counters = {}
        for operation in correction_operations:
            # get the genes in the higher coverage path
            fw_higher_coverage_genes = self.get_genes_in_unitig(list(operation[1]))
            # get the genes in the lower coverage path
            lower_coverage_genes = self.get_genes_in_unitig(list(operation[0]))
            # get the alignment of the high coverage and low coverage path
            fw_alignment, rv_alignment, SNP_count, INDEL_count = self.compare_paths(
                lower_coverage_genes, fw_higher_coverage_genes
            )
            # make sure that we do not delete AMR genes
            if any(
                c[1][1:] in genesOfInterest and c[0][1:] not in genesOfInterest
                for c in fw_alignment
            ):
                continue
            fw_alignments[operation] = fw_alignment
            bw_alignments[operation] = rv_alignment
            # get the k-mers in the low coverage path
            gene_mers = []
            reverse_gene_mers = []
            for i in range(len(lower_coverage_genes) - (self.get_kmerSize() - 1)):
                # take a slice of the list of Genes from index i to i + kmerSize
                gene_mer = lower_coverage_genes[i : i + self.get_kmerSize()]
                gene_mers.append(tuple(gene_mer))
                reverse_gene_mers.append(tuple(self.reverse_list_of_genes(gene_mer)))
            # make counters for the gene-mer tuples
            fw_counter = Counter(gene_mers)
            bw_counter = Counter(reverse_gene_mers)
            fw_counters[operation] = fw_counter
            bw_counters[operation] = bw_counter
        # iterate through the reads
        for read_id in reads_to_correct:
            if reads_to_correct[read_id] not in fw_alignments:
                continue
            # get the forward and backward alignments
            fw_alignment = fw_alignments[reads_to_correct[read_id]]
            rv_alignment = bw_alignments[reads_to_correct[read_id]]
            # get the list of genes on this read
            genes_on_read = self.get_reads()[read_id][:]
            # get the gene-mers on the read
            gene_mers_on_read = self.get_gene_mer_strings(genes_on_read)
            # Determine which replacement dict to use based on shared counts
            read_alignment = self.reorient_alignment(
                gene_mers_on_read,
                fw_counters[reads_to_correct[read_id]],
                bw_counters[reads_to_correct[read_id]],
                fw_alignment,
                rv_alignment,
            )
            if read_alignment is None:
                continue
            # get the mappings of path gene to alignment gene
            higher_mapping, lower_mapping = self.get_path_to_alignment_mapping(read_alignment)
            # get the genes in the lower coverage path
            genes_in_lower_coverage_path = [a[1] for a in read_alignment if not a[1] == "*"]
            # get the longest common sublist
            (
                common_sublist,
                (start_path, end_path),
                (first_shared_read_index, last_shared_read_index),
            ) = self.longest_common_sublist(genes_in_lower_coverage_path, genes_on_read)
            # get the alignment columns for the first and last elements
            alignment_subset = read_alignment[
                lower_mapping[start_path] : lower_mapping[end_path] + 1
            ]
            # modify the alignment subset
            alignment_subset = self.modify_alignment_subset(
                alignment_subset,
                genes_on_read[first_shared_read_index : last_shared_read_index + 1],
            )
            # Correct the read using the alignment subset
            if len(alignment_subset) != 0:
                # correct the genes on the read
                self.correct_genes_on_read(
                    genes_on_read,
                    first_shared_read_index,
                    last_shared_read_index,
                    alignment_subset,
                    read_id,
                )
                # correct the gene positions on the read
                self.correct_gene_positions_on_read(
                    first_shared_read_index,
                    last_shared_read_index,
                    alignment_subset,
                    read_id,
                    fastq_data,
                )
        return path_coverages

    def find_sublist_indices(self, main_list, sublist):
        indices = []
        sublist_length = len(sublist)
        # Loop through the main list from start
        for i in range(len(main_list) - sublist_length + 1):
            # Check if the slice of main_list starting at i matches the sublist
            if main_list[i : i + sublist_length] == sublist:
                # Append the start and stop indices to the results
                indices.append((i, i + sublist_length - 1))
        return indices

    def all_sublist_combinations(self, lst):
        sublists = []
        # Iterate over all possible starting indices
        for start in range(len(lst)):
            # Iterate over all possible ending indices for each start
            for end in range(start + 1, len(lst) + 1):
                sublists.append(lst[start:end])
        return sublists

    def get_path_to_alignment_mapping(self, alignment):
        higher_index = 0
        lower_index = 0
        higher_mapping = {}
        lower_mapping = {}
        for i in range(len(alignment)):
            col = alignment[i]
            if col[1] != "*":
                lower_mapping[lower_index] = i
                lower_index += 1
            if col[0] != "*":
                higher_mapping[higher_index] = i
                higher_index += 1
        return higher_mapping, lower_mapping

    def longest_common_sublist(self, a, b):
        """
        Returns the longest contiguous common sublist between two lists `a` and `b`.
        """
        len_a, len_b = len(a), len(b)
        dp = [[0] * (len_b + 1) for _ in range(len_a + 1)]
        max_len = 0
        end_a = 0
        end_b = 0

        for i in range(1, len_a + 1):
            for j in range(1, len_b + 1):
                if a[i - 1] == b[j - 1]:
                    dp[i][j] = dp[i - 1][j - 1] + 1
                    if dp[i][j] > max_len:
                        max_len = dp[i][j]
                        end_a = i
                        end_b = j

        start_a = end_a - max_len
        start_b = end_b - max_len
        sublist = a[start_a:end_a]
        return sublist, (start_a, end_a - 1), (start_b, end_b - 1)

    def mp_get_all_paths_between_junctions_in_component(
        self, potential_bubble_starts_component, max_distance, cores
    ):
        # store the paths
        unique_paths = set()
        # iterate through the bubble nodes and correct the pairs to process
        pairs_to_process = set()
        for start_hash, start_direction in potential_bubble_starts_component:
            for stop_hash, stop_direction in potential_bubble_starts_component:
                if start_hash == stop_hash:
                    continue
                pair_tuple = tuple(
                    sorted([(start_hash, start_direction), (stop_hash, stop_direction)])
                )
                pairs_to_process.add(pair_tuple)

        def find_paths_between_pairs(batch, max_distance):
            paths_between_these_nodes = []
            for pair in batch:
                start_hash = pair[0][0]
                stop_hash = pair[1][0]
                start_direction = pair[0][1]
                paths_between_nodes = self.new_find_paths_between_nodes(
                    start_hash, stop_hash, max_distance, start_direction
                )
                # remove paths that do not start and stop in the same direction
                valid_paths_between_nodes = [
                    p
                    for p in paths_between_nodes
                    if p[0] == (start_hash, start_direction)
                    and (p[-1][0], self.get_direction_between_two_nodes(p[-2][0], p[-1][0]))
                    == (stop_hash, stop_direction)
                ]
                # check if there is more than one option
                if len(valid_paths_between_nodes) > 1:
                    for p in valid_paths_between_nodes:
                        paths_between_these_nodes.append(
                            tuple(sorted([p, list(reversed([(t[0], t[1] * -1) for t in p]))])[0])
                        )
            return paths_between_these_nodes

        pairs_to_process = list(pairs_to_process)
        batches = [pairs_to_process[i::cores] for i in range(cores)]
        all_paths = Parallel(n_jobs=cores)(
            delayed(find_paths_between_pairs)(batch, max_distance) for batch in tqdm(batches)
        )
        for batch_result in all_paths:
            unique_paths.update(batch_result)
        return unique_paths

    def get_all_paths_between_junctions_in_component(
        self, potential_bubble_starts_component, max_distance, cores=1
    ):
        # store the paths
        unique_paths = set()
        seen_pairs = set()
        # iterate through the bubble nodes
        for start_hash, start_direction in tqdm(potential_bubble_starts_component):
            for stop_hash, stop_direction in potential_bubble_starts_component:
                if start_hash == stop_hash:
                    continue
                # if tuple(sorted([start_hash, stop_hash])) in seen_pairs:
                #    continue
                # seen_pairs.add(tuple(sorted([start_hash, stop_hash])))
                paths_between_nodes = self.new_find_paths_between_nodes(
                    start_hash, stop_hash, max_distance, start_direction
                )
                # remove paths that do not start and stop in the same direction
                valid_paths_between_nodes = [
                    p
                    for p in paths_between_nodes
                    if p[0] == (start_hash, start_direction)
                    and (p[-1][0], self.get_direction_between_two_nodes(p[-2][0], p[-1][0]))
                    == (stop_hash, stop_direction)
                ]
                # check if there is more than one option
                if len(valid_paths_between_nodes) > 1:
                    for p in valid_paths_between_nodes:
                        unique_paths.add(
                            tuple(sorted([p, list(reversed([(t[0], t[1] * -1) for t in p]))])[0])
                        )
                seen_pairs.add(tuple(sorted([start_hash, stop_hash])))
        return list(unique_paths)

    def separate_paths_by_terminal_nodes(self, sorted_filtered_paths):
        paired_sorted_filtered_paths = {}
        for p in sorted_filtered_paths:
            sorted_terminals = tuple(sorted([p[0][0][0], p[0][-1][0]]))
            if sorted_terminals not in paired_sorted_filtered_paths:
                paired_sorted_filtered_paths[sorted_terminals] = []
            paired_sorted_filtered_paths[sorted_terminals].append(p)
        sorted_filtered_paths = {
            key: value
            for key, value in sorted(
                paired_sorted_filtered_paths.items(),
                key=lambda x: max(len(path[0]) for path in x[1]),
                reverse=True,
            )
        }
        return sorted_filtered_paths

    def is_sublist(self, long_list, sub_list):
        """Check if list is a sublist of long_list."""
        assert isinstance(long_list, list) and isinstance(sub_list, list)
        len_sub = len(sub_list)
        return any(
            sub_list == long_list[i : i + len_sub] for i in range(len(long_list) - len_sub + 1)
        )

    def filter_paths_between_bubble_starts(self, unique_paths):
        tree = Tree({i: p for i, p in enumerate(unique_paths)})
        # sort by ascending length
        unique_paths = sorted(list(unique_paths), key=len)
        filtered_paths = []
        for i in range(len(unique_paths)):
            p = unique_paths[i]
            p_list = list(p)
            rev_p_list = list(reversed(p_list))
            res = [path_id for path_id, path in tree.find_all(p_list)]
            rv_res = [path_id for path_id, path in tree.find_all(rev_p_list)]
            if len(res) == 1 and len(rv_res) == 0 and len(p) > 2:
                filtered_paths.append((p, self.calculate_path_coverage(p)))
        return filtered_paths

    def get_minhash_of_nodes(self, batch, node_minhashes, fastq_data):
        for node_hash in batch:
            node = self.get_node_by_hash(node_hash)
            minhash = sourmash.MinHash(n=0, ksize=13, scaled=10)
            for read in node.get_reads():
                indices = [i for i, n in enumerate(self.get_readNodes()[read]) if n == node_hash]
                positions = [self.get_readNodePositions()[read][i] for i in indices]
                entire_read_sequence = fastq_data[read]["sequence"]
                for p in positions:
                    minhash.add_sequence(entire_read_sequence[p[0] : p[1] + 1], force=True)
            node_minhashes[node_hash] = minhash

    def get_minhash_of_path(self, batch, path_minimizers, node_minhashes):
        for path_tuple in batch:
            for node_hash in path_tuple:
                path_minimizers[path_tuple].append(node_minhashes[node_hash])

    def get_minhashes_for_paths(self, sorted_filtered_paths, fastq_data, cores):
        path_minimizers = defaultdict(set)
        node_minhashes = {}

        for path_tuple, path_coverage in sorted_filtered_paths:
            path = [p[0] for p in path_tuple]
            for node_hash in path:
                if node_hash not in node_minhashes:
                    node_minhashes[node_hash] = None

            path_minimizers[tuple(path)] = []

        # Parallel computation of node minhashes
        keys = list(node_minhashes.keys())
        batches = [keys[i::cores] for i in range(cores)]
        Parallel(n_jobs=cores, prefer="threads")(
            delayed(self.get_minhash_of_nodes)(batch, node_minhashes, fastq_data)
            for batch in batches
        )

        # Parallel computation of path minhashes
        keys = list(path_minimizers.keys())
        batches = [keys[i::cores] for i in range(cores)]
        Parallel(n_jobs=cores, prefer="threads")(
            delayed(self.get_minhash_of_path)(batch, path_minimizers, node_minhashes)
            for batch in batches
        )
        # Ensure that all path minimizers are populated
        assert not any(v is None for v in path_minimizers.values())
        return path_minimizers

    def correct_low_coverage_paths(
        self,
        fastq_data,
        genesOfInterest,
        cores,
        min_path_coverage,
        components_to_skip,
        use_minimizers=False,
    ):
        # ensure there are gene positions in the input
        assert self.get_gene_positions()
        # get the potential start nodes
        potential_bubble_starts = self.identify_potential_bubble_starts()
        # define a maximum distance
        max_distance = self.get_kmerSize() * 3
        # iterate through the components
        path_coverages = []
        # iterate through the components in the graph
        for component in self.components():
            message = "\n\tAmira: popping bubbles using 1 CPU for component "
            message += f"{component} / {len(self.components()) - 1}\n"
            sys.stderr.write(message)
            # skip this component if specified
            if component in components_to_skip:
                continue
            # convert the bubbles nodes to a list
            if component not in potential_bubble_starts:
                continue
            potential_bubble_starts_component = potential_bubble_starts[component]
            # get all the paths of length <= max distance between all pairs of junctions
            unique_paths = self.get_all_paths_between_junctions_in_component(
                potential_bubble_starts_component, max_distance, cores
            )
            # filter the paths out if they are subpaths of longer paths
            filtered_paths = self.filter_paths_between_bubble_starts(unique_paths)
            # sort by path length
            sorted_filtered_paths = sorted(filtered_paths, key=lambda x: len(x[0]), reverse=True)
            # get a minhash object for each path
            if use_minimizers:
                path_minimizers = self.get_minhashes_for_paths(
                    sorted_filtered_paths, fastq_data, cores
                )
            else:
                path_minimizers = None
            # bin the paths based on their terminal nodes
            sorted_filtered_paths = self.separate_paths_by_terminal_nodes(sorted_filtered_paths)
            # clean the paths
            path_coverages += self.correct_bubble_paths(
                sorted_filtered_paths,
                fastq_data,
                path_minimizers,
                genesOfInterest,
                min_path_coverage,
            )
        return self.get_reads(), self.get_gene_positions(), path_coverages, min_path_coverage

    def identify_potential_bubble_starts(self):
        potential_bubble_starts = {}
        for node in self.all_nodes():
            fw_edges = node.get_forward_edge_hashes()  # Forward edges
            bw_edges = node.get_backward_edge_hashes()  # Backward edges
            if len(fw_edges) > 1:  # and len(bw_edges) < 2:
                if not node.get_component() in potential_bubble_starts:
                    potential_bubble_starts[node.get_component()] = []
                potential_bubble_starts[node.get_component()].append((node.__hash__(), 1))
            if len(bw_edges) > 1:  # and len(fw_edges) < 2:
                if not node.get_component() in potential_bubble_starts:
                    potential_bubble_starts[node.get_component()] = []
                potential_bubble_starts[node.get_component()].append((node.__hash__(), -1))
        return potential_bubble_starts

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
            sliced = p[: index + 1]
            if not len(sliced) == 0:
                valid_paths.add(tuple(sliced))
        paths_from_start = {}
        for p in valid_paths:
            p = list(p)
            if 2 < len(p):
                terminals = (p[0][0], p[-1][0])
                if terminals not in paths_from_start:
                    paths_from_start[terminals] = []
                path_coverage = self.calculate_path_coverage(p)
                paths_from_start[terminals].append((([n[0] for n in p]), path_coverage))
        return paths_from_start

    def new_find_paths_between_nodes(
        self, start_hash, end_hash, distance, current_direction, path=None, seen_nodes=None
    ):
        if path is None:
            path = []
        if seen_nodes is None:
            seen_nodes = set()

        path.append(
            (start_hash, current_direction)
        )  # Append the current node and direction to the path
        seen_nodes.add(start_hash)  # Add the current node to the seen nodes to avoid revisiting
        # Check if we've reached the end node or if end_hash is None
        if (end_hash and start_hash == end_hash and len(path) - 1 <= distance) or (
            end_hash is None and len(path) - 1 == distance
        ):
            path_copy = path.copy()  # Make a copy of the path to return
            path.pop()  # Remove the last element before returning to ensure correct backtracking
            return [path_copy]  # Return the path copy

        if (
            len(path) - 1 > distance
        ):  # Check if the current path length exceeds the allowed distance
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
                # Make a copy of seen_nodes for each recursive call
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

    def reverse_complement(self, seq):
        reversed_seq = list(reversed(list(seq)))
        replacement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return "".join([replacement[b] for b in reversed_seq])

    def merge_dict(self, dict1, dict2):
        merged_dict = {}
        # Iterate over each dictionary and merge
        for d in [dict1, dict2]:
            for key, value in d.items():
                if key in merged_dict:
                    merged_dict[key].update(value)  # Merge sets if key already exists
                else:
                    merged_dict[key] = value.copy()
        return merged_dict

    def split_into_subpaths(
        self,
        geneOfInterest,
        pathsOfinterest,
        path_coverages,
        path_reads,
        mean_node_coverage=None,
    ):
        allele_count = 1
        gene_clusters = {}
        # get the mean node coverage for all nodes
        if mean_node_coverage is None:
            mean_node_coverage = self.get_mean_node_coverage()
        # track the start and end points of genes, no alleles should overlap in their sequence
        read_tracking = {}
        # iterate through the paths
        for path in pathsOfinterest:
            # keep track of a modified version of the path
            modified_path = list(path)
            # get the genes in the path
            genes_in_path = list(path)
            reverse_genes_in_path = self.reverse_list_of_genes(genes_in_path)
            # make a separate cluster for each allele
            fw_indices_in_path = {}
            rv_indices_in_path = {}
            for g in range(len(genes_in_path)):
                if genes_in_path[g][1:] == geneOfInterest:
                    fw_indices_in_path[g] = f"{geneOfInterest}_{allele_count}"
                    rv_indices_in_path[len(genes_in_path) - g - 1] = (
                        f"{geneOfInterest}_{allele_count}"
                    )
                    # add the allele to the cluster dictionary
                    gene_clusters[f"{geneOfInterest}_{allele_count}"] = []
                    read_tracking[f"{geneOfInterest}_{allele_count}"] = set()
                    # modify the path with the allele name
                    modified_path[g] = f"{genes_in_path[g][0]}{geneOfInterest}_{allele_count}"
                    # increment the allele count
                    allele_count += 1
            # convert the modified path to a tuple
            modified_path = tuple(modified_path)
            # iterate through the reads in this path
            for read_id in self.get_reads():
                # get the genes on the read
                genes_on_read = self.get_reads()[read_id]
                # get the positions of the path on the read
                if self.is_sublist(genes_on_read, genes_in_path):
                    positions_of_path = self.find_sublist_indices(genes_on_read, genes_in_path)
                    indices_in_path = fw_indices_in_path
                elif self.is_sublist(genes_on_read, reverse_genes_in_path):
                    positions_of_path = self.find_sublist_indices(
                        genes_on_read, reverse_genes_in_path
                    )
                    indices_in_path = rv_indices_in_path
                else:
                    continue
                if len(positions_of_path) == 1:
                    # add the path to the path reads
                    if modified_path not in path_reads:
                        path_reads[modified_path] = set()
                    # add the reads to the path reads
                    path_reads[modified_path].add(read_id)
                    # iterate through the path indices
                    for path_start, path_end in positions_of_path:
                        for gene_index in indices_in_path:
                            assert (
                                self.get_reads()[read_id][path_start + gene_index][1:]
                                == geneOfInterest
                            )
                            # get the sequence start and end of the allele on the read
                            sequence_start, sequence_end = self.get_gene_positions()[read_id][
                                path_start + gene_index
                            ]
                            # add the read to the allele
                            gene_clusters[indices_in_path[gene_index]].append(
                                f"{read_id}_{sequence_start}_{sequence_end}"
                            )
                            # add the read to the tracking dictionary
                            read_tracking[indices_in_path[gene_index]].add(
                                f"{read_id}_{sequence_start}_{sequence_end}"
                            )
        sorted_alleles = sorted(
            [a for a in read_tracking], key=lambda x: len(read_tracking[x]), reverse=True
        )
        clusters_to_delete = set()
        for i in range(len(sorted_alleles)):
            a1 = sorted_alleles[i]
            if a1 in clusters_to_delete:
                continue
            for a2 in sorted_alleles[i + 1 :]:
                if a1 == a2:
                    continue
                if len(read_tracking[a1] & read_tracking[a2]) > 0:
                    clusters_to_delete.add(a2)
        for d in clusters_to_delete:
            del gene_clusters[d]
        return gene_clusters, path_reads

    def new_get_minhashes_for_paths(self, pathsOfInterest, fastq_dict):
        path_minhashes = {}
        for path in pathsOfInterest:
            # make a sourmash minhash object for the path
            minhash = sourmash.MinHash(n=0, ksize=9, scaled=1)
            for read_id in pathsOfInterest[path]:
                # split the read identifier
                read = "_".join(read_id.split("_")[:-2])
                start = read_id.split("_")[-2]
                end = read_id.split("_")[-1]
                # get the sequence of the path
                read_sequence = fastq_dict[read]["sequence"][int(start) : int(end) + 1]
                # add the sequence to the minhash object
                minhash.add_sequence(read_sequence, force=True)
            path_minhashes[path] = minhash
        return path_minhashes

    def find(self, parent, item):
        if parent[item] != item:
            parent[item] = self.find(parent, parent[item])  # Path compression
        return parent[item]

    def union(self, parent, rank, set1, set2):
        root1 = self.find(parent, set1)
        root2 = self.find(parent, set2)
        if root1 != root2:
            # Union by rank
            if rank[root1] > rank[root2]:
                parent[root2] = root1
            elif rank[root1] < rank[root2]:
                parent[root1] = root2
            else:
                parent[root2] = root1
                rank[root1] += 1

    def cluster_paths(self, clusters):
        parent = {}
        rank = {}
        # Initialize parent and rank
        for node in clusters:
            parent[node] = node
            rank[node] = 0
            for connected_node in clusters[node]:
                parent[connected_node] = connected_node
                rank[connected_node] = 0
        # Perform unions
        for node in clusters:
            for connected_node in clusters[node]:
                self.union(parent, rank, node, connected_node)
        # Gather clusters
        result = {}
        for node in parent:
            root = self.find(parent, node)
            if root not in result:
                result[root] = set()
            result[root].add(node)
        return result

    def assess_connectivity(self, pathsOfInterest, minhash_for_paths, threshold):
        cluster_pairs = {}
        for i, p1 in enumerate(pathsOfInterest):
            if p1 not in cluster_pairs:
                cluster_pairs[p1] = set()
            for j in range(i + 1, len(pathsOfInterest)):
                p2 = list(pathsOfInterest.keys())[j]
                # Calculate Jaccard containment in both directions and take the maximum
                jaccard_containment = max(
                    minhash_for_paths[p1].contained_by(minhash_for_paths[p2]),
                    minhash_for_paths[p2].contained_by(minhash_for_paths[p1]),
                )
                if jaccard_containment >= threshold:
                    cluster_pairs[p1].add(p2)
                    # Ensure bidirectional association
                    if p2 not in cluster_pairs:
                        cluster_pairs[p2] = set()
                    cluster_pairs[p2].add(p1)
        return cluster_pairs

    def merge_read_clusters(self, merged_paths, pathsOfInterest):
        merged_clusters = {}
        for cluster in merged_paths:
            merged_clusters[cluster] = set()
            for path in merged_paths[cluster]:
                # add the reads to the cluster
                merged_clusters[cluster].update(pathsOfInterest[path])
        return merged_clusters

    def new_merge_clusters(self, pathsOfInterest, fastq_dict):
        # get the minhash objetcs for the paths
        minhash_for_paths = self.new_get_minhashes_for_paths(pathsOfInterest, fastq_dict)
        # make a dictionary that tells us which paths are similar to one another
        cluster_pairs = self.assess_connectivity(pathsOfInterest, minhash_for_paths, 0.85)
        # merged the paths using a Union-Find algorithm
        merged_paths = self.cluster_paths(cluster_pairs)
        # make the key the highest coverage path
        final_merged_paths = {}
        for key in merged_paths:
            paths_in_cluster = merged_paths[key]
            selected = None
            coverage = 0
            for p in paths_in_cluster:
                if len(self.collect_reads_in_path(list(p))) > coverage:
                    selected = p
            final_merged_paths[selected] = paths_in_cluster
        # merge the clusters of reads for each merged path
        merged_clusters = self.merge_read_clusters(merged_paths, pathsOfInterest)
        return merged_clusters

    def calculate_mean_node_coverage(self):
        # get a list of all node coverages
        all_node_coverages = self.get_all_node_coverages()
        # return the mean of the node covewrages
        return statistics.mean(all_node_coverages)

    def make_intersection_matrix(self):
        node_hashes = [h for h in self.get_nodes()]  # Retrieve all nodes from the graph
        num_nodes = len(node_hashes)
        read_sets = [
            set(self.get_node_by_hash(n_hash).get_list_of_reads()) for n_hash in node_hashes
        ]  # List of sets of reads for each node
        # Initialize the matrix with zeros
        matrix = [[0] * num_nodes for _ in range(num_nodes)]
        # Compute the number of intersecting reads for each pair of nodes
        for i in range(num_nodes):
            for j in range(i, num_nodes):  # Only compute for each pair once
                if i == j:
                    # For the same node, we consider all reads as intersecting with themselves
                    matrix[i][j] = len(read_sets[i])
                else:
                    # Intersection count for different nodes
                    intersecting_reads = read_sets[i].intersection(read_sets[j])
                    matrix[i][j] = matrix[j][i] = len(intersecting_reads)
        return matrix, node_hashes  # iterate through the components

    def get_node_with_highest_subthreshold_connections(self, matrix, threshold):
        # Start by assuming no nodes are below the threshold sufficiently
        highest_count = -1
        node_index = None
        for i, row in enumerate(matrix):
            if not np.any(np.isnan(row)):  # Check if the row is already nullified
                # Count the connections below the threshold
                count_below_threshold = np.sum(row < threshold)
                if count_below_threshold > highest_count:
                    highest_count = count_below_threshold
                    node_index = i
        return node_index

    def filter_nodes_by_intersection(self, matrix, node_hashes, threshold=5):
        matrix = np.array(matrix, dtype=float)  # Ensure matrix is of float type to handle NaN
        nodes_to_check = True

        while nodes_to_check:
            lowest = self.get_node_with_highest_subthreshold_connections(matrix, threshold)
            if lowest is not None:
                # Nullify the selected row and column
                matrix[lowest, :] = np.nan
                matrix[:, lowest] = np.nan
            else:
                nodes_to_check = False  # No more nodes to process
        return

    def trim_fringe_nodes(self, number_of_intersecting_reads, intersection_matrix, node_hashes):
        # get a list of nodes to remove
        nodes_to_delete = []
        for i, n_hash in enumerate(node_hashes):
            if all(v < number_of_intersecting_reads for v in intersection_matrix[i]):
                nodes_to_delete.append(self.get_node_by_hash(n_hash))
        # delete the nodes
        for node in nodes_to_delete:
            self.remove_node(node)
        return self

    def get_AMR_anchors(self, AMRNodes):
        # Initialize the node anchor set
        nodeAnchors = set()
        terminals = {}
        # Iterate over all AMR node hashes
        for nodeHash in AMRNodes:
            terminals[nodeHash] = []
            node = self.get_node_by_hash(nodeHash)
            is_anchor = False  # Assume the node is not an anchor until proven otherwise
            singletons = []
            # double check that there are no fw or reverse AMR connections
            forward_neighbors = self.get_forward_neighbors(node)
            # get the backward neighbors for this node
            backward_neighbors = self.get_backward_neighbors(node)
            # get the forward neighbors that contain an AMR gene
            fw_non_self = [n for n in forward_neighbors if n.__hash__() != nodeHash]
            # get the backward neighbors that contain an AMR gene
            bw_non_self = [n for n in forward_neighbors if n.__hash__() != nodeHash]
            if len(fw_non_self) == 0 or len(bw_non_self) == 0:
                nodeAnchors.add(nodeHash)
            # Iterate through all reads associated with the current AMR node
            for r in node.get_reads():
                read_nodes = self.get_readNodes()[r]
                # If the read has only one node, it is isolated
                if len(read_nodes) == 1 and read_nodes[0] == nodeHash:
                    singletons.append(True)
                    terminals[nodeHash].append(True)
                    break
                else:
                    singletons.append(False)
                    # Identify all AMR node positions in the read
                    AMR_indices = [1 if n in AMRNodes else 0 for n in read_nodes]
                    # Check all occurrences of nodeHash in the read
                    for index in [i for i, n in enumerate(read_nodes) if n == nodeHash]:
                        # Check if the nodeHash is isolated in this read
                        if index != 0 and index != len(read_nodes) - 1:
                            if (index != 0 and AMR_indices[index - 1] == 0) or (
                                index != len(read_nodes) - 1 and AMR_indices[index + 1] == 0
                            ):
                                is_anchor = True
                                break
                            terminals[nodeHash].append(False)
                        else:
                            terminals[nodeHash].append(True)
                if is_anchor:
                    nodeAnchors.add(nodeHash)
                    break
            if all(s is True for s in singletons) or all(t is True for t in terminals[nodeHash]):
                # double check that there are no fw or reverse AMR connections
                forward_neighbors = self.get_forward_neighbors(node)
                # get the backward neighbors for this node
                backward_neighbors = self.get_backward_neighbors(node)
                # get the forward neighbors that contain an AMR gene
                forwardAMRNodes = [n for n in forward_neighbors if n.__hash__() in AMRNodes]
                # get the backward neighbors that contain an AMR gene
                backwardAMRNodes = [n for n in backward_neighbors if n.__hash__() in AMRNodes]
                if len(backwardAMRNodes) == 0 or len(forwardAMRNodes) == 0:
                    nodeAnchors.add(nodeHash)
        for nodeHash in terminals:
            if len(terminals[nodeHash]) > 0:
                if terminals[nodeHash].count(True) / len(terminals[nodeHash]) > 0.3:
                    nodeAnchors.add(nodeHash)
        return nodeAnchors

    def get_singleton_paths(self, all_seen_nodes, nodeAnchors, final_paths, final_path_coverages):
        for a in nodeAnchors:
            if a not in all_seen_nodes:
                final_paths[tuple(self.get_genes_in_unitig([a]))] = len(
                    set(self.get_node_by_hash(a).get_list_of_reads())
                )
                final_path_coverages[tuple(self.get_genes_in_unitig([a]))] = [
                    self.get_node_by_hash(a).get_node_coverage()
                ]

    def get_reads_supporting_path(self, path, suffix_tree):
        # shortest_differentiating_path = f_min_up + full_paths[f]["core"] + f_min_down
        reads_with_full_path = set()
        for read_id, path in suffix_tree.find_all(list(path)):
            new_read_id = read_id.replace("_reverse", "")
            reads_with_full_path.add(new_read_id)
        return reads_with_full_path

    def get_all_sublists(self, lst, gene_call_subset, threshold, geneOfInterest, cores):
        """Generate all sublists meeting the threshold criterion."""
        sublists = {}
        args_list = [
            (i, threshold, geneOfInterest, lst, gene_call_subset) for i in range(1, len(lst) + 1)
        ]
        with Pool(processes=cores) as pool:
            results = pool.map(process_combinations_for_i, args_list)
        for i_result in results:
            for sub_list in i_result:
                if sub_list:
                    sublists[sub_list] = i_result[sub_list]
        return sublists

    def get_full_paths(
        self, node_tree, reads, nodeAnchors, threshold, gene_call_subset, geneOfInterest, cores
    ):
        # Store full and partial blocks using the suffix tree
        full_blocks = {}
        # Iterate through the anchors
        for a1 in nodeAnchors:
            # Extract paths and build a subtree
            suffixes = get_suffixes_from_initial_tree(node_tree, a1)
            reversed_suffixes = {}
            for read in suffixes:
                reversed_suffixes[read] = list(reversed(suffixes[read]))
            sub_tree = Tree(reversed_suffixes)
            process_anchors(sub_tree, nodeAnchors, a1, full_blocks, reads, node_tree, threshold)
        # get all of the subpath options
        gene_blocks = {}
        for f in full_blocks:
            genes_in_path = self.get_genes_in_unitig(f)
            all_sublists_with_coverages = self.get_all_sublists(
                genes_in_path, gene_call_subset, threshold, geneOfInterest, cores
            )
            if len(all_sublists_with_coverages) > 0:
                gene_blocks[f] = all_sublists_with_coverages
        # Filter and return the blocks
        filtered_blocks = filter_blocks({f: full_blocks[f] for f in gene_blocks})
        final_paths = {}
        final_path_coverages = {}
        seen_nodes = set()
        for f1 in filtered_blocks:
            seen_nodes.update(f1)
            differentiating_paths = set()
            if f1 not in gene_blocks:
                continue
            for o1 in gene_blocks[f1]:
                if not any(
                    self.is_sublist(self.get_genes_in_unitig(list(f2)), list(o1))
                    or self.is_sublist(
                        self.get_genes_in_unitig(list(f2)), self.reverse_list_of_genes(list(o1))
                    )
                    for f2 in filtered_blocks
                    if f1 != f2
                ):
                    differentiating_paths.add(o1)
            if len(differentiating_paths) > 0:
                selected_path = sorted(
                    list(differentiating_paths),
                    key=lambda x: (
                        x.count(f"+{geneOfInterest}") + x.count(f"-{geneOfInterest}"),
                        gene_blocks[f1][x],
                        len(x),
                    ),
                    reverse=True,
                )[0]
                final_paths[selected_path] = gene_blocks[f1][selected_path]
                final_path_coverages[selected_path] = [
                    self.get_node_by_hash(n).get_node_coverage() for n in list(f1)
                ]
        return final_paths, seen_nodes, final_path_coverages

    def assign_final_alleles_to_components(
        self, finalAllelesOfInterest, clustered_reads, allele_counts, geneOfInterest
    ):
        # iterate through the alleles
        for allele in finalAllelesOfInterest:
            # get the component of this allele
            for read_id in finalAllelesOfInterest[allele]:
                for node_hash in self.get_readNodes()["_".join(read_id.split("_")[:-2])]:
                    component = self.get_node_by_hash(node_hash).get_component()
                    break
                break
            underscore_split = allele.split("_")
            gene_name = "_".join(underscore_split[:-1])
            if gene_name not in allele_counts:
                allele_counts[gene_name] = 1
            if component not in clustered_reads:
                clustered_reads[component] = {}
            if geneOfInterest not in clustered_reads[component]:
                clustered_reads[component][geneOfInterest] = {}
            clustered_reads[component][geneOfInterest][
                f"{gene_name}_{allele_counts[gene_name]}"
            ] = finalAllelesOfInterest[allele]
            # increment the allele count
            allele_counts[gene_name] += 1

    def get_paths_for_gene(
        self,
        node_suffix_tree,
        gene_call_subset,
        nodeHashesOfInterest,
        threshold,
        geneOfInterest,
        cores,
    ):
        nodeAnchors = self.get_AMR_anchors(nodeHashesOfInterest)
        final_paths, seen_nodes, final_path_coverages = self.get_full_paths(
            node_suffix_tree,
            self.get_readNodes(),
            nodeAnchors,
            threshold,
            gene_call_subset,
            geneOfInterest,
            cores,
        )
        self.get_singleton_paths(seen_nodes, nodeAnchors, final_paths, final_path_coverages)
        return final_paths, final_path_coverages

    def collect_component_missed_genes(
        self,
        component_nodeHashesOfInterest,
        clustered_reads,
        allele_counts,
        geneOfInterest,
        path_reads,
    ):
        # iterate through the components
        for component in component_nodeHashesOfInterest:
            nodeHashesOfInterest = component_nodeHashesOfInterest[component]
            # get the reads containing the nodes
            reads = self.collect_reads_in_path(nodeHashesOfInterest)
            # add the component to the clustered reads
            if component not in clustered_reads:
                clustered_reads[component] = {}
            if geneOfInterest not in clustered_reads[component]:
                clustered_reads[component][geneOfInterest] = {}
            # collect all of the reads containing the gene
            if (
                component in clustered_reads
                and len(clustered_reads[component][geneOfInterest]) == 0
            ):
                # add the allele
                if geneOfInterest not in allele_counts:
                    allele_counts[geneOfInterest] = 1
                allele_name = f"{geneOfInterest}_{allele_counts[geneOfInterest]}"
                allele_name_tuple = tuple([f"+{allele_name}"])
                clustered_reads[component][geneOfInterest][allele_name] = []
                # get the reads
                reads = self.collect_reads_in_path(nodeHashesOfInterest)
                # get the positions of all occurrences of the AMR gene
                for read_id in reads:
                    # get the genes on the read
                    genes = self.get_reads()[read_id]
                    # get the index of each occurrence of the gene on the read
                    indices = [i for i, gene in enumerate(genes) if gene[1:] == geneOfInterest]
                    # get the positions of the gene on each read
                    for i in indices:
                        gene_start, gene_end = self.get_gene_positions()[read_id][i]
                        clustered_reads[component][geneOfInterest][allele_name].append(
                            f"{read_id}_{gene_start}_{gene_end}"
                        )
                    # add the read to the path reads
                    if allele_name_tuple not in path_reads:
                        path_reads[allele_name_tuple] = set()
                    path_reads[allele_name_tuple].add(read_id)
                allele_counts[geneOfInterest] += 1

    def assign_reads_to_genes(
        self, listOfGenes, cores, allele_counts={}, mean_node_coverage=None, path_threshold=5
    ):
        # initialise dictionaries to track the alleles
        clustered_reads = {}
        path_reads = {}
        # get the mean node coverage if it is not specified
        if mean_node_coverage is None:
            mean_node_coverage = self.get_mean_node_coverage()
        # iterate through the genes we are interested in
        for geneOfInterest in tqdm(listOfGenes):
            # get the graph nodes containing this gene
            nodesOfInterest = self.get_nodes_containing(geneOfInterest)
            # get the node hashes containing this gene
            nodeHashesOfInterest = [n.__hash__() for n in nodesOfInterest]
            # get the reads containing the gene
            reads_with_gene = self.collect_reads_in_path(nodeHashesOfInterest)
            # construct a suffix tree of the nodes on the reads
            node_suffix_tree = construct_suffix_tree(
                {r: self.get_readNodes()[r] for r in reads_with_gene}
            )
            # build a gene-based suffix tree
            gene_call_subset = {r: self.get_reads()[r] for r in reads_with_gene}
            rc_reads = {}
            for r in gene_call_subset:
                rc_reads[r + "_reverse"] = self.reverse_list_of_genes(gene_call_subset[r])
            gene_call_subset.update(rc_reads)
            # get the paths containing this gene
            pathsOfInterest, pathCoverages = self.get_paths_for_gene(
                node_suffix_tree,
                gene_call_subset,
                nodeHashesOfInterest,
                mean_node_coverage / 20,
                geneOfInterest,
                cores,
            )
            # split the paths into subpaths
            finalAllelesOfInterest, path_reads = self.split_into_subpaths(
                geneOfInterest, pathsOfInterest, pathCoverages, path_reads, mean_node_coverage
            )
            # assign the clusters to components
            self.assign_final_alleles_to_components(
                finalAllelesOfInterest, clustered_reads, allele_counts, geneOfInterest
            )
            # assign node hashes to components
            component_nodeHashesOfInterest = {}
            for n in nodeHashesOfInterest:
                node_component = self.get_node_by_hash(n).get_component()
                if node_component not in component_nodeHashesOfInterest:
                    component_nodeHashesOfInterest[node_component] = set()
                component_nodeHashesOfInterest[node_component].add(n)
            # collect AMR genes in the component but without a cluster
            self.collect_component_missed_genes(
                component_nodeHashesOfInterest,
                clustered_reads,
                allele_counts,
                geneOfInterest,
                path_reads,
            )
        return clustered_reads, path_reads

    def remove_non_AMR_associated_nodes(self, genesOfInterest):
        # get the reads containing AMR genes
        readsOfInterest = set()
        for geneOfInterest in tqdm(genesOfInterest):
            # get the graph nodes containing this gene
            nodesOfInterest = self.get_nodes_containing(geneOfInterest)
            # get the reads containing this gene
            for node in nodesOfInterest:
                for read in node.get_reads():
                    readsOfInterest.add(read)
        # go through all of the nodes to check for AMR reads
        nodes_to_remove = []
        for node_hash in self.get_nodes():
            node = self.get_node_by_hash(node_hash)
            node_reads = set(node.get_list_of_reads())
            if not len(readsOfInterest.intersection(node_reads)) > 0:
                nodes_to_remove.append(node)
        for node in nodes_to_remove:
            self.remove_node(node)
