import json
import os
import shutil
import statistics
import subprocess
import sys
from collections import Counter, defaultdict, deque
from itertools import product

import numpy as np
import pandas as pd
import pysam
import sourmash
from joblib import Parallel, delayed
from tqdm import tqdm

from amira.construct_edge import Edge
from amira.construct_gene import convert_int_strand_to_string
from amira.construct_gene_mer import GeneMer
from amira.construct_node import Node
from amira.construct_read import Read
from amira.path_finding_operations import construct_suffix_tree, get_full_paths

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
        for readId in tqdm(self.get_reads()):
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
                nodes_in_component = set(self.get_nodes_in_component(component))
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

    def is_sublist(self, long_list, sub_list):
        """Check if list is a sublist of long_list."""
        assert isinstance(long_list, list) and isinstance(sub_list, list)
        len_sub = len(sub_list)
        return any(
            sub_list == long_list[i : i + len_sub] for i in range(len(long_list) - len_sub + 1)
        )

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
        for read_id in tqdm(readNodes):
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
        entire_read_sequence = fastq_data[read_id.split("_")[0]]["sequence"]
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
                    if "_" in read_id:
                        read_id = "_".join(read_id.split("_")[:-1])
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
        self, first_shared_read_index, last_shared_read_index, alignment_subset, read_id, fastq_data
    ):
        # get the gene positions
        gene_positions = self.get_gene_positions()[read_id]
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

    def correct_low_coverage_path(
        self,
        lower_coverage_path,
        fw_gene_mer_counter,
        bw_gene_mer_counter,
        fw_alignment,
        rv_alignment,
        fastq_data,
        corrected_reads,
    ):
        # get the reads in the path
        nodes_to_correct = lower_coverage_path[1:-1]
        reads_in_path = self.collect_reads_in_path(nodes_to_correct)
        for read_id in reads_in_path:
            # get the list of genes on this read
            genes_on_read = self.get_reads()[read_id][:]
            # get the gene-mers on the read
            gene_mers_on_read = self.get_gene_mer_strings(genes_on_read)
            # Determine which replacement dict to use based on shared counts
            alignment = self.reorient_alignment(
                gene_mers_on_read,
                fw_gene_mer_counter,
                bw_gene_mer_counter,
                fw_alignment,
                rv_alignment,
            )
            if alignment is None:
                continue
            # get the alignment columns for the first and last elements
            (
                alignment_subset,
                first_shared_read_index,
                last_shared_read_index,
            ) = self.slice_alignment_by_shared_elements(alignment, genes_on_read)
            # Correct the read using the alignment subset
            if len(alignment_subset) != 0:
                # add the read to the corrected reads dictionary
                if read_id not in corrected_reads:
                    corrected_reads[read_id] = set()
                # get the gene-mers in the low coverage path
                gene_mers_on_lcp = self.get_gene_mer_strings(
                    [col[1] for col in alignment_subset if col[1] != "*"]
                )
                # skip the correction if the path we are correcting comes from a previous correction
                if len(corrected_reads[read_id].intersection(gene_mers_on_lcp)) > 0:
                    continue
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
                # get the k-mers in the high coverage path
                gene_mers_on_hcp = self.get_gene_mer_strings(
                    [col[0] for col in alignment_subset if col[0] != "*"]
                )
                # prevent correcting a correct path
                corrected_reads[read_id].update(gene_mers_on_hcp)
                # check that the number of genes and number of gene positions are equal
                assert len(self.get_reads()[read_id]) == len(self.get_gene_positions()[read_id])
            else:
                pass
            assert None not in self.get_reads()[read_id]
        return corrected_reads

    def compare_paths(self, lower_coverage_genes, fw_higher_coverage_genes):
        # get the alignment of this path
        fw_alignment = self.needleman_wunsch(fw_higher_coverage_genes, lower_coverage_genes)
        rv_alignment = self.reverse_gene_alignment(fw_alignment)
        # get the SNP differences
        SNP_count = self.count_snps_in_alignment(fw_alignment)
        # get the INDEL differences
        INDEL_count = self.count_indels_in_alignment(fw_alignment)
        return fw_alignment, rv_alignment, SNP_count, INDEL_count

    def correct_bubble_paths(
        self,
        bubbles,
        fastq_data,
        path_minimizers,
        genesOfInterest,
        min_path_coverage,
        threshold=0.80,
    ):
        corrected_reads = {}
        # iterate through the path terminals
        path_coverages = []
        all_containments = []
        for pair in tqdm(bubbles):
            if len(bubbles[pair]) > 1:
                # sort the paths from highest to lower coverage
                paths = sorted(list(bubbles[pair]), key=lambda x: x[1], reverse=True)
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
                    # check if the current high coverage path has been corrected
                    if not tuple(higher_coverage_path) in corrected_paths:
                        # get the genes in the higher coverage path
                        fw_higher_coverage_genes = self.get_genes_in_unitig(higher_coverage_path)
                        # get the minimzers if this is the final round of bubble popping
                        if path_minimizers is not None:
                            high_coverage_minimizers = path_minimizers[tuple(higher_coverage_path)]
                        # iterate through the remaining paths
                        for lower_coverage_path, lower_coverage in paths[i + 1 :]:
                            # skip the correction if we have already corrected this path
                            if tuple(lower_coverage_path) in corrected_paths:
                                continue
                            # get the nodes in the path
                            lower_coverage_path = [n[0] for n in lower_coverage_path]
                            # see if this is the final round of correction
                            if path_minimizers is None:
                                # if not, only correct paths with a mean coverage < 15
                                # if lower_coverage >= min_path_coverage:
                                #    continue
                                if all(
                                    self.get_node_by_hash(n).get_node_coverage()
                                    >= min_path_coverage
                                    for n in lower_coverage_path
                                ):
                                    continue
                            else:
                                # if so, get the minimizers of the lower coverage path
                                low_coverage_minimizers = path_minimizers[
                                    tuple(lower_coverage_path)
                                ]
                            # get the genes in the lower coverage path
                            lower_coverage_genes = self.get_genes_in_unitig(lower_coverage_path)
                            # get the k-mers in the low coverage path
                            gene_mers = []
                            reverse_gene_mers = []
                            for i in range(len(lower_coverage_genes) - (self.get_kmerSize() - 1)):
                                # take a slice of the list of Genes from index i to i + kmerSize
                                gene_mer = lower_coverage_genes[i : i + self.get_kmerSize()]
                                gene_mers.append(tuple(gene_mer))
                                reverse_gene_mers.append(
                                    tuple(self.reverse_list_of_genes(gene_mer))
                                )
                            # make counters for the gene-mer tuples
                            fw_counter = Counter(gene_mers)
                            bw_counter = Counter(reverse_gene_mers)
                            # get the alignment of the high coverage and low coverage path
                            fw_alignment, rv_alignment, SNP_count, INDEL_count = self.compare_paths(
                                lower_coverage_genes, fw_higher_coverage_genes
                            )
                            # placeholder to decide if we are correcting this low coverage path
                            correct_path = False
                            # see if we are using the minimizers to decide on correction
                            if path_minimizers is None:
                                # correct the lower coverage path
                                if SNP_count <= 2 and INDEL_count <= 2:
                                    correct_path = True
                            else:
                                # see how similar the paths are based on their minimizers
                                containment = max(
                                    [
                                        len(high_coverage_minimizers & low_coverage_minimizers)
                                        / len(low_coverage_minimizers),
                                        len(high_coverage_minimizers & low_coverage_minimizers)
                                        / len(high_coverage_minimizers),
                                    ]
                                )
                                all_containments.append(str(containment))
                                # if the containment is greater than the threshold
                                if containment > threshold:
                                    correct_path = True
                            if correct_path is True:
                                # make sure that we do not delete AMR genes
                                if any(
                                    c[1][1:] in genesOfInterest and c[0][1:] not in genesOfInterest
                                    for c in fw_alignment
                                ):
                                    continue
                                # correct the lower coverage path to the higher coverage path
                                corrected_reads = self.correct_low_coverage_path(
                                    lower_coverage_path,
                                    fw_counter,
                                    bw_counter,
                                    fw_alignment,
                                    rv_alignment,
                                    fastq_data,
                                    corrected_reads,
                                )
        # with open("all_containments.txt", "w") as o:
        #     o.write("\n".join(sorted(all_containments)))
        return path_coverages

    def split_into_contiguous_chunks(
        self, shared_read_indices_with_lcp, shared_read_indices_with_hcp
    ):
        if not shared_read_indices_with_lcp:
            return []
        # Initialize the list of chunks with the first number
        chunks = []
        current_chunk = [shared_read_indices_with_lcp[0]]
        # Iterate through the numbers starting from the second element
        for i in range(1, len(shared_read_indices_with_lcp)):
            # Check if the current number is contiguous with the last number in the current chunk
            if shared_read_indices_with_lcp[i] == shared_read_indices_with_lcp[i - 1] + 1:
                current_chunk.append(shared_read_indices_with_lcp[i])
            else:
                if shared_read_indices_with_lcp[i] - 1 in shared_read_indices_with_hcp:
                    current_chunk += [
                        shared_read_indices_with_lcp[i] - 1,
                        shared_read_indices_with_lcp[i],
                    ]
                else:
                    chunks.append(current_chunk)
                    current_chunk = [shared_read_indices_with_lcp[i]]
        # Add the last chunk to the list
        chunks.append(current_chunk)
        return chunks

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

    def slice_alignment_by_shared_elements(self, alignment, genes_on_read):
        # get the mappings of path gene to alignment gene
        higher_mapping, lower_mapping = self.get_path_to_alignment_mapping(alignment)
        # get the indices of all genes on the read that are shared with the alignment
        genes_in_lower_coverage_path = [a[1] for a in alignment if not a[1] == "*"]
        shared_read_indices_with_lcp = [
            i for i in range(len(genes_on_read)) if genes_on_read[i] in genes_in_lower_coverage_path
        ]
        genes_in_higher_coverage_path = [a[0] for a in alignment if not a[0] == "*"]
        shared_read_indices_with_hcp = [
            i
            for i in range(len(genes_on_read))
            if genes_on_read[i] in genes_in_higher_coverage_path
        ]
        # split the shared gene indices into contiguous chunks
        read_chunks = self.split_into_contiguous_chunks(
            shared_read_indices_with_lcp, shared_read_indices_with_hcp
        )
        if not len(read_chunks) == 0:
            if len(read_chunks) > 1 and all(len(c) == 1 for c in read_chunks):
                return [], None, None
            current_chunk = []
            current_sorted_sublist_positions = []
            for chunk in read_chunks:
                # get all possible sublist combinations
                sublist_combinations = self.all_sublist_combinations(chunk)
                # filter out the sublists that are not found in the path
                valid_sublists = [
                    s
                    for s in sublist_combinations
                    if self.is_sublist(
                        genes_in_lower_coverage_path,
                        [
                            genes_on_read[i]
                            for i in s
                            if genes_on_read[i] in genes_in_lower_coverage_path
                            and genes_on_read[i] not in genes_in_higher_coverage_path
                        ],
                    )
                ]
                # get the locations of the sublists in the path
                sublist_positions = [
                    (
                        s,
                        self.find_sublist_indices(
                            genes_in_lower_coverage_path,
                            [
                                genes_on_read[i]
                                for i in s
                                if genes_on_read[i] in genes_in_lower_coverage_path
                            ],
                        ),
                    )
                    for s in valid_sublists
                ]
                # filter sublists that do not occur in the lower coverage path
                filtered_sublist_positions = [t for t in sublist_positions if len(t[1]) != 0]
                # sort the sublist positions by length
                sorted_sublist_positions = sorted(
                    filtered_sublist_positions,
                    key=lambda x: max([s[1] - s[0] for s in x[1]]),
                    reverse=True,
                )
                # choose the longest sublist match as the true subset
                if len(chunk) > len(current_chunk):
                    current_chunk = chunk
                    current_sorted_sublist_positions = sorted_sublist_positions
            if len(current_sorted_sublist_positions[0][1]) == 1:
                alignment_indices = [
                    lower_mapping[i] for i in list(current_sorted_sublist_positions[0][1][0])
                ]
                return (
                    alignment[alignment_indices[0] : alignment_indices[-1] + 1],
                    current_sorted_sublist_positions[0][0][0],
                    current_sorted_sublist_positions[0][0][-1],
                )
            else:
                return [], None, None
        else:
            return [], None, None

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

    def filter_paths_between_bubble_starts(self, unique_paths):
        # convert to list
        unique_paths = list(unique_paths)
        # Assuming unique_paths is a list of lists or a similar iterable of iterables.
        filtered_paths = []
        for p in unique_paths:
            p_list = list(p)
            # Skip self-comparison and check for sublists in both forward and reversed directions
            if not any(
                self.is_sublist(list(q), list(p)) or self.is_sublist(list(q), list(reversed(p)))
                for q in unique_paths
                if q != p
            ):
                if len(p_list) > 2:
                    filtered_paths.append((p_list, self.calculate_path_coverage(p_list)))
        return filtered_paths

    def get_minhash_of_nodes(self, batch, node_minhashes, fastq_data):
        for node_hash in batch:
            node = self.get_node_by_hash(node_hash)
            minhash = sourmash.MinHash(n=0, ksize=13, scaled=10)
            for read in node.get_reads():
                indices = [i for i, n in enumerate(self.get_readNodes()[read]) if n == node_hash]
                positions = [self.get_readNodePositions()[read][i] for i in indices]
                entire_read_sequence = fastq_data[read.split("_")[0]]["sequence"]
                for p in positions:
                    minhash.add_sequence(entire_read_sequence[p[0] : p[1] + 1], force=True)
            node_minhashes[node_hash] = minhash

    def get_minhash_of_path(self, batch, path_minimizers, node_minhashes):
        for path_tuple in batch:
            for node_hash in path_tuple:
                path_minimizers[path_tuple].update(node_minhashes[node_hash].hashes)

    def get_minhashes_for_paths(self, sorted_filtered_paths, fastq_data, cores):
        path_minimizers = defaultdict(set)
        node_minhashes = {}

        for path_tuple, path_coverage in sorted_filtered_paths:
            path = [p[0] for p in path_tuple]
            for node_hash in path:
                if node_hash not in node_minhashes:
                    node_minhashes[node_hash] = None

            path_minimizers[tuple(path)] = set()

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
        for component in self.components():
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
        mean_node_coverage=None,
    ):
        allele_count = 1
        gene_clusters = {}
        copy_numbers = {}
        # get the mean node coverage for all nodes
        if mean_node_coverage is None:
            mean_node_coverage = self.get_mean_node_coverage()
        # iterate through the paths
        for path in pathsOfinterest:
            # get the genes in the path
            genes_in_path = self.get_genes_in_unitig(list(path))
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
                    allele_count += 1
            # iterate through the reads in this path
            for read_id in self.get_reads():  # pathsOfinterest[path]:
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
                            # estimate the copy number of the allele
                            copy_numbers[indices_in_path[gene_index]] = max(
                                1.0, min(path_coverages[path]) / mean_node_coverage
                            )
        return gene_clusters, copy_numbers

    def new_get_minhashes_for_paths(self, pathsOfInterest, fastq_dict):
        path_minhashes = {}
        for path in pathsOfInterest:
            # make a sourmash minhash object for the path
            minhash = sourmash.MinHash(n=0, ksize=9, scaled=1)
            for read_id in pathsOfInterest[path]:
                # split the read identifier
                read, start, end = read_id.split("_")
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
        return matrix, node_hashes

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
        # Iterate over all AMR node hashes
        for nodeHash in AMRNodes:
            node = self.get_node_by_hash(nodeHash)
            is_anchor = False  # Assume the node is not an anchor until proven otherwise
            singletons = []
            terminals = []
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
                    terminals.append(True)
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
                            terminals.append(False)
                        else:
                            terminals.append(True)
                if is_anchor:
                    nodeAnchors.add(nodeHash)
                    break
            if all(s is True for s in singletons) or all(t is True for t in terminals):
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
        return nodeAnchors

    def get_singleton_paths(self, final_paths, nodeAnchors):
        all_seen_nodes = set()
        for p in final_paths:
            all_seen_nodes.update(p)
        for a in nodeAnchors:
            if a not in all_seen_nodes:
                final_paths[tuple([a])] = set(self.get_node_by_hash(a).get_list_of_reads())

    def get_paths_for_gene(
        self,
        suffix_tree,
        nodeHashesOfInterest,
        threshold,
    ):
        nodeAnchors = self.get_AMR_anchors(nodeHashesOfInterest)
        final_paths = get_full_paths(suffix_tree, self.get_readNodes(), nodeAnchors, threshold)
        self.get_singleton_paths(final_paths, nodeAnchors)
        final_path_coverages = {}
        for path in final_paths:
            final_path_coverages[path] = [
                self.get_node_by_hash(n).get_node_coverage() for n in path
            ]
        return final_paths, final_path_coverages

    def assign_reads_to_genes(
        self, listOfGenes, fastq_dict, allele_counts={}, mean_node_coverage=None, path_threshold=5
    ):
        # initialise dictionaries to track the alleles
        clustered_reads = {}
        cluster_copy_numbers = {}
        # get the mean node coverage if it is not specified
        if mean_node_coverage is None:
            mean_node_coverage = self.get_mean_node_coverage()
        # construct a suffix tree of the nodes on the reads
        suffix_tree = construct_suffix_tree(self.get_readNodes())
        # iterate through the genes we are interested in
        for geneOfInterest in tqdm(listOfGenes):
            # get the graph nodes containing this gene
            nodesOfInterest = self.get_nodes_containing(geneOfInterest)
            # get the node hashes containing this gene
            nodeHashesOfInterest = [n.__hash__() for n in nodesOfInterest]
            # assign node hashes to components
            component_nodeHashesOfInterest = {}
            for n in nodeHashesOfInterest:
                node_component = self.get_node_by_hash(n).get_component()
                if node_component not in component_nodeHashesOfInterest:
                    component_nodeHashesOfInterest[node_component] = []
                component_nodeHashesOfInterest[node_component].append(n)
            # iterate through the components
            for component in component_nodeHashesOfInterest:
                nodeHashesOfInterest = component_nodeHashesOfInterest[component]
                # get the reads containing the nodes
                reads = self.collect_reads_in_path(nodeHashesOfInterest)
                # get the paths containing this gene
                pathsOfInterest, pathCoverages = self.get_paths_for_gene(
                    suffix_tree,
                    nodeHashesOfInterest,
                    mean_node_coverage / 20,
                    # max(5, mean_node_coverage / 20),
                )
                # split the paths into subpaths
                finalAllelesOfInterest, copy_numbers = self.split_into_subpaths(
                    geneOfInterest, pathsOfInterest, pathCoverages, mean_node_coverage
                )
                # add the component to the clustered reads
                if component not in clustered_reads:
                    clustered_reads[component] = {}
                    cluster_copy_numbers[component] = {}
                if geneOfInterest not in clustered_reads[component]:
                    clustered_reads[component][geneOfInterest] = {}
                    cluster_copy_numbers[component][geneOfInterest] = {}
                # iterate through the alleles
                for allele in finalAllelesOfInterest:
                    underscore_split = allele.split("_")
                    gene_name = "_".join(underscore_split[:-1])
                    if gene_name not in allele_counts:
                        allele_counts[gene_name] = 1
                    # add this allele if there are 5 or more reads
                    if len(finalAllelesOfInterest[allele]) >= mean_node_coverage / 20:
                        clustered_reads[component][geneOfInterest][
                            f"{gene_name}_{allele_counts[gene_name]}"
                        ] = finalAllelesOfInterest[allele]
                        # get the copy number estimates
                        cluster_copy_numbers[component][geneOfInterest][
                            f"{gene_name}_{allele_counts[gene_name]}"
                        ] = copy_numbers[allele]
                        # increment the allele count
                        allele_counts[gene_name] += 1
                    else:
                        message = f"Amira: allele {allele} in component {component} "
                        message += "filtered due to an insufficient number of reads "
                        message += f"({len(finalAllelesOfInterest[allele])}).\n"
                        sys.stderr.write(message)
                # collect all of the reads containing the gene
                if (
                    component in clustered_reads
                    and len(clustered_reads[component][geneOfInterest]) == 0
                ):
                    # add the allele
                    if geneOfInterest not in allele_counts:
                        allele_counts[geneOfInterest] = 1
                    allele_name = f"{geneOfInterest}_{allele_counts[geneOfInterest]}"
                    clustered_reads[component][geneOfInterest][allele_name] = []
                    cluster_copy_numbers[component][geneOfInterest][allele_name] = max(
                        1.0,
                        min(
                            [
                                self.get_node_by_hash(h).get_node_coverage()
                                for h in nodeHashesOfInterest
                            ]
                        ),
                    )
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
                    allele_counts[geneOfInterest] += 1
        return clustered_reads, cluster_copy_numbers, allele_counts

    def get_closest_allele(self, bam_file_path, mapping_type, ref_cov_proportion=None):
        valid_references = []
        invalid_references = []
        ref_covered = {}
        ref_matching = {}
        ref_lengths = {}
        unique_reads = set()
        # Open the BAM file
        with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
            for read in bam_file.fetch():
                # check if the read is mapped
                if read.is_unmapped:
                    continue
                # add the read name to a set
                unique_reads.add(read.query_name)
                # get the length of the reference
                total_length = bam_file.get_reference_length(read.reference_name)
                # add the reference to the coverage dictionary
                if read.reference_name not in ref_covered:
                    ref_covered[read.reference_name] = 0
                    ref_matching[read.reference_name] = 0
                    ref_lengths[read.reference_name] = total_length
                # get the proportion of bases matching the reference
                matching_bases = 0
                for op, length in read.cigartuples:
                    if op == 7:
                        matching_bases += length
                if mapping_type == "reads":
                    prop_matching = (matching_bases / total_length) * 100
                    prop_covered = ref_cov_proportion[read.reference_name]
                if mapping_type == "allele":
                    prop_matching = (matching_bases / read.infer_read_length()) * 100
                    # get the proportion of the reference covered by the query
                    prop_covered = read.query_alignment_length / total_length
                # add the stats to the dictionaries
                if prop_matching > ref_matching[read.reference_name]:
                    ref_matching[read.reference_name] = prop_matching
                if prop_covered > ref_covered[read.reference_name]:
                    ref_covered[read.reference_name] = prop_covered
        for ref in ref_matching:
            if ref_covered[ref] >= 0.85:
                valid_references.append(
                    (ref, ref_matching[ref], ref_lengths[ref], ref_covered[ref])
                )
            else:
                invalid_references.append(
                    (ref, ref_matching[ref], ref_lengths[ref], ref_covered[ref])
                )
        valid_references = sorted(valid_references, key=lambda x: (x[1], x[2]), reverse=True)
        if len(valid_references) != 0:
            return True, valid_references, unique_reads
        else:
            invalid_references = sorted(
                invalid_references, key=lambda x: (x[2], x[1]), reverse=True
            )
            return False, invalid_references, unique_reads

    def create_output_dir(self, output_dir, allele_name):
        output_dir = os.path.join(output_dir, allele_name)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        return output_dir

    def write_fasta(self, file_path, sequences):
        with open(file_path, "w") as outFasta:
            outFasta.write("\n".join(sequences))

    def run_subprocess(self, command):
        subprocess.run(command, shell=True, check=True)

    def map_reads(
        self, output_dir, reference_fasta, input_fasta, output_prefix, minimap_opts, minimap2_path
    ):
        sam_file = os.path.join(output_dir, f"{output_prefix}.mapped.sam")
        bam_file = os.path.join(output_dir, f"{output_prefix}.mapped.bam")
        map_command = f"{minimap2_path} {minimap_opts} {reference_fasta} {input_fasta} > {sam_file}"
        self.run_subprocess(map_command)
        sort_and_index_command = (
            f"samtools sort {sam_file} > {bam_file} && samtools index {bam_file}"
        )
        self.run_subprocess(sort_and_index_command)
        return bam_file

    def racon_polish(
        self,
        output_dir,
        racon_path,
        input_fasta,
        sam_file,
        reference_fasta,
        output_fasta,
        window_length,
    ):
        racon_command = f"{racon_path} -t 1 --no-trimming -w {window_length} {input_fasta} "
        racon_command += f"{sam_file} {reference_fasta} > {output_fasta}"
        self.run_subprocess(racon_command)

    def prepare_reference_sequences(self, gene_name, reference_genes, output_dir):
        sequences = [f">{allele}\n{seq}" for allele, seq in reference_genes[gene_name].items()]
        self.write_fasta(os.path.join(output_dir, "01.reference_alleles.fasta"), sequences)
        return [len(seq) for seq in reference_genes[gene_name].values()]

    def get_allele_sequence(self, gene_name, allele, reference_genes):
        return reference_genes[gene_name][allele]

    def racon_one_iteration(
        self,
        output_dir,
        racon_path,
        read_file,
        sam_file,
        sequence_to_polish,
        polished_sequence,
        window_size,
        minimap2_path,
    ):
        # run minimap2
        bam_file = self.map_reads(
            output_dir,
            os.path.join(output_dir, sequence_to_polish),
            read_file,
            sam_file.replace(".mapped.sam", ""),
            "-a --MD -t 1 -x map-ont --eqx",
            minimap2_path,
        )
        # run racon
        self.racon_polish(
            output_dir,
            racon_path,
            read_file,
            os.path.join(output_dir, sam_file),
            os.path.join(output_dir, sequence_to_polish),
            os.path.join(output_dir, polished_sequence),
            window_size,
        )
        # rename the polished sequence
        shutil.copy(
            os.path.join(output_dir, polished_sequence),
            os.path.join(output_dir, sequence_to_polish),
        )
        return bam_file

    def get_ref_allele_pileups(self, bam_file, output_dir):
        read_depths = []
        ref_allele_positions = {}
        cov_proportion = {}
        # Open the BAM file
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            # Iterate over each reference (chromosome/contig)
            for ref in bam.references:
                ref_length = bam.get_reference_length(ref)
                # Initialize a list to store read depths for each position
                depth_list = [0] * ref_length
                # Fetch all reads mapped to this reference
                for read in bam.fetch(ref):
                    if read.is_unmapped:
                        continue
                    # Iterate over the aligned positions of the read
                    for pos in read.get_reference_positions():
                        if pos < ref_length:
                            depth_list[pos] += 1
                # Convert the depth list to a string format
                depth_list_str = ",".join(map(str, depth_list))
                read_depths.append(f">{ref}\n{depth_list_str}")
                # get the first and last non-zero element
                try:
                    first_index = depth_list.index(next(x for x in depth_list if x != 0))
                    last_index = (
                        len(depth_list)
                        - 1
                        - depth_list[::-1].index(next(x for x in reversed(depth_list) if x != 0))
                    )
                except StopIteration:
                    first_index, last_index = None, None
                ref_allele_positions[ref] = (first_index, last_index)
                cov_proportion[ref] = len([d for d in depth_list if d != 0]) / len(depth_list)
        # Write the read depths to the output file
        output_file_path = os.path.join(output_dir, "reference_allele_coverages.txt")
        with open(output_file_path, "w") as output_file:
            output_file.write("\n".join(read_depths))
        return ref_allele_positions, cov_proportion

    def compare_reads_to_references(
        self,
        allele_file,
        output_dir,
        reference_genes,
        racon_path,
        fastqContent,
        phenotypes,
        debug,
        minimap2_path,
    ):
        allele_name = os.path.basename(allele_file).replace(".fastq.gz", "")
        gene_name = "_".join(allele_name.split("_")[:-1])
        output_dir = self.create_output_dir(output_dir, allele_name)
        self.prepare_reference_sequences(gene_name, reference_genes, output_dir)
        bam_file = self.map_reads(
            output_dir,
            os.path.join(output_dir, "01.reference_alleles.fasta"),
            allele_file,
            "02.read",
            "-a --MD -t 1 -x map-ont --eqx",
            minimap2_path,
        )
        ref_allele_positions, ref_cov_proportion = self.get_ref_allele_pileups(bam_file, output_dir)
        if debug is False:
            os.remove(os.path.join(output_dir, "reference_allele_coverages.txt"))
        validity, references, unique_reads = self.get_closest_allele(
            bam_file, "reads", ref_cov_proportion
        )
        if validity is True:
            valid_allele, match_proportion, match_length, coverage_proportion = references[0]
            valid_allele_sequence = self.get_allele_sequence(
                gene_name, valid_allele, reference_genes
            )
            first_base, last_base = ref_allele_positions[valid_allele]
            self.write_fasta(
                os.path.join(output_dir, "03.sequence_to_polish.fasta"),
                [f">{valid_allele}\n{valid_allele_sequence[first_base: last_base+1]}"],
            )
            for _ in range(5):
                try:
                    self.racon_one_iteration(
                        output_dir,
                        racon_path,
                        allele_file,
                        "02.read.mapped.sam",
                        "03.sequence_to_polish.fasta",
                        "04.polished_sequence.fasta",
                        str(len(valid_allele_sequence)),
                        minimap2_path,
                    )
                except subprocess.CalledProcessError:
                    try:
                        gene_name = valid_allele.split(".")[0]
                        closest_ref = valid_allele.split(".")[1]
                    except IndexError:
                        gene_name = "_".join(allele_name.split("_")[:-1])
                        closest_ref = valid_allele
                    return {
                        "Gene name": gene_name,
                        "Sequence name": phenotypes[valid_allele],
                        "Closest reference": closest_ref,
                        "Reference length": match_length,
                        "Identity (%)": 0,
                        "Coverage (%)": 0,
                        "Amira allele": allele_name,
                        "Number of reads": len(unique_reads),
                    }

            bam_file = self.map_reads(
                output_dir,
                os.path.join(output_dir, "01.reference_alleles.fasta"),
                os.path.join(output_dir, "04.polished_sequence.fasta"),
                "05.read",
                "-a --MD -t 1 --eqx",
                minimap2_path,
            )
            validity, references, _ = self.get_closest_allele(bam_file, "allele")
            max_similarity = references[0][1]
            references = [r for r in references if r[1] == max_similarity]
            if len(references) == 1:
                closest_allele, match_proportion, match_length, coverage_proportion = references[0]
                with open(os.path.join(output_dir, "04.polished_sequence.fasta")) as i:
                    final_allele_sequence = "".join(i.read().split("\n")[1:])
                self.write_fasta(
                    os.path.join(output_dir, "06.final_sequence.fasta"),
                    [f">{closest_allele}\n{final_allele_sequence}"],
                )
                try:
                    gene_name = closest_allele.split(".")[0]
                    closest_ref = closest_allele.split(".")[1]
                except IndexError:
                    gene_name = "_".join(allele_name.split("_")[:-1])
                    closest_ref = closest_allele
                return {
                    "Gene name": gene_name,
                    "Sequence name": phenotypes[closest_allele],
                    "Closest reference": closest_ref,
                    "Reference length": match_length,
                    "Identity (%)": round(match_proportion, 1),
                    "Coverage (%)": min(100.0, round(coverage_proportion * 100, 1)),
                    "Amira allele": allele_name,
                    "Number of reads": len(unique_reads),
                }
            if len(references) > 1:
                closest_allele, match_proportion, match_length, coverage_proportion = [], [], [], []
                for r in references:
                    closest_allele.append(r[0])
                    match_proportion.append(r[1])
                    match_length.append(r[2])
                    coverage_proportion.append(r[3])
                with open(os.path.join(output_dir, "04.polished_sequence.fasta")) as i:
                    final_allele_sequence = "".join(i.read().split("\n")[1:])
                self.write_fasta(
                    os.path.join(output_dir, "06.final_sequence.fasta"),
                    [f">{'/'.join(closest_allele)}\n{final_allele_sequence}"],
                )
                try:
                    gene_names = "/".join(
                        sorted(list(set([c.split(".")[0] for c in closest_allele])))
                    )
                    closest_refs = "/".join([c.split(".")[1] for c in closest_allele])
                    phenotypes = "/".join([phenotypes[c] for c in closest_allele])
                except IndexError:
                    gene_names = "_".join(allele_name.split("_")[:-1])
                    closest_refs = "/".join(closest_allele)
                    phenotypes = "/".join([phenotypes[c] for c in closest_allele])
                return {
                    "Gene name": gene_names,
                    "Sequence name": phenotypes,
                    "Closest reference": closest_refs,
                    "Reference length": "/".join([str(m) for m in match_length]),
                    "Identity (%)": "/".join([str(round(p, 1)) for p in match_proportion]),
                    "Coverage (%)": "/".join(
                        [str(min(100.0, round(p * 100, 1))) for p in coverage_proportion]
                    ),
                    "Amira allele": allele_name,
                    "Number of reads": len(unique_reads),
                }
        else:
            if len(references) != 0:
                invalid_allele, match_proportion, match_length, coverage_proportion = references[0]
                try:
                    gene_name = invalid_allele.split(".")[0]
                    closest_ref = invalid_allele.split(".")[1]
                except IndexError:
                    gene_name = "_".join(allele_name.split("_")[:-1])
                    closest_ref = invalid_allele
                return {
                    "Gene name": gene_name,
                    "Sequence name": phenotypes[invalid_allele],
                    "Closest reference": closest_ref,
                    "Reference length": match_length,
                    "Identity (%)": round(match_proportion, 1),
                    "Coverage (%)": min(100.0, round(coverage_proportion * 100, 1)),
                    "Amira allele": allele_name,
                    "Number of reads": len(unique_reads),
                }
            else:
                return {
                    "Gene name": "",
                    "Sequence name": "",
                    "Closest reference": "",
                    "Reference length": 0,
                    "Identity (%)": 0,
                    "Coverage (%)": 0,
                    "Amira allele": allele_name,
                    "Number of reads": len(unique_reads),
                }

    def get_alleles(
        self,
        readFiles,
        racon_path,
        threads,
        output_dir,
        reference_genes,
        fastqContent,
        phenotypes_path,
        debug,
        minimap2_path,
    ):
        # import the phenotypes
        with open(phenotypes_path) as i:
            phenotypes = json.load(i)
        # batch the read files for multi processing
        job_list = [readFiles[i : i + threads] for i in range(0, len(readFiles), threads)]
        gene_data = []
        for subset in tqdm(job_list):
            gene_data += Parallel(n_jobs=threads)(
                delayed(self.compare_reads_to_references)(
                    r,
                    output_dir,
                    reference_genes,
                    racon_path,
                    fastqContent.get("_".join(os.path.basename(r).split("_")[:-1]), None),
                    phenotypes,
                    debug,
                    minimap2_path,
                )
                for r in subset
            )
        gene_data = sorted(gene_data, key=lambda x: x["Gene name"])
        return pd.DataFrame(gene_data)

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
