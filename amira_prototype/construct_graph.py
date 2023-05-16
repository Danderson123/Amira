from itertools import combinations
import os
from tqdm import tqdm

from construct_node import Node
from construct_edge import Edge
from construct_gene_mer import GeneMer
from construct_read import Read
from construct_gene import convert_int_strand_to_string

class GeneMerGraph:

    def __init__(self,
                readDict,
                kmerSize):
        self._reads = readDict
        self._kmerSize = kmerSize
        self._minNodeCoverage = 1
        self._minEdgeCoverage = 1
        self._nodes = {}
        self._edges = {}
        self._readNodes = {}
        # initialise the graph
        for readId in tqdm(self.get_reads()):
            read = Read(readId,
                        self.get_reads()[readId])
            # get the gene mers on this read
            geneMers = [g for g in read.get_geneMers(self.get_kmerSize())]
            # if there are no gene mers there is nothing to add to the graph
            if len(geneMers) == 0:
                continue
            if not len(geneMers) == 1:
                # iterate through all the gene-mers except the last
                for g in range(len(geneMers) - 1):
                    # get the read id for this read
                    readId = read.get_readId()
                    # add a source node to the graph
                    sourceNode = self.add_node(geneMers[g],
                                            [read])
                    # add the sourceNode to the read
                    self.add_node_to_read(sourceNode,
                                    geneMers[g].get_geneMerDirection(),
                                    readId)
                    # increase the source node coverage by 1
                    sourceNode.increment_node_coverage()
                    # add the target node to the graph
                    targetNode = self.add_node(geneMers[g+1],
                                            [read])
                    # add an edge from source to target and target to source
                    sourceToTargetEdge, reverseTargetToSourceEdge = self.add_edge(geneMers[g],
                                                                                geneMers[g+1])
                    # increment the edge coverages by 1
                    sourceToTargetEdge.increment_edge_coverage()
                    reverseTargetToSourceEdge.increment_edge_coverage()
                targetNodeHash = geneMers[-1].__hash__()
                targetNode = self.get_node_by_hash(targetNodeHash)
                # increment the coverage of the target if it is the last gene mer in the read
                targetNode.increment_node_coverage()
                self.add_node_to_read(targetNode,
                                    geneMers[-1].get_geneMerDirection(),
                                    readId)
            else:
                # add a single node to the graph if there is only 1 gene mer
                readId = read.get_readId()
                sourceNode = self.add_node(geneMers[0],
                                        [read])
                self.add_node_to_read(sourceNode,
                                    geneMers[0].get_geneMerDirection(),
                                    readId)
                sourceNode.increment_node_coverage()

    def get_reads(self):
        """ return a dictionary of all reads and their genes """
        return self._reads
    def get_readNodes(self):
        """ return a dictionary of all reads and their node hashes """
        return self._readNodes
    def get_kmerSize(self):
        """ return an integer of the gene-mer size """
        return self._kmerSize
    def get_minEdgeCoverage(self):
        """ return an integer of the minimum coverage for edges """
        return self._minEdgeCoverage
    def get_nodes(self):
        """ return the node dictionary """
        return self._nodes
    def get_edges(self):
        """ return the edge dictionary """
        return self._edges
    def get_minNodeCoverage(self):
        """ return the minimum node coverage """
        return self._minNodeCoverage
    def all_nodes(self):
        """ return a generator for all nodes in the graph and their attributes """
        for nodeHash in self.get_nodes():
            yield self.get_nodes()[nodeHash]
    def add_node_to_read(self,
                        node: Node,
                        nodeDirection: int,
                        readId: str):
        # add the read ID to the read dict if it is not present
        if not readId in self.get_readNodes():
            self.get_readNodes()[readId] = []
        # add the hash for the node as attributes for this read
        self.get_readNodes()[readId].append((node.__hash__(), nodeDirection))
        # each node occurrence will occur in the list (including duplicates)
        return self.get_readNodes()[readId]
    def get_nodes_containing_read(self,
                                readId):
        """ return a list of nodes that contain a read of interest """
        # get the node hashes that contain this read
        listOfNodeHashes = [nodeTuple[0] for nodeTuple in list(self.get_readNodes()[readId])]
        # this function gets the node for a node hash if it has not been filtered from the graph
        listOfNodes = [self.get_node_by_hash(h) for h in listOfNodeHashes if h in self.get_nodes()]
        return listOfNodes
    def add_node_to_nodes(self,
                        node,
                        nodeHash):
        """ add a node to the dictionary of nodes """
        self._nodes[nodeHash] = node
    def get_node_by_hash(self,
                        nodeHash):
        """ return a node object corresponding to a node hash """
        return self.get_nodes()[nodeHash]
    def add_node(self,
                geneMer: GeneMer,
                reads: list) -> Node:
        """
        Add a gene mer to the graph if it does not exist, else get the node.
        Returns the Node itself
        """
        nodeHash = geneMer.__hash__()
        if not nodeHash in self.get_nodes():
            # convert the gene-mer to a node
            node = Node(geneMer)
            # add the node to the graph
            self.add_node_to_nodes(node,
                                nodeHash)
        else:
            node = self.get_node_by_hash(nodeHash)
        # add the reads for this node to the dictionary of reads
        for r in reads:
            node.add_read(r.get_readId())
        return node
    def get_node(self,
                geneMer: GeneMer) -> Node:
        """ get a node given a GeneMer """
        # the gene-mer and node hashes are the same
        nodeHash = geneMer.__hash__()
        # confirm the node is in the node dictionary
        assert nodeHash in self.get_nodes(), "This gene-mer is not in the graph"
        # return the node
        return self.get_nodes()[nodeHash]
    def get_nodes_containing(self,
                            geneOfInterest: str) -> list:
        """
        Return all nodes that contain a given gene.
        This is useful for later, when we want to get nodes that contain AMR genes.
        """
        # ensure there is no strand information present on the requested genes
        assert not (geneOfInterest[0] == "+" or geneOfInterest[0] == "-"), "Strand information cannot be present for any specified genes"
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
    def create_edges(self,
                    sourceNode: Node,
                    targetNode: Node,
                    sourceGeneMerDirection: int,
                    targetGeneMerDirection: int):
        """ create and return the forward and reverse edges from source to target and target to source """
        # define the edge that will go from the source gene mer to the target gene mer
        sourceToTargetEdge = Edge(sourceNode,
                                targetNode,
                                sourceGeneMerDirection,
                                targetGeneMerDirection)
        # define the edge that will go from the reverse source gene mer to the reverse target gene mer
        reverseTargetToSourceEdge = Edge(targetNode,
                                        sourceNode,
                                        targetGeneMerDirection * -1,
                                        sourceGeneMerDirection * -1)
        return sourceToTargetEdge, reverseTargetToSourceEdge
    def get_edge_by_hash(self,
                        edgeHash):
        """ return an edge object corresponding to a edge hash """
        return self.get_edges()[edgeHash]
    def add_edge_to_edges(self,
                        edge):
        """Add the edge to the graph edges if it is not present, else increment the edge coverage by 1.
        Returns the edge.
        """
        # get the hash for the edge we want
        edgeHash = edge.__hash__()
        # see if the edge is in the edge dictionary
        if not edgeHash in self.get_edges():
            self._edges[edgeHash] = edge
        return self.get_edge_by_hash(edgeHash)
    def add_edges_to_graph(self,
                        sourceToTargetEdge,
                        reverseTargetToSourceEdge):
        """ add the edges to the dictionary of edges and return the updated edges """
        sourceToTargetEdge = self.add_edge_to_edges(sourceToTargetEdge)
        reverseTargetToSourceEdge = self.add_edge_to_edges(reverseTargetToSourceEdge)
        return sourceToTargetEdge, reverseTargetToSourceEdge
    def add_edge_to_node(self,
                        node,
                        edge):
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
    def add_edge(self,
                sourceGeneMer: GeneMer,
                targetGeneMer: GeneMer):
        """
        Add a forward and backward edge if they do not exist, else increment the edge coverage by 1.
        Also add the nodes to the graph if they are not present already.
        Returns the updated source and target lists of forward and reverse edges.
        """
        # if the source is not in the graph then add it, else get the node
        sourceNode = self.add_node(sourceGeneMer,
                                [])
        # if the target is not in the graph then add it, else get the node
        targetNode = self.add_node(targetGeneMer,
                                [])
        # create the edge objects
        sourceToTargetEdge, reverseTargetToSourceEdge = self.create_edges(sourceNode,
                                                                        targetNode,
                                                                        sourceGeneMer.get_geneMerDirection(),
                                                                        targetGeneMer.get_geneMerDirection())
        # add the edges to the graph
        sourceToTargetEdge, reverseTargetToSourceEdge = self.add_edges_to_graph(sourceToTargetEdge,
                                                                                reverseTargetToSourceEdge)
        # add the hashes of the edges to the source and target nodes
        sourceNode = self.add_edge_to_node(sourceNode,
                                        sourceToTargetEdge)
        targetNode = self.add_edge_to_node(targetNode,
                                        reverseTargetToSourceEdge)
        return sourceToTargetEdge, reverseTargetToSourceEdge
    def get_degree(self,
                node: Node) -> int:
        """ return an integer of the number of neighbours for this node """
        degrees = len(node.get_forward_edge_hashes()) + len(node.get_backward_edge_hashes())
        return degrees
    def get_forward_neighbors(self,
                            node: Node) -> list:
        """ return a list of nodes corresponding to the forward neighbors for this node """
        return [self.get_edge_by_hash(edgeHash).get_targetNode() for edgeHash in node.get_forward_edge_hashes()]
    def get_backward_neighbors(self,
                            node: Node) -> list:
        """ return a list of nodes corresponding to the backward neighbors for this node """
        return [self.get_edge_by_hash(edgeHash).get_targetNode() for edgeHash in node.get_backward_edge_hashes()]
    def get_all_neighbors(self,
                        node: Node):
        """ return a set of combined forward and reverse node hashes for this node """
        return set([neighbor.__hash__() for neighbor in self.get_forward_neighbors(node) + self.get_backward_neighbors(node)])
    def get_forward_edges(self,
                        node: Node) -> list:
        """ return a list of integers of node identifiers connected to this node by a forward edge """
        return node.get_forward_edge_hashes()
    def get_backward_edges(self,
                        node: Node) -> list:
        """ return a list of integers of node identifiers connected to this node by a backward edge """
        return node.get_backward_edge_hashes()
    def get_edge_hashes_between_nodes(self,
                                    sourceNode: Node,
                                    targetNode: Node) -> tuple:
        """ return a tuple of the source to target and target to source edge hash that are adjacent """
        # check that the two nodes are adjacent
        assert targetNode.__hash__() in self.get_all_neighbors(sourceNode) and sourceNode.__hash__() in self.get_all_neighbors(targetNode)
        # get the edge hash from source to target
        for edgeHash in self.get_forward_edges(sourceNode) + self.get_backward_edges(sourceNode):
            if self.get_edge_by_hash(edgeHash).get_targetNode() == targetNode:
                sourceNodeEdgeHash = edgeHash
                break
        # get the edge hash from target to source
        for edgeHash in self.get_forward_edges(targetNode) + self.get_backward_edges(targetNode):
            if self.get_edge_by_hash(edgeHash).get_targetNode() == sourceNode:
                targetNodeEdgeHash = edgeHash
                break
        assert not (sourceNodeEdgeHash == [] or targetNodeEdgeHash == []), "There are edges missing from the source and target nodes"
        return (sourceNodeEdgeHash, targetNodeEdgeHash)
    def remove_edge_from_edges(self,
                            edgeHash):
        """ remove an edge object from the dictionary of all edges by hash """
        del self._edges[edgeHash]
        assert not edgeHash in self.get_edges(), "This edge was not removed from the graph successfully"
    def remove_edge(self,
                edgeHash: int):
        """ remove an edge from all graph edges and the node edges """
        # check that the edge exists in the graph
        assert edgeHash in self.get_edges(), "This edge is not in the graph"
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
        """ return an integer of the total number of nodes in the filtered graph """
        return len(self.get_nodes())
    def get_total_number_of_edges(self) -> int:
        """ return an integer of the total number of edges in the filtered graph """
        return len(self.get_edges())
    def get_total_number_of_reads(self) -> int:
        """ return an integer of the total number of reads contributing to the filtered graph """
        return len(self.get_reads())
    def remove_node_from_reads(self,
                            node_to_remove):
        """ remove a node from the list of nodes annotated on each read and return the modified node list"""
        for readId in node_to_remove.get_reads():
            self.get_readNodes()[readId] = [nodeTuple for nodeTuple in self.get_readNodes()[readId] if nodeTuple[0] != node_to_remove.__hash__()]
        return self.get_readNodes()[readId]
    def remove_node(self,
                node: Node):
        """ remove a node from the graph and all of its edges """
        # get the node hash
        nodeHash = node.__hash__()
        # confirm the node to remove is in the node dictionary
        assert nodeHash in self.get_nodes(), "This node is not in the graph"
        node_to_remove = self.get_nodes()[nodeHash]
        # remove the nodeHash from the readNodes
        self.remove_node_from_reads(node_to_remove)
        # get the forward edge hashes of the node to remove
        forward_edges_to_remove = node_to_remove.get_forward_edge_hashes()
        # get the backward edge hashes of the node to remove
        backward_edges_to_remove = node_to_remove.get_backward_edge_hashes()
        # remove the forward and backward edges from the edge dictionary
        for edgeHash in forward_edges_to_remove + backward_edges_to_remove:
            targetNode = self.get_edge_by_hash(edgeHash).get_targetNode()
            edgesToTarget = self.get_edge_hashes_between_nodes(node, targetNode)
            for e in edgesToTarget:
                # confirm the edge is in the graph
                self.remove_edge(e)
        # remove the node from the node dictionary
        del self.get_nodes()[nodeHash]
    def set_minNodeCoverage(self,
                            minNodeCoverage):
        """ set the minimum node coverage for the graph and return the new node coverage """
        self._minNodeCoverage = minNodeCoverage
        return self.get_minNodeCoverage()
    def set_minEdgeCoverage(self,
                            minEdgeCoverage):
        """ set the minimum edge coverage for the graph and return the new edge coverage """
        self._minEdgeCoverage = minEdgeCoverage
        return self.get_minEdgeCoverage()
    def filter_graph(self,
                    minNodeCoverage,
                    minEdgeCoverage):
        """ filters the nodes and edges in the graph and then returns the filtered graph """
        # set the new node coverage
        minNodeCoverage = self.set_minNodeCoverage(minNodeCoverage)
        # set the new edge coverage
        minEdgeCoverage = self.set_minEdgeCoverage(minEdgeCoverage)
        # filter the nodes
        nodesToRemove = set()
        for nodeHash in tqdm(self.get_nodes()):
            # mark a node for removal if the coverage is less than the specified coverage
            if not self.get_nodes()[nodeHash].get_node_coverage() > minNodeCoverage - 1:
                nodesToRemove.add(self.get_nodes()[nodeHash])
        # filter the edges
        edgesToRemove = set()
        for edgeHash in tqdm(self.get_edges()):
            # mark an edge for removal if the coverage is less than the specified coverage
            if not self.get_edges()[edgeHash].get_edge_coverage() > minEdgeCoverage - 1:
                edgesToRemove.add(edgeHash)
            # mark an edge for removal if either of the nodes are marked for removal
            if any(n in nodesToRemove for n in [self.get_edges()[edgeHash].get_sourceNode(),
                                                self.get_edges()[edgeHash].get_targetNode()]):
                edgesToRemove.add(edgeHash)
        for e in edgesToRemove:
            self.remove_edge(e)
        for n in nodesToRemove:
            self.remove_node(n)
        return self
    def write_node_entry(self,
                        node_id,
                        node_string,
                        node_coverage,
                        reads,
                        nodeColor):
        """ return a string of a gml node entry """
        if node_coverage > 100:
            node_coverage = 100
        node_entry = "\tnode\t[\n"
        node_entry += "\t\tid\t" + str(node_id) + "\n"
        node_entry += '\t\tlabel\t"' + node_string + '"\n'
        node_entry += "\t\tcoverage\t" + str(node_coverage) + "\n"
        node_entry += '\t\treads\t"' + ",".join(reads) + '"\n'
        if nodeColor:
            node_entry += '\t\tcolor\t"' + str(nodeColor) + '"\n'
        node_entry += "\t]"
        return node_entry
    def write_edge_entry(self,
                        source_node,
                        target_node,
                        edge_coverage):
        """ return a string of a gml edge entry """
        edge_entry = "\tedge\t[\n"
        edge_entry += "\t\tsource\t" + str(source_node) + "\n"
        edge_entry += "\t\ttarget\t" + str(target_node) + "\n"
        edge_entry += "\t\tweight\t" + str(edge_coverage) + "\n"
        edge_entry += "\t]"
        return edge_entry
    def assign_Id_to_nodes(self):
        """ assign an integer ID to each node in the graph """
        currentNodeId = 0
        for node in self.all_nodes():
            assignedId = node.assign_node_Id(currentNodeId)
            assert assignedId == currentNodeId, "This node was assigned an incorrect ID"
            currentNodeId += 1
    def write_gml_to_file(self,
                        output_file,
                        gml_content):
        """ Writes the gml to a file titled <output_file>.gml. Returns nothing."""
        # check if there is a directory to make
        if os.path.dirname(output_file) != "":
            # make the directory if it does not exist
            if not os.path.exists(os.path.dirname(output_file)):
                os.mkdir(os.path.dirname(output_file))
        # write the gml to the output file
        with open(output_file + ".gml", "w") as outGml:
            outGml.write("\n".join(gml_content))
    def get_gene_mer_genes(self,
                        sourceNode: Node) -> list:
        """ return a list of genes in the canonical gene mer and their strands"""
        # get the canonical gene mer for the node
        geneMer = sourceNode.get_canonical_geneMer()
        return [convert_int_strand_to_string(g.get_strand()) + g.get_name() for g in geneMer]
    def get_reverse_gene_mer_genes(self,
                        sourceNode: Node) -> list:
        """ return a list of genes in the canonical gene mer and their strands"""
        # get the reversed gene mer for the node
        geneMer = sourceNode.get_reverse_geneMer()
        return [convert_int_strand_to_string(g.get_strand()) + g.get_name() for g in geneMer]
    def get_gene_mer_label(self,
                        sourceNode: Node) -> str:
        """ returns a string of the genes in the canonical gene mer for a node and their strands """
        # get a list of the gene strand + gene name for the canonical gene mer
        geneMerGenes = self.get_gene_mer_genes(sourceNode)
        # return a string of the gene mer genes and strands
        return "~~~".join(geneMerGenes)
    def get_nodes_with_degree(self,
                            degree: int):
        """ return a list of node objects with the specified degree """
        assert type(degree) == int, "The input degree must be an integer."
        nodesWithDegree = []
        for node in self.all_nodes():
            nodeDegree = self.get_degree(node)
            if nodeDegree == degree:
                nodesWithDegree.append(node)
        return nodesWithDegree
    def insert_between_nodes_on_read(self,
                            readNodes,
                            preInsertion,
                            postInsertion,
                            node_to_insert):
        """ returns a modified list where "*" has been inserted at the sit where a node will be added later """
        if not len(readNodes) == 1:
            newReadNodes = []
            i = 0
            while i < len(readNodes):
                newReadNodes.append(readNodes[i])
                if i + 1 < len(readNodes) and ((readNodes[i] == preInsertion and readNodes[i + 1] == postInsertion) \
                        or (readNodes[i] == postInsertion and readNodes[i + 1] == preInsertion)):
                    newReadNodes.append(node_to_insert)
                i += 1
            # check if either of the nodes are at the ends but the node adjacent is missing
            if (readNodes[0] == preInsertion and not readNodes[1] == postInsertion) \
                or (readNodes[0] == postInsertion and not readNodes[1] == preInsertion):
                newReadNodes.insert(0, node_to_insert)
            if (readNodes[-1] == preInsertion and not readNodes[-2] == postInsertion) \
                or (readNodes[-1] == postInsertion and not readNodes[-2] == preInsertion):
                newReadNodes.append(node_to_insert)
        else:
            newReadNodes = readNodes
        return newReadNodes
    def make_replacement_dict(self,
                        old_path,
                        new_path):
        # insert a star to indicate where the new node will go
        assert not len(new_path) % 2 == 0, "The gene-mer size must be an odd number"
        if len(old_path) == 0:
            old_path.append('*')
        else:
            assert old_path[0] == new_path[0] and old_path[-1] == new_path[-1], [old_path[0], new_path[0], old_path[-1], new_path[-1]]
        # If the list length is even, insert a star at the midpoint
        if len(old_path) % 2 == 0:
            midpoint = len(old_path) // 2
            old_path.insert(midpoint, '*')
        # make a dictionary that will be used to replace nodes later
        replacementDict = {}
        for n in range(len(old_path)):
            replacementDict[old_path[n]] = new_path[n]
        return replacementDict
    def modify_readNode_list(self,
                            readNodes,
                            replacementDict):
        modifiedReadNodes = readNodes[:]
        to_replace = list(replacementDict.keys())
        for i in range(len(to_replace)):
            if to_replace[i] == "*":
                preInsertion = to_replace[i-1]
                postInsertion = to_replace[i+1]
                # get the indices of the node hashes immediately adjacent to the node to the insertion site
                modifiedReadNodes = self.insert_between_nodes_on_read(modifiedReadNodes,
                                                            preInsertion,
                                                            postInsertion,
                                                            "*")
        return modifiedReadNodes
    def replace_nodes_on_read(self,
                            modifiedReadNodes,
                            replacementDict):
        corrected_list = []
        for nodeHash in modifiedReadNodes:
            if not nodeHash in replacementDict:
                corrected_list.append((nodeHash,
                            self.get_node_by_hash(nodeHash).get_geneMer().get_geneMerDirection()))
            else:
                corrected_list.append((replacementDict[nodeHash],
                            self.get_node_by_hash(replacementDict[nodeHash]).get_geneMer().get_geneMerDirection()))
        return corrected_list
    def correct_read_nodes(self,
                        readId,
                        replacementDict):
        """ replace a single node hash in a list with the elements in a list and return the new readNode dictionary """
        # get the nodes on this read
        readNodes = [n[0] for n in self.get_readNodes()[readId]]
        # get the keys of nodes that are adjacent to the insertion
        newReadNodes = self.modify_readNode_list(readNodes,
                                                replacementDict)
        corrected_list = self.replace_nodes_on_read(newReadNodes,
                                                replacementDict)
        # replace the read nodes with the new list
        self.get_readNodes()[readId] = corrected_list
        return self.get_readNodes()[readId]
    def get_new_node_annotations(self,
                                readId):
        """ traverse a read and return a list of strings representing the new gene annotations and their strands """
        readNodes = self.get_readNodes()[readId]
        # store the new annotations in a list
        newAnnotations = []
        # iterate through the new readNodes
        for n in range(len(readNodes)-1):
            sourceNode = self.get_node_by_hash(readNodes[n][0])
            targetNode = self.get_node_by_hash(readNodes[n+1][0])
            # get the edge between the source and target node
            edgeHash = self.get_edge_hashes_between_nodes(sourceNode,
                                                        targetNode)
            edge = self.get_edge_by_hash(edgeHash[0])
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
    def get_new_read_annotations(self):
        """ return a dictionary of post-error cleaning read annotations """
        # get the corrected read node annotations
        readNodeDict = self.get_readNodes()
        # convert the nodes to genes by traversing the nodes from start to finish
        annotatedReads = {}
        for readId in readNodeDict:
            annotatedReads[readId] = self.get_new_node_annotations(readId)
        return annotatedReads
    def remove_short_linear_paths(self,
                                min_length):
        """ remove nodeHash annotations on reads if they belong to a linear path of length < min_length """
        paths_to_remove = []
        for node in self.all_nodes():
            if self.get_degree(node) == 1:
                path = self.get_linear_path_for_node(node)
                if len(path) > 0 and (len(path) < min_length or all(self.get_node_by_hash(n).get_node_coverage() < 2 for n in path)):
                    paths_to_remove.append(path)
        removed = set()
        for path in paths_to_remove:
            for nodeHash in path:
                if not nodeHash in removed:
                    self.remove_node(self.get_node_by_hash(nodeHash))
                    removed.add(nodeHash)
    def group_paths(self,
                paths):
        """ group paths into tuples if they share a start and end node """
        path_dict = {}
        for path in paths:
            if len(path) > 0:
                key = (path[0], path[-1])
                if not key in path_dict:
                    path_dict[key] = []
                path_dict[key].append(path)
        return [tuple(paths) for paths in path_dict.values()]
    def find_paths_between_two_nodes(self,
                                start,
                                nodesOfInterest,
                                k):

        def dfs(path, current_length):
            node = path[-1]
            if node.__hash__() in nodesOfInterest and not node == path[0]:
                if len(path) < k + 3 and len(path) > k:
                    paths.append(path)
                return
            if len(path) > k + 2:
                return

            for neighbor in [self.get_node_by_hash(n) for n in self.get_all_neighbors(node)]:
                if neighbor not in path and neighbor:
                    dfs(path + [neighbor], current_length + 1)
        paths = []
        dfs([start], 0)
        grouped_paths = self.group_paths(paths)
        return grouped_paths
    def pop_bubbles(self):
        # start by getting a set of all node hashes with exactly 3 or 4 neighbors
        threeOrFourNeighbors = set([node.__hash__() for node in self.get_nodes_with_degree(3) + self.get_nodes_with_degree(4)])
        # keep track of the bubbles we have already corrected
        seenBubbles = set()
        # iterate through the nodes with three or 4 neighbors
        for nodeHash in tqdm(threeOrFourNeighbors):
            # get all paths of length k or k-1 that end with another node that has three or 4 neighbors
            grouped_paths = self.find_paths_between_two_nodes(self.get_node_by_hash(nodeHash),
                                                        threeOrFourNeighbors,
                                                        self.get_kmerSize())
            # iterate through the pairs of paths
            for pair_paths in grouped_paths:
                # convert the nodes in the paths to hashes
                pair_paths = [[n.__hash__() for n in p] for p in pair_paths]
                # remove paths where not all nodes have a degree of 2
                pair_paths = [p for p in pair_paths if all(self.get_degree(self.get_node_by_hash(n)) == 2 for n in p[1:-1])]
                # remove bubbles where the gene annotation tool has missed a gene
                if len(pair_paths) == 2 and not len(pair_paths[0]) == len(pair_paths[1]):
                    # when a gene has been missed, one path will be shorter than the other by 1 node
                    nodes_to_replace, replacement_list = min(pair_paths, key=len), max(pair_paths, key=len)
                    # skip this correction if we have already corrected it
                    bubbleTerminals = tuple(sorted([nodes_to_replace[0], nodes_to_replace[-1]]))
                    if not bubbleTerminals in seenBubbles:
                        # make a dictionary to decide how we are modifying the nodes on the read we are correcting
                        replacementDict = self.make_replacement_dict(nodes_to_replace[:],
                                                            replacement_list[:])
                        # get the reads we are going to correct
                        pathReads = set()
                        for nodeHash in nodes_to_replace:
                            for readId in self.get_node_by_hash(nodeHash).get_reads():
                                pathReads.add(readId)
                        # correct the reads
                        for readId in list(pathReads):
                            old_readNodes = self.get_readNodes()[readId][:]
                            newReadNodes = self.correct_read_nodes(readId,
                                                                replacementDict)
                        # keep track of this pair so that we don't try to correct it again
                        seenBubbles.add(tuple(sorted([nodes_to_replace[0], nodes_to_replace[-1]])))
    def correct_errors(self,
                    min_linearPathLength):
        """ return a dictionary of corrected read annotation to build a new, cleaner gene-mer graph """
        # clean up short linear paths that look like hairs on the gene-mer graph
        self.remove_short_linear_paths(min_linearPathLength)
        # pop bubbles in the graph where a single gene annotation has been missed
        self.pop_bubbles()
        # return a dictionary of new read annotations
        return self.get_new_read_annotations()
    def get_forward_node_from_node(self,
                            sourceNode) -> list:
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
                # if the degree of the target node is 1 or 2 we can extend the linear path to the next node
                if targetNodeDegree == 2 or targetNodeDegree == 1:
                    return True, targetNode, targetNodeDirection
                else:
                    return False, targetNode, targetNodeDirection
        return False, None, None
    def get_backward_node_from_node(self,
                                sourceNode) -> list:
        # get the list of forward edge hashes for this node
        nodeBackwardEdges = sourceNode.get_backward_edge_hashes()
        if len(nodeBackwardEdges) == 1:
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
                # if the degree of the target node is 1 or 2 we can extend the linear path to the next node
                if targetNodeDegree == 2 or targetNodeDegree == 1:
                    # get the direction we are going into the target node
                    return True, targetNode, targetNodeDirection
                else:
                    # else we cannot extend the linear path
                    return False, targetNode, targetNodeDirection
        return False, None, None
    def get_forward_path_from_node(self,
                                node: Node,
                                wantBranchedNode = False) -> list:
        """ return a list of node hashes in the forward direction from this node """
        forward_nodes_from_node = [node.__hash__()]
        # get the next node in the forward direction
        forwardExtend, forwardNode, forwardNodeDirection = self.get_forward_node_from_node(node)
        # if we are extending further in the forward direction, get the next canonical gene mer
        while forwardExtend:
            # if we enter the next node in the forward direction, we get the next forward node
            if forwardNodeDirection == 1:
                forward_nodes_from_node.append(forwardNode.__hash__())
                forwardExtend, forwardNode, forwardNodeDirection = self.get_forward_node_from_node(forwardNode)
            # if we enter the next node in the backward direction, we get the next backward node
            else:
                forward_nodes_from_node.append(forwardNode.__hash__())
                forwardExtend, forwardNode, forwardNodeDirection = self.get_backward_node_from_node(forwardNode)
        # if we want the final branched node too, then we add this to the path
        if wantBranchedNode:
            forward_nodes_from_node.append(forwardNode.__hash__())
        return forward_nodes_from_node
    def get_backward_path_from_node(self,
                                node: Node,
                                wantBranchedNode = False) -> list:
        """ return a list of node hashes in the backward direction from this node """
        backward_nodes_from_node = [node.__hash__()]
        # get the next node in the backward direction
        backwardExtend, backwardNode, backwardNodeDirection = self.get_backward_node_from_node(node)
        # if we are extending further in the backward direction, get the next canonical gene mer
        while backwardExtend:
            if backwardNodeDirection == -1:
                backward_nodes_from_node.insert(0, backwardNode.__hash__())
                backwardExtend, backwardNode, backwardNodeDirection = self.get_backward_node_from_node(backwardNode)
            # if we enter the next node in the forward direction, we get the next forward node
            else:
                backward_nodes_from_node.insert(0, backwardNode.__hash__())
                backwardExtend, backwardNode, backwardNodeDirection = self.get_forward_node_from_node(backwardNode)
        # if we want the final branched node too, then we add this to the path
        if wantBranchedNode:
            backward_nodes_from_node.insert(0, backwardNode.__hash__())
        return backward_nodes_from_node
    def get_linear_path_for_node(self,
                                node: Node) -> list:
        """ return a list of nodes that correspond to the linear path that contains the specified node and does not include the terminal nodes with a degree of more than 2"""
        linear_path = self.get_backward_path_from_node(node)[:-1] + [node.__hash__()] + self.get_forward_path_from_node(node)[1:]
        return linear_path
    def generate_gml(self,
                    output_file: str,
                    geneMerSize: int,
                    min_node_coverage: int,
                    min_edge_coverage: int):
        """ Write a gml of the filtered graph to the output directory. Returns the written content as a list """
        graph_data = ["graph\t["]
        self.assign_Id_to_nodes()
        # iterate through the nodes in the graph
        for sourceNode in self.all_nodes():
            # add the node entry
            nodeEntry = self.write_node_entry(sourceNode.get_node_Id(),
                                            self.get_gene_mer_label(sourceNode),
                                            sourceNode.get_node_coverage(),
                                            [read for read in sourceNode.get_reads()],
                                            sourceNode.get_color())
            graph_data.append(nodeEntry)
            # get the forward edges for this node
            nodeForwardEdgeHashes = sourceNode.get_forward_edge_hashes()
            nodeBackwardEdgeHashes = sourceNode.get_backward_edge_hashes()
            # get the edge objects corresponding to the forward edge hashes
            nodeEdges = [self.get_edge_by_hash(edgeHash) for edgeHash in nodeForwardEdgeHashes + nodeBackwardEdgeHashes]
            for edge in nodeEdges:
                targetNode = edge.get_targetNode()
                edgeEntry = self.write_edge_entry(sourceNode.get_node_Id(),
                                                targetNode.get_node_Id(),
                                                edge.get_edge_coverage())
                graph_data.append(edgeEntry)
        graph_data.append("]")
        output_file = ".".join([output_file,
                                str(geneMerSize),
                                str(min_node_coverage),
                                str(min_edge_coverage)])
        self.write_gml_to_file(output_file,
                            graph_data)
        return graph_data
