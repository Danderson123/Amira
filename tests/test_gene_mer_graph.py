import unittest
import sys
sys.path.insert(0, "amira_prototype")

from construct_graph import GeneMerGraph
from construct_read import Read
from construct_node import Node

class TestGeneMerConstructor(unittest.TestCase):

    def test___init_GeneMerGraph(self):
        # setup
        graph = GeneMerGraph({},
                            0,
                            1,
                            2)
        # execution
        actual_reads = graph.get_reads()
        actual_kmerSize = graph.get_kmerSize()
        actual_minGeneMerCoverage = graph.get_minGeneMerCoverage()
        actual_minEdgeCoverage = graph.get_minEdgeCoverage()
        actual_nodes = graph.get_nodes()
        actual_edges = graph.get_edges()
        # assertion
        expected_reads = {}
        expected_kmerSize = 0
        expected_minGeneMerCoverage = 1
        expected_minEdgeCoverage = 2
        expected_nodes = {}
        expected_edges = {}
        self.assertEqual(actual_reads, expected_reads)
        self.assertEqual(actual_kmerSize, expected_kmerSize)
        self.assertEqual(actual_minGeneMerCoverage, expected_minGeneMerCoverage)
        self.assertEqual(actual_minEdgeCoverage, expected_minEdgeCoverage)
        self.assertEqual(actual_nodes, expected_nodes)
        self.assertEqual(actual_edges, expected_edges)

    def test___all_nodes(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "-gene3", "+gene2", "-gene1"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        for g in geneMers:
            graph.add_node(g)
        # execution
        actual_all_nodes = [x for x in graph.all_nodes()]
        actual_number_of_nodes = len(actual_all_nodes)
        # assertion
        expected_number_of_nodes = 6
        self.assertEqual(actual_number_of_nodes, expected_number_of_nodes)

    def test___add_node_to_empty_graph(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes)
        geneMer = [x for x in read1.get_geneMers(3)][0]
        graph = GeneMerGraph({},
                        3,
                        1,
                        1)
        # execution
        actual_returned_node = graph.add_node(geneMer)
        actual_returned_node.add_read("read1")
        actual_node_read_list = [x for x in actual_returned_node.get_reads()]
        actual_node_hash = actual_returned_node.__hash__()
        actual_node_coverage = actual_returned_node.get_node_coverage()
        # assertion
        expected_node_read_list = ["read1"]
        expected_node_hash = geneMer.__hash__()
        expected_node_coverage = 1
        self.assertEqual(actual_node_read_list, expected_node_read_list)
        self.assertEqual(actual_node_hash, expected_node_hash)
        self.assertEqual(actual_node_coverage, expected_node_coverage)

    def test___add_node_to_non_empty_graph(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes1)
        geneMer1 = [x for x in read1.get_geneMers(3)][0]
        graph = GeneMerGraph([],
                        3,
                        1,
                        1)
        graph.add_node(geneMer1)
        genes2 = ["+gene4", "-gene3", "+gene2"]
        read2 = Read("read2",
                    genes2)
        geneMer2 = [x for x in read2.get_geneMers(3)][0]
        # execution
        actual_returned_node2 = graph.add_node(geneMer2)
        actual_returned_node2.add_read("read2")
        actual_node2_read_list = [x for x in actual_returned_node2.get_reads()]
        actual_node2_hash = actual_returned_node2.__hash__()
        actual_node2_coverage = actual_returned_node2.get_node_coverage()
        # assertion
        expected_node2_read_list = ["read2"]
        expected_node2_hash = geneMer2.__hash__()
        expected_node2_coverage = 1
        self.assertEqual(len(graph.get_nodes()), 2)
        self.assertEqual(actual_node2_read_list, expected_node2_read_list)
        self.assertEqual(actual_node2_hash, expected_node2_hash)
        self.assertEqual(actual_node2_coverage, expected_node2_coverage)

    def test___add_same_node_to_non_empty_graph(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes)
        geneMer = [x for x in read1.get_geneMers(3)][0]
        graph = GeneMerGraph([],
                        3,
                        1,
                        1)
        graph.add_node(geneMer)
        # execution
        actual_returned_node = graph.add_node(geneMer)
        actual_returned_node.add_read("read1")
        actual_returned_node.add_read("read2")
        actual_returned_node.increment_node_coverage()
        actual_returned_node = graph.get_node(geneMer)
        actual_node_read_list = [x for x in actual_returned_node.get_reads()]
        actual_node_hash = actual_returned_node.__hash__()
        actual_node_coverage = actual_returned_node.get_node_coverage()
        # assertion
        expected_node_read_list = ["read1", "read2"]
        expected_node_hash = geneMer.__hash__()
        expected_node_coverage = 2
        self.assertEqual(len(graph.get_nodes()), 1)
        self.assertEqual(actual_node_read_list, expected_node_read_list)
        self.assertEqual(actual_node_hash, expected_node_hash)
        self.assertEqual(actual_node_coverage, expected_node_coverage)

    def test___get_node_in_graph(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes1)
        geneMer1 = [x for x in read1.get_geneMers(3)][0]
        graph = GeneMerGraph([],
                        3,
                        1,
                        1)
        graph.add_node(geneMer1)
        genes2 = ["+gene4", "-gene3", "+gene2"]
        read2 = Read("read2",
                    genes2)
        geneMer2 = [x for x in read2.get_geneMers(3)][0]
        graph.add_node(geneMer2)
        # execution
        actual_node1 = graph.get_node(geneMer1)
        actual_node1.add_read("read1")
        actual_node1 = graph.get_node(geneMer1)
        actual_node1_read_list = [x for x in actual_node1.get_reads()]
        actual_node1_hash = actual_node1.__hash__()
        actual_node1_coverage = actual_node1.get_node_coverage()
        actual_node2 = graph.get_node(geneMer2)
        actual_node2.add_read("read2")
        actual_node2 = graph.get_node(geneMer2)
        actual_node2_read_list = [x for x in actual_node2.get_reads()]
        actual_node2_hash = actual_node2.__hash__()
        actual_node2_coverage = actual_node2.get_node_coverage()
        # assertion
        expected_node1_read_list = ["read1"]
        expected_node1_hash = geneMer1.__hash__()
        expected_node1_coverage = 1
        expected_node2_read_list = ["read2"]
        expected_node2_hash = geneMer2.__hash__()
        expected_node2_coverage = 1
        self.assertEqual(actual_node1_read_list, expected_node1_read_list)
        self.assertEqual(actual_node1_hash, expected_node1_hash)
        self.assertEqual(actual_node1_coverage, expected_node1_coverage)
        self.assertEqual(actual_node2_read_list, expected_node2_read_list)
        self.assertEqual(actual_node2_hash, expected_node2_hash)
        self.assertEqual(actual_node2_coverage, expected_node2_coverage)

    def test___get_node_not_in_graph(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes1)
        geneMer1 = [x for x in read1.get_geneMers(3)][0]
        graph = GeneMerGraph([],
                        3,
                        1,
                        1)
        graph.add_node(geneMer1)
        genes2 = ["+gene4", "-gene3", "+gene2"]
        read2 = Read("read2",
                    genes2)
        geneMer2 = [x for x in read2.get_geneMers(3)][0]
        # assertion
        self.assertRaises(AssertionError, graph.get_node, geneMer2)

    def test___get_nodes_containing_subset(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "-gene3", "+gene2", "-gene1"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        for g in geneMers:
            graph.add_node(g)
        selectedGenes = ["gene2", "gene6"]
        # execution
        actual_selectedNodes = [n for n in graph.get_nodes_containing(selectedGenes)]
        actual_selectedGenes = [[g.get_name() for g in n.get_canonical_geneMer()] for n in actual_selectedNodes]
        actual_geneMerCount = len(actual_selectedGenes)
        # assertion
        expected_geneMerCount = 5
        self.assertTrue(all(any(g in ng for g in selectedGenes) for ng in actual_selectedGenes))
        self.assertEqual(actual_geneMerCount, expected_geneMerCount)

    def test___get_nodes_containing_all(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "-gene3", "+gene2", "-gene1"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        for g in geneMers:
            graph.add_node(g)
        selectedGenes = ["gene1", "gene2", "gene3", "gene4", "gene5", "gene6", "gene3", "gene2", "gene1"]
        # execution
        actual_selectedNodes = [n for n in graph.get_nodes_containing(selectedGenes)]
        actual_selectedGenes = [[g.get_name() for g in n.get_canonical_geneMer()] for n in actual_selectedNodes]
        actual_geneMerCount = len(actual_selectedGenes)
        # assertion
        expected_geneMerCount = 6
        self.assertTrue(all(any(g in ng for g in selectedGenes) for ng in actual_selectedGenes))
        self.assertEqual(actual_geneMerCount, expected_geneMerCount)

    def test___get_nodes_containing_gene_not_in_graph(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "-gene3", "+gene2", "-gene1"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        for g in geneMers:
            graph.add_node(g)
        selectedGenes = ["gene10"]
        # execution
        actual_selectedNodes = [n for n in graph.get_nodes_containing(selectedGenes)]
        actual_selectedGenes = [[g.get_name() for g in n.get_canonical_geneMer()] for n in actual_selectedNodes]
        # assertion
        expected_selectedGenes = []
        self.assertEqual(actual_selectedGenes, expected_selectedGenes)

    def test___get_nodes_containing_none(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "-gene3", "+gene2", "-gene1"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        for g in geneMers:
            graph.add_node(g)
        selectedGenes = []
        # execution
        actual_selectedNodes = [n for n in graph.get_nodes_containing(selectedGenes)]
        actual_selectedGenes = [[g.get_name() for g in n.get_canonical_geneMer()] for n in actual_selectedNodes]
        # assertion
        expected_selectedGenes = []
        self.assertEqual(actual_selectedGenes, expected_selectedGenes)

    def test___check_no_strand_in_genes_no_strands(self):
        # setup
        listOfGenes = ["gene1", "gene2", "gene3"]
        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        # execution
        actual_bool = graph.check_no_strand_in_genes(listOfGenes)
        # assertion
        expected_bool = True
        self.assertEqual(actual_bool, expected_bool)

    def test___check_no_strand_in_genes_all_strands(self):
        # setup
        listOfGenes = ["+gene1", "-gene2", "+gene3"]
        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        # execution
        actual_bool = graph.check_no_strand_in_genes(listOfGenes)
        # assertion
        expected_bool = False
        self.assertEqual(actual_bool, expected_bool)

    def test___check_no_strand_in_genes_one_strand(self):
        # setup
        listOfGenes = ["gene1", "-gene2", "gene3"]
        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        # execution
        actual_bool = graph.check_no_strand_in_genes(listOfGenes)
        # assertion
        expected_bool = False
        self.assertEqual(actual_bool, expected_bool)

    def test___get_nodes_containing_strand(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "-gene3", "+gene2", "-gene1"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        for g in geneMers:
            graph.add_node(g)
        selectedGenes = ["gene2", "+gene6"]
        # assertion
        self.assertRaises(AssertionError, graph.get_nodes_containing, selectedGenes)

    def test___create_edges_positive_to_positive(self):

        class fakeNode:
            def __init__(self,
                        direction):
                self.direction = direction
            def get_direction(self):
                return self.direction

        # setup
        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        mock_sourceNode = fakeNode(1)
        mock_targetNode = fakeNode(1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        # execution
        actual_sourceToTargetEdge, actual_reverseTargetToSourceEdge = graph.create_edges(mock_sourceNode,
                                                                                        mock_targetNode,
                                                                                        mock_sourceDirection,
                                                                                        mock_targetDirection)
        actual_sourceToTargetEdge_sourceNode = actual_sourceToTargetEdge.get_sourceNode()
        actual_sourceToTargetEdge_targetNode = actual_sourceToTargetEdge.get_targetNode()
        actual_sourceToTargetEdge_sourceDirection = actual_sourceToTargetEdge.get_sourceNodeDirection()
        actual_sourceToTargetEdge_targetDirection = actual_sourceToTargetEdge.get_targetNodeDirection()
        actual_sourceToTargetEdge_coverage = actual_sourceToTargetEdge.get_edge_coverage()
        expected_sourceToTargetEdge_sourceNode = mock_sourceNode
        expected_sourceToTargetEdge_targetNode = mock_targetNode
        expected_sourceToTargetEdge_sourceDirection = mock_sourceDirection
        expected_sourceToTargetEdge_targetDirection = mock_targetDirection
        expected_sourceToTargetEdge_coverage = 0

        actual_reverseTargetToSourceEdge_sourceNode = actual_reverseTargetToSourceEdge.get_sourceNode()
        actual_reverseTargetToSourceEdge_targetNode = actual_reverseTargetToSourceEdge.get_targetNode()
        actual_reverseTargetToSourceEdge_sourceDirection = actual_reverseTargetToSourceEdge.get_sourceNodeDirection()
        actual_reverseTargetToSourceEdge_targetDirection = actual_reverseTargetToSourceEdge.get_targetNodeDirection()
        actual_reverseTargetToSourceEdge_coverage = actual_reverseTargetToSourceEdge.get_edge_coverage()
        expected_reverseTargetToSourceEdge_sourceNode = mock_targetNode
        expected_reverseTargetToSourceEdge_targetNode = mock_sourceNode
        expected_reverseTargetToSourceEdge_sourceDirection = mock_targetDirection * -1
        expected_reverseTargetToSourceEdge_targetDirection = mock_sourceDirection * -1
        expected_reverseTargetToSourceEdge_coverage = 0

        actual_sourceToTargetEdge_hash = actual_sourceToTargetEdge.__hash__()
        actual_reverseTargetToSourceEdge_hash = actual_reverseTargetToSourceEdge.__hash__()
        # assertion
        self.assertEqual(actual_sourceToTargetEdge_sourceNode, expected_sourceToTargetEdge_sourceNode)
        self.assertEqual(actual_sourceToTargetEdge_targetNode, expected_sourceToTargetEdge_targetNode)
        self.assertEqual(actual_sourceToTargetEdge_sourceDirection, expected_sourceToTargetEdge_sourceDirection)
        self.assertEqual(actual_sourceToTargetEdge_targetDirection, expected_sourceToTargetEdge_targetDirection)
        self.assertEqual(actual_sourceToTargetEdge_coverage, expected_sourceToTargetEdge_coverage)

        self.assertEqual(actual_reverseTargetToSourceEdge_sourceNode, expected_reverseTargetToSourceEdge_sourceNode)
        self.assertEqual(actual_reverseTargetToSourceEdge_targetNode, expected_reverseTargetToSourceEdge_targetNode)
        self.assertEqual(actual_reverseTargetToSourceEdge_sourceDirection, expected_reverseTargetToSourceEdge_sourceDirection)
        self.assertEqual(actual_reverseTargetToSourceEdge_targetDirection, expected_reverseTargetToSourceEdge_targetDirection)
        self.assertEqual(actual_reverseTargetToSourceEdge_coverage, expected_reverseTargetToSourceEdge_coverage)

        self.assertNotEqual(actual_sourceToTargetEdge_hash, actual_reverseTargetToSourceEdge_hash)

    def test___create_edges_negative_to_negative(self):

        class fakeNode:
            def __init__(self,
                        direction):
                self.direction = direction
            def get_direction(self):
                return self.direction

        # setup
        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        mock_sourceNode = fakeNode(-1)
        mock_targetNode = fakeNode(-1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        # execution
        actual_sourceToTargetEdge, actual_reverseTargetToSourceEdge = graph.create_edges(mock_sourceNode,
                                                                                        mock_targetNode,
                                                                                        mock_sourceDirection,
                                                                                        mock_targetDirection)
        actual_sourceToTargetEdge_sourceNode = actual_sourceToTargetEdge.get_sourceNode()
        actual_sourceToTargetEdge_targetNode = actual_sourceToTargetEdge.get_targetNode()
        actual_sourceToTargetEdge_sourceDirection = actual_sourceToTargetEdge.get_sourceNodeDirection()
        actual_sourceToTargetEdge_targetDirection = actual_sourceToTargetEdge.get_targetNodeDirection()
        actual_sourceToTargetEdge_coverage = actual_sourceToTargetEdge.get_edge_coverage()
        expected_sourceToTargetEdge_sourceNode = mock_sourceNode
        expected_sourceToTargetEdge_targetNode = mock_targetNode
        expected_sourceToTargetEdge_sourceDirection = mock_sourceDirection
        expected_sourceToTargetEdge_targetDirection = mock_targetDirection
        expected_sourceToTargetEdge_coverage = 0

        actual_reverseTargetToSourceEdge_sourceNode = actual_reverseTargetToSourceEdge.get_sourceNode()
        actual_reverseTargetToSourceEdge_targetNode = actual_reverseTargetToSourceEdge.get_targetNode()
        actual_reverseTargetToSourceEdge_sourceDirection = actual_reverseTargetToSourceEdge.get_sourceNodeDirection()
        actual_reverseTargetToSourceEdge_targetDirection = actual_reverseTargetToSourceEdge.get_targetNodeDirection()
        actual_reverseTargetToSourceEdge_coverage = actual_reverseTargetToSourceEdge.get_edge_coverage()
        expected_reverseTargetToSourceEdge_sourceNode = mock_targetNode
        expected_reverseTargetToSourceEdge_targetNode = mock_sourceNode
        expected_reverseTargetToSourceEdge_sourceDirection = mock_targetDirection * -1
        expected_reverseTargetToSourceEdge_targetDirection = mock_sourceDirection * -1
        expected_reverseTargetToSourceEdge_coverage = 0

        actual_sourceToTargetEdge_hash = actual_sourceToTargetEdge.__hash__()
        actual_reverseTargetToSourceEdge_hash = actual_reverseTargetToSourceEdge.__hash__()
        # assertion
        self.assertEqual(actual_sourceToTargetEdge_sourceNode, expected_sourceToTargetEdge_sourceNode)
        self.assertEqual(actual_sourceToTargetEdge_targetNode, expected_sourceToTargetEdge_targetNode)
        self.assertEqual(actual_sourceToTargetEdge_sourceDirection, expected_sourceToTargetEdge_sourceDirection)
        self.assertEqual(actual_sourceToTargetEdge_targetDirection, expected_sourceToTargetEdge_targetDirection)
        self.assertEqual(actual_sourceToTargetEdge_coverage, expected_sourceToTargetEdge_coverage)

        self.assertEqual(actual_reverseTargetToSourceEdge_sourceNode, expected_reverseTargetToSourceEdge_sourceNode)
        self.assertEqual(actual_reverseTargetToSourceEdge_targetNode, expected_reverseTargetToSourceEdge_targetNode)
        self.assertEqual(actual_reverseTargetToSourceEdge_sourceDirection, expected_reverseTargetToSourceEdge_sourceDirection)
        self.assertEqual(actual_reverseTargetToSourceEdge_targetDirection, expected_reverseTargetToSourceEdge_targetDirection)
        self.assertEqual(actual_reverseTargetToSourceEdge_coverage, expected_reverseTargetToSourceEdge_coverage)

        self.assertNotEqual(actual_sourceToTargetEdge_hash, actual_reverseTargetToSourceEdge_hash)

    def test___create_edges_positive_to_negative(self):

        class fakeNode:
            def __init__(self,
                        direction):
                self.direction = direction
            def get_direction(self):
                return self.direction

        # setup
        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        mock_sourceNode = fakeNode(1)
        mock_targetNode = fakeNode(-1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        # execution
        actual_sourceToTargetEdge, actual_reverseTargetToSourceEdge = graph.create_edges(mock_sourceNode,
                                                                                        mock_targetNode,
                                                                                        mock_sourceDirection,
                                                                                        mock_targetDirection)
        actual_sourceToTargetEdge_sourceNode = actual_sourceToTargetEdge.get_sourceNode()
        actual_sourceToTargetEdge_targetNode = actual_sourceToTargetEdge.get_targetNode()
        actual_sourceToTargetEdge_sourceDirection = actual_sourceToTargetEdge.get_sourceNodeDirection()
        actual_sourceToTargetEdge_targetDirection = actual_sourceToTargetEdge.get_targetNodeDirection()
        actual_sourceToTargetEdge_coverage = actual_sourceToTargetEdge.get_edge_coverage()
        expected_sourceToTargetEdge_sourceNode = mock_sourceNode
        expected_sourceToTargetEdge_targetNode = mock_targetNode
        expected_sourceToTargetEdge_sourceDirection = mock_sourceDirection
        expected_sourceToTargetEdge_targetDirection = mock_targetDirection
        expected_sourceToTargetEdge_coverage = 0

        actual_reverseTargetToSourceEdge_sourceNode = actual_reverseTargetToSourceEdge.get_sourceNode()
        actual_reverseTargetToSourceEdge_targetNode = actual_reverseTargetToSourceEdge.get_targetNode()
        actual_reverseTargetToSourceEdge_sourceDirection = actual_reverseTargetToSourceEdge.get_sourceNodeDirection()
        actual_reverseTargetToSourceEdge_targetDirection = actual_reverseTargetToSourceEdge.get_targetNodeDirection()
        actual_reverseTargetToSourceEdge_coverage = actual_reverseTargetToSourceEdge.get_edge_coverage()
        expected_reverseTargetToSourceEdge_sourceNode = mock_targetNode
        expected_reverseTargetToSourceEdge_targetNode = mock_sourceNode
        expected_reverseTargetToSourceEdge_sourceDirection = mock_targetDirection * -1
        expected_reverseTargetToSourceEdge_targetDirection = mock_sourceDirection * -1
        expected_reverseTargetToSourceEdge_coverage = 0

        actual_sourceToTargetEdge_hash = actual_sourceToTargetEdge.__hash__()
        actual_reverseTargetToSourceEdge_hash = actual_reverseTargetToSourceEdge.__hash__()
        # assertion
        self.assertEqual(actual_sourceToTargetEdge_sourceNode, expected_sourceToTargetEdge_sourceNode)
        self.assertEqual(actual_sourceToTargetEdge_targetNode, expected_sourceToTargetEdge_targetNode)
        self.assertEqual(actual_sourceToTargetEdge_sourceDirection, expected_sourceToTargetEdge_sourceDirection)
        self.assertEqual(actual_sourceToTargetEdge_targetDirection, expected_sourceToTargetEdge_targetDirection)
        self.assertEqual(actual_sourceToTargetEdge_coverage, expected_sourceToTargetEdge_coverage)

        self.assertEqual(actual_reverseTargetToSourceEdge_sourceNode, expected_reverseTargetToSourceEdge_sourceNode)
        self.assertEqual(actual_reverseTargetToSourceEdge_targetNode, expected_reverseTargetToSourceEdge_targetNode)
        self.assertEqual(actual_reverseTargetToSourceEdge_sourceDirection, expected_reverseTargetToSourceEdge_sourceDirection)
        self.assertEqual(actual_reverseTargetToSourceEdge_targetDirection, expected_reverseTargetToSourceEdge_targetDirection)
        self.assertEqual(actual_reverseTargetToSourceEdge_coverage, expected_reverseTargetToSourceEdge_coverage)

        self.assertNotEqual(actual_sourceToTargetEdge_hash, actual_reverseTargetToSourceEdge_hash)

    def test___create_edges_negative_to_positive(self):

        class fakeNode:
            def __init__(self,
                        direction):
                self.direction = direction
            def get_direction(self):
                return self.direction

        # setup
        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        mock_sourceNode = fakeNode(1)
        mock_targetNode = fakeNode(-1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        # execution
        actual_sourceToTargetEdge, actual_reverseTargetToSourceEdge = graph.create_edges(mock_sourceNode,
                                                                                        mock_targetNode,
                                                                                        mock_sourceDirection,
                                                                                        mock_targetDirection)
        actual_sourceToTargetEdge_sourceNode = actual_sourceToTargetEdge.get_sourceNode()
        actual_sourceToTargetEdge_targetNode = actual_sourceToTargetEdge.get_targetNode()
        actual_sourceToTargetEdge_sourceDirection = actual_sourceToTargetEdge.get_sourceNodeDirection()
        actual_sourceToTargetEdge_targetDirection = actual_sourceToTargetEdge.get_targetNodeDirection()
        actual_sourceToTargetEdge_coverage = actual_sourceToTargetEdge.get_edge_coverage()
        expected_sourceToTargetEdge_sourceNode = mock_sourceNode
        expected_sourceToTargetEdge_targetNode = mock_targetNode
        expected_sourceToTargetEdge_sourceDirection = mock_sourceDirection
        expected_sourceToTargetEdge_targetDirection = mock_targetDirection
        expected_sourceToTargetEdge_coverage = 0

        actual_reverseTargetToSourceEdge_sourceNode = actual_reverseTargetToSourceEdge.get_sourceNode()
        actual_reverseTargetToSourceEdge_targetNode = actual_reverseTargetToSourceEdge.get_targetNode()
        actual_reverseTargetToSourceEdge_sourceDirection = actual_reverseTargetToSourceEdge.get_sourceNodeDirection()
        actual_reverseTargetToSourceEdge_targetDirection = actual_reverseTargetToSourceEdge.get_targetNodeDirection()
        actual_reverseTargetToSourceEdge_coverage = actual_reverseTargetToSourceEdge.get_edge_coverage()
        expected_reverseTargetToSourceEdge_sourceNode = mock_targetNode
        expected_reverseTargetToSourceEdge_targetNode = mock_sourceNode
        expected_reverseTargetToSourceEdge_sourceDirection = mock_targetDirection * -1
        expected_reverseTargetToSourceEdge_targetDirection = mock_sourceDirection * -1
        expected_reverseTargetToSourceEdge_coverage = 0

        actual_sourceToTargetEdge_hash = actual_sourceToTargetEdge.__hash__()
        actual_reverseTargetToSourceEdge_hash = actual_reverseTargetToSourceEdge.__hash__()
        # assertion
        self.assertEqual(actual_sourceToTargetEdge_sourceNode, expected_sourceToTargetEdge_sourceNode)
        self.assertEqual(actual_sourceToTargetEdge_targetNode, expected_sourceToTargetEdge_targetNode)
        self.assertEqual(actual_sourceToTargetEdge_sourceDirection, expected_sourceToTargetEdge_sourceDirection)
        self.assertEqual(actual_sourceToTargetEdge_targetDirection, expected_sourceToTargetEdge_targetDirection)
        self.assertEqual(actual_sourceToTargetEdge_coverage, expected_sourceToTargetEdge_coverage)

        self.assertEqual(actual_reverseTargetToSourceEdge_sourceNode, expected_reverseTargetToSourceEdge_sourceNode)
        self.assertEqual(actual_reverseTargetToSourceEdge_targetNode, expected_reverseTargetToSourceEdge_targetNode)
        self.assertEqual(actual_reverseTargetToSourceEdge_sourceDirection, expected_reverseTargetToSourceEdge_sourceDirection)
        self.assertEqual(actual_reverseTargetToSourceEdge_targetDirection, expected_reverseTargetToSourceEdge_targetDirection)
        self.assertEqual(actual_reverseTargetToSourceEdge_coverage, expected_reverseTargetToSourceEdge_coverage)

        self.assertNotEqual(actual_sourceToTargetEdge_hash, actual_reverseTargetToSourceEdge_hash)

    def test___add_missing_edge_to_edges(self):
        # setup

        class fakeNode:
            def __init__(self,
                        direction):
                self.direction = direction
            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        mock_sourceNode = fakeNode(-1)
        mock_targetNode = fakeNode(1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(mock_sourceNode,
                                                                        mock_targetNode,
                                                                        mock_sourceDirection,
                                                                        mock_targetDirection)
        # execution
        actual_edge = graph.add_edge_to_edges(sourceToTargetEdge)
        actual_reverse_edge = graph.add_edge_to_edges(reverseTargetToSourceEdge)
        actual_graph_edges = graph.get_edges()
        actual_number_of_edges = len(actual_graph_edges)
        actual_edge_coverage = actual_edge.get_edge_coverage()
        actual_reverse_edge_coverage = actual_reverse_edge.get_edge_coverage()
        # assertion
        expected_number_of_edges = 2
        expected_edge = actual_graph_edges[sourceToTargetEdge.__hash__()]
        expected_reverse_edge = actual_graph_edges[reverseTargetToSourceEdge.__hash__()]
        expected_edge_coverage = 1
        self.assertEqual(actual_number_of_edges, expected_number_of_edges)
        self.assertEqual(actual_edge, expected_edge)
        self.assertEqual(actual_reverse_edge, expected_reverse_edge)
        self.assertEqual(actual_edge_coverage, expected_edge_coverage)
        self.assertEqual(actual_reverse_edge_coverage, expected_edge_coverage)

    def test___add_duplicate_edge_to_edges(self):
        # setup

        class fakeNode:
            def __init__(self,
                        direction):
                self.direction = direction
            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        mock_sourceNode = fakeNode(-1)
        mock_targetNode = fakeNode(1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(mock_sourceNode,
                                                                        mock_targetNode,
                                                                        mock_sourceDirection,
                                                                        mock_targetDirection)
        graph.add_edge_to_edges(sourceToTargetEdge)
        graph.add_edge_to_edges(reverseTargetToSourceEdge)
        # execution
        actual_edge = graph.add_edge_to_edges(sourceToTargetEdge)
        actual_reverse_edge = graph.add_edge_to_edges(reverseTargetToSourceEdge)
        actual_graph_edges = graph.get_edges()
        actual_number_of_edges = len(actual_graph_edges)
        actual_edge_coverage = actual_edge.get_edge_coverage()
        actual_reverse_edge_coverage = actual_reverse_edge.get_edge_coverage()
        # assertion
        expected_number_of_edges = 2
        expected_edge = actual_graph_edges[sourceToTargetEdge.__hash__()]
        expected_reverse_edge = actual_graph_edges[reverseTargetToSourceEdge.__hash__()]
        expected_edge_coverage = 2
        self.assertEqual(actual_number_of_edges, expected_number_of_edges)
        self.assertEqual(actual_edge, expected_edge)
        self.assertEqual(actual_reverse_edge, expected_reverse_edge)
        self.assertEqual(actual_edge_coverage, expected_edge_coverage)
        self.assertEqual(actual_reverse_edge_coverage, expected_edge_coverage)

    def test___add_third_edge_to_edges(self):
        # setup

        class fakeNode:
            def __init__(self,
                        direction):
                self.direction = direction
            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        mock_sourceNode = fakeNode(-1)
        mock_targetNode = fakeNode(1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(mock_sourceNode,
                                                                        mock_targetNode,
                                                                        mock_sourceDirection,
                                                                        mock_targetDirection)
        graph.add_edge_to_edges(sourceToTargetEdge)
        graph.add_edge_to_edges(reverseTargetToSourceEdge)
        graph.add_edge_to_edges(sourceToTargetEdge)
        graph.add_edge_to_edges(reverseTargetToSourceEdge)
        # execution
        actual_edge = graph.add_edge_to_edges(sourceToTargetEdge)
        actual_reverse_edge = graph.add_edge_to_edges(reverseTargetToSourceEdge)
        actual_graph_edges = graph.get_edges()
        actual_number_of_edges = len(actual_graph_edges)
        actual_edge_coverage = actual_edge.get_edge_coverage()
        actual_reverse_edge_coverage = actual_reverse_edge.get_edge_coverage()
        # assertion
        expected_number_of_edges = 2
        expected_edge = actual_graph_edges[sourceToTargetEdge.__hash__()]
        expected_reverse_edge = actual_graph_edges[reverseTargetToSourceEdge.__hash__()]
        expected_edge_coverage = 3
        self.assertEqual(actual_number_of_edges, expected_number_of_edges)
        self.assertEqual(actual_edge, expected_edge)
        self.assertEqual(actual_reverse_edge, expected_reverse_edge)
        self.assertEqual(actual_edge_coverage, expected_edge_coverage)
        self.assertEqual(actual_reverse_edge_coverage, expected_edge_coverage)

    def test___add_new_edge_to_graph(self):
        # setup

        class fakeNode:
            def __init__(self,
                        direction):
                self.direction = direction
            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        mock_sourceNode = fakeNode(-1)
        mock_targetNode = fakeNode(1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(mock_sourceNode,
                                                                        mock_targetNode,
                                                                        mock_sourceDirection,
                                                                        mock_targetDirection)
        # execution
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(sourceToTargetEdge,
                                                                                reverseTargetToSourceEdge)
        actual_graph_edges = graph.get_edges()
        actual_number_of_edges = len(actual_graph_edges)
        # assertion
        expected_number_of_edges = 2
        expected_edge_hashes = [sourceToTargetEdge.__hash__(), reverseTargetToSourceEdge.__hash__()]
        expected_edge_coverage = 1
        self.assertEqual(actual_number_of_edges, expected_number_of_edges)
        self.assertTrue(all(h in list(actual_graph_edges.keys()) for h in expected_edge_hashes))
        self.assertTrue(all(actual_graph_edges[h].get_edge_coverage() == expected_edge_coverage for h in expected_edge_hashes))

    def test___add_duplicate_edges_to_graph(self):
        # setup

        class fakeNode:
            def __init__(self,
                        direction):
                self.direction = direction
            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        mock_sourceNode = fakeNode(-1)
        mock_targetNode = fakeNode(1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(mock_sourceNode,
                                                                        mock_targetNode,
                                                                        mock_sourceDirection,
                                                                        mock_targetDirection)
        graph.add_edge_to_edges(sourceToTargetEdge)
        graph.add_edge_to_edges(reverseTargetToSourceEdge)
        # execution
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(sourceToTargetEdge,
                                                                                reverseTargetToSourceEdge)
        actual_graph_edges = graph.get_edges()
        actual_number_of_edges = len(actual_graph_edges)
        # assertion
        expected_number_of_edges = 2
        expected_edge_hashes = [sourceToTargetEdge.__hash__(), reverseTargetToSourceEdge.__hash__()]
        expected_edge_coverage = 2
        self.assertEqual(actual_number_of_edges, expected_number_of_edges)
        self.assertTrue(all(h in list(actual_graph_edges.keys()) for h in expected_edge_hashes))
        self.assertTrue(all(actual_graph_edges[h].get_edge_coverage() == expected_edge_coverage for h in expected_edge_hashes))

    def test___add_one_new_one_duplicate_edges_to_graph(self):
        # setup

        class fakeNode:
            def __init__(self,
                        id,
                        direction):
                self.id = id
                self.direction = direction
            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        mock_sourceNode = fakeNode(0, -1)
        mock_targetNode = fakeNode(1, 1)
        mock_thirdNode = fakeNode(1, -1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        mock_thirdDirection = mock_thirdNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(mock_sourceNode,
                                                                        mock_targetNode,
                                                                        mock_sourceDirection,
                                                                        mock_targetDirection)
        targetToThirdEdge, reverseThirdToTargetEdge = graph.create_edges(mock_targetNode,
                                                                        mock_thirdNode,
                                                                        mock_targetDirection,
                                                                        mock_thirdDirection)
        graph.add_edge_to_edges(sourceToTargetEdge)
        graph.add_edge_to_edges(reverseTargetToSourceEdge)
        # execution
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(sourceToTargetEdge,
                                                                                reverseTargetToSourceEdge)
        targetToThirdEdge, reverseThirdToTargetEdge = graph.add_edges_to_graph(targetToThirdEdge,
                                                                                reverseThirdToTargetEdge)
        actual_graph_edges = graph.get_edges()
        actual_number_of_edges = len(actual_graph_edges)
        # assertion
        expected_number_of_edges = 4
        expected_edge_hashes = [sourceToTargetEdge.__hash__(), reverseTargetToSourceEdge.__hash__(), targetToThirdEdge.__hash__(), reverseThirdToTargetEdge.__hash__()]
        expected_source_target_edge_coverage = 2
        expected_target_third_edge_coverage = 1
        self.assertEqual(actual_number_of_edges, expected_number_of_edges)
        self.assertTrue(all(h in list(actual_graph_edges.keys()) for h in expected_edge_hashes))
        self.assertEqual(sourceToTargetEdge.get_edge_coverage(), expected_source_target_edge_coverage)
        self.assertEqual(reverseTargetToSourceEdge.get_edge_coverage(), expected_source_target_edge_coverage)
        self.assertEqual(targetToThirdEdge.get_edge_coverage(), expected_target_third_edge_coverage)
        self.assertEqual(reverseThirdToTargetEdge.get_edge_coverage(), expected_target_third_edge_coverage)

    def test___add_two_duplicate_reverse_edges_to_graph_all_positive(self):
        # setup

        class fakeNode:
            def __init__(self,
                        id,
                        direction):
                self.id = id
                self.direction = direction
            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        mock_sourceNode = fakeNode(0, 1)
        mock_targetNode = fakeNode(1, 1)
        mock_thirdNode = fakeNode(1, 1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        mock_thirdDirection = mock_thirdNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(mock_sourceNode,
                                                                        mock_targetNode,
                                                                        mock_sourceDirection,
                                                                        mock_targetDirection)
        targetToThirdEdge, reverseThirdToTargetEdge = graph.create_edges(mock_targetNode,
                                                                        mock_thirdNode,
                                                                        mock_targetDirection,
                                                                        mock_thirdDirection)
        graph.add_edge_to_edges(sourceToTargetEdge)
        graph.add_edge_to_edges(reverseTargetToSourceEdge)
        graph.add_edge_to_edges(targetToThirdEdge)
        graph.add_edge_to_edges(reverseThirdToTargetEdge)
        # execution
        rc_sourceToTargetEdge = sourceToTargetEdge
        rc_sourceToTargetEdge.set_sourceNodeDirection(sourceToTargetEdge.get_sourceNodeDirection() * -1)
        rc_sourceToTargetEdge.set_targetNodeDirection(sourceToTargetEdge.get_targetNodeDirection() * -1)
        rc_reverseTargetToSourceEdge = reverseTargetToSourceEdge
        rc_reverseTargetToSourceEdge.set_sourceNodeDirection(reverseTargetToSourceEdge.get_sourceNodeDirection() * -1)
        rc_reverseTargetToSourceEdge.set_targetNodeDirection(reverseTargetToSourceEdge.get_targetNodeDirection() * -1)
        rc_targetToThirdEdge = targetToThirdEdge
        rc_targetToThirdEdge.set_sourceNodeDirection(targetToThirdEdge.get_sourceNodeDirection() * -1)
        rc_targetToThirdEdge.set_targetNodeDirection(targetToThirdEdge.get_targetNodeDirection() * -1)
        rc_reverseThirdToTargetEdge = reverseThirdToTargetEdge
        rc_reverseThirdToTargetEdge.set_sourceNodeDirection(reverseThirdToTargetEdge.get_sourceNodeDirection() * -1)
        reverseThirdToTargetEdge.set_targetNodeDirection(reverseThirdToTargetEdge.get_targetNodeDirection() * -1)
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(rc_sourceToTargetEdge,
                                                                                rc_reverseTargetToSourceEdge)
        targetToThirdEdge, reverseThirdToTargetEdge = graph.add_edges_to_graph(rc_targetToThirdEdge,
                                                                            rc_reverseThirdToTargetEdge)
        actual_graph_edges = graph.get_edges()
        actual_number_of_edges = len(actual_graph_edges)
        # assertion
        expected_number_of_edges = 4
        expected_edge_hashes = [sourceToTargetEdge.__hash__(), reverseTargetToSourceEdge.__hash__(), targetToThirdEdge.__hash__(), reverseThirdToTargetEdge.__hash__()]
        expected_source_target_edge_coverage = 2
        expected_target_third_edge_coverage = 2
        self.assertEqual(actual_number_of_edges, expected_number_of_edges)
        self.assertTrue(all(h in list(actual_graph_edges.keys()) for h in expected_edge_hashes))
        self.assertEqual(sourceToTargetEdge.get_edge_coverage(), expected_source_target_edge_coverage)
        self.assertEqual(reverseTargetToSourceEdge.get_edge_coverage(), expected_source_target_edge_coverage)
        self.assertEqual(targetToThirdEdge.get_edge_coverage(), expected_target_third_edge_coverage)
        self.assertEqual(reverseThirdToTargetEdge.get_edge_coverage(), expected_target_third_edge_coverage)

    def test___add_two_duplicate_reverse_edges_to_graph_all_negative(self):
        # setup

        class fakeNode:
            def __init__(self,
                        id,
                        direction):
                self.id = id
                self.direction = direction
            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        mock_sourceNode = fakeNode(0,-1)
        mock_targetNode = fakeNode(1, -1)
        mock_thirdNode = fakeNode(1, -1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        mock_thirdDirection = mock_thirdNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(mock_sourceNode,
                                                                        mock_targetNode,
                                                                        mock_sourceDirection,
                                                                        mock_targetDirection)
        targetToThirdEdge, reverseThirdToTargetEdge = graph.create_edges(mock_targetNode,
                                                                        mock_thirdNode,
                                                                        mock_targetDirection,
                                                                        mock_thirdDirection)
        graph.add_edge_to_edges(sourceToTargetEdge)
        graph.add_edge_to_edges(reverseTargetToSourceEdge)
        graph.add_edge_to_edges(targetToThirdEdge)
        graph.add_edge_to_edges(reverseThirdToTargetEdge)
        # execution
        rc_sourceToTargetEdge = sourceToTargetEdge
        rc_sourceToTargetEdge.set_sourceNodeDirection(sourceToTargetEdge.get_sourceNodeDirection() * -1)
        rc_sourceToTargetEdge.set_targetNodeDirection(sourceToTargetEdge.get_targetNodeDirection() * -1)
        rc_reverseTargetToSourceEdge = reverseTargetToSourceEdge
        rc_reverseTargetToSourceEdge.set_sourceNodeDirection(reverseTargetToSourceEdge.get_sourceNodeDirection() * -1)
        rc_reverseTargetToSourceEdge.set_targetNodeDirection(reverseTargetToSourceEdge.get_targetNodeDirection() * -1)
        rc_targetToThirdEdge = targetToThirdEdge
        rc_targetToThirdEdge.set_sourceNodeDirection(targetToThirdEdge.get_sourceNodeDirection() * -1)
        rc_targetToThirdEdge.set_targetNodeDirection(targetToThirdEdge.get_targetNodeDirection() * -1)
        rc_reverseThirdToTargetEdge = reverseThirdToTargetEdge
        rc_reverseThirdToTargetEdge.set_sourceNodeDirection(reverseThirdToTargetEdge.get_sourceNodeDirection() * -1)
        reverseThirdToTargetEdge.set_targetNodeDirection(reverseThirdToTargetEdge.get_targetNodeDirection() * -1)
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(rc_sourceToTargetEdge,
                                                                                rc_reverseTargetToSourceEdge)
        targetToThirdEdge, reverseThirdToTargetEdge = graph.add_edges_to_graph(rc_targetToThirdEdge,
                                                                            rc_reverseThirdToTargetEdge)
        actual_graph_edges = graph.get_edges()
        actual_number_of_edges = len(actual_graph_edges)
        # assertion
        expected_number_of_edges = 4
        expected_edge_hashes = [sourceToTargetEdge.__hash__(), reverseTargetToSourceEdge.__hash__(), targetToThirdEdge.__hash__(), reverseThirdToTargetEdge.__hash__()]
        expected_source_target_edge_coverage = 2
        expected_target_third_edge_coverage = 2
        self.assertEqual(actual_number_of_edges, expected_number_of_edges)
        self.assertTrue(all(h in list(actual_graph_edges.keys()) for h in expected_edge_hashes))
        self.assertEqual(sourceToTargetEdge.get_edge_coverage(), expected_source_target_edge_coverage)
        self.assertEqual(reverseTargetToSourceEdge.get_edge_coverage(), expected_source_target_edge_coverage)
        self.assertEqual(targetToThirdEdge.get_edge_coverage(), expected_target_third_edge_coverage)
        self.assertEqual(reverseThirdToTargetEdge.get_edge_coverage(), expected_target_third_edge_coverage)

    def test___add_two_duplicate_reverse_edges_to_graph_one_positive_one_negative(self):
        # setup

        class fakeNode:
            def __init__(self,
                        id,
                        direction):
                self.id = id
                self.direction = direction
            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        mock_sourceNode = fakeNode(0,1)
        mock_targetNode = fakeNode(1, 1)
        mock_thirdNode = fakeNode(1, -1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        mock_thirdDirection = mock_thirdNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(mock_sourceNode,
                                                                        mock_targetNode,
                                                                        mock_sourceDirection,
                                                                        mock_targetDirection)
        targetToThirdEdge, reverseThirdToTargetEdge = graph.create_edges(mock_targetNode,
                                                                        mock_thirdNode,
                                                                        mock_targetDirection,
                                                                        mock_thirdDirection)
        graph.add_edge_to_edges(sourceToTargetEdge)
        graph.add_edge_to_edges(reverseTargetToSourceEdge)
        graph.add_edge_to_edges(targetToThirdEdge)
        graph.add_edge_to_edges(reverseThirdToTargetEdge)
        # execution
        rc_targetToThirdEdge = targetToThirdEdge
        rc_targetToThirdEdge.set_sourceNodeDirection(targetToThirdEdge.get_sourceNodeDirection() * -1)
        rc_targetToThirdEdge.set_targetNodeDirection(targetToThirdEdge.get_targetNodeDirection() * -1)
        rc_reverseThirdToTargetEdge = reverseThirdToTargetEdge
        rc_reverseThirdToTargetEdge.set_sourceNodeDirection(reverseThirdToTargetEdge.get_sourceNodeDirection() * -1)
        reverseThirdToTargetEdge.set_targetNodeDirection(reverseThirdToTargetEdge.get_targetNodeDirection() * -1)
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(sourceToTargetEdge,
                                                                                reverseTargetToSourceEdge)
        targetToThirdEdge, reverseThirdToTargetEdge = graph.add_edges_to_graph(rc_targetToThirdEdge,
                                                                            rc_reverseThirdToTargetEdge)
        actual_graph_edges = graph.get_edges()
        actual_number_of_edges = len(actual_graph_edges)
        # assertion
        expected_number_of_edges = 4
        expected_edge_hashes = [sourceToTargetEdge.__hash__(), reverseTargetToSourceEdge.__hash__(), targetToThirdEdge.__hash__(), reverseThirdToTargetEdge.__hash__()]
        expected_source_target_edge_coverage = 2
        expected_target_third_edge_coverage = 2
        self.assertEqual(actual_number_of_edges, expected_number_of_edges)
        self.assertTrue(all(h in list(actual_graph_edges.keys()) for h in expected_edge_hashes))
        self.assertEqual(sourceToTargetEdge.get_edge_coverage(), expected_source_target_edge_coverage)
        self.assertEqual(reverseTargetToSourceEdge.get_edge_coverage(), expected_source_target_edge_coverage)
        self.assertEqual(targetToThirdEdge.get_edge_coverage(), expected_target_third_edge_coverage)
        self.assertEqual(reverseThirdToTargetEdge.get_edge_coverage(), expected_target_third_edge_coverage)

    def test___add_two_duplicate_reverse_edges_to_graph_one_negative_one_positive(self):
        # setup

        class fakeNode:
            def __init__(self,
                        id,
                        direction):
                self.id = id
                self.direction = direction
            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({},
                            3,
                            1,
                            1)
        mock_sourceNode = fakeNode(0,-1)
        mock_targetNode = fakeNode(1, 1)
        mock_thirdNode = fakeNode(1, 1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        mock_thirdDirection = mock_thirdNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(mock_sourceNode,
                                                                        mock_targetNode,
                                                                        mock_sourceDirection,
                                                                        mock_targetDirection)
        targetToThirdEdge, reverseThirdToTargetEdge = graph.create_edges(mock_targetNode,
                                                                        mock_thirdNode,
                                                                        mock_targetDirection,
                                                                        mock_thirdDirection)
        graph.add_edge_to_edges(sourceToTargetEdge)
        graph.add_edge_to_edges(reverseTargetToSourceEdge)
        graph.add_edge_to_edges(targetToThirdEdge)
        graph.add_edge_to_edges(reverseThirdToTargetEdge)
        # execution
        rc_sourceToTargetEdge = sourceToTargetEdge
        rc_sourceToTargetEdge.set_sourceNodeDirection(sourceToTargetEdge.get_sourceNodeDirection() * -1)
        rc_sourceToTargetEdge.set_targetNodeDirection(sourceToTargetEdge.get_targetNodeDirection() * -1)
        rc_reverseTargetToSourceEdge = reverseTargetToSourceEdge
        rc_reverseTargetToSourceEdge.set_sourceNodeDirection(reverseTargetToSourceEdge.get_sourceNodeDirection() * -1)
        rc_reverseTargetToSourceEdge.set_targetNodeDirection(reverseTargetToSourceEdge.get_targetNodeDirection() * -1)
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(rc_sourceToTargetEdge,
                                                                                rc_reverseTargetToSourceEdge)
        targetToThirdEdge, reverseThirdToTargetEdge = graph.add_edges_to_graph(targetToThirdEdge,
                                                                            reverseThirdToTargetEdge)
        actual_graph_edges = graph.get_edges()
        actual_number_of_edges = len(actual_graph_edges)
        # assertion
        expected_number_of_edges = 4
        expected_edge_hashes = [sourceToTargetEdge.__hash__(), reverseTargetToSourceEdge.__hash__(), targetToThirdEdge.__hash__(), reverseThirdToTargetEdge.__hash__()]
        expected_source_target_edge_coverage = 2
        expected_target_third_edge_coverage = 2
        self.assertEqual(actual_number_of_edges, expected_number_of_edges)
        self.assertTrue(all(h in list(actual_graph_edges.keys()) for h in expected_edge_hashes))
        self.assertEqual(sourceToTargetEdge.get_edge_coverage(), expected_source_target_edge_coverage)
        self.assertEqual(reverseTargetToSourceEdge.get_edge_coverage(), expected_source_target_edge_coverage)
        self.assertEqual(targetToThirdEdge.get_edge_coverage(), expected_target_third_edge_coverage)
        self.assertEqual(reverseThirdToTargetEdge.get_edge_coverage(), expected_target_third_edge_coverage)