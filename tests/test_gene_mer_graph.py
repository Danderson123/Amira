import unittest
import sys
sys.path.insert(0, "amira_prototype")

from construct_graph import GeneMerGraph
from construct_read import Read
from construct_node import Node
from construct_edge import Edge

class TestGeneMerGraphConstructor(unittest.TestCase):

    def test___init_empty_GeneMerGraph(self):
        # setup
        graph = GeneMerGraph({},
                            0)
        # execution
        actual_reads = graph.get_reads()
        actual_kmerSize = graph.get_kmerSize()
        actual_minGeneMerCoverage = graph.get_minNodeCoverage()
        actual_minEdgeCoverage = graph.get_minEdgeCoverage()
        actual_nodes = graph.get_nodes()
        actual_edges = graph.get_edges()
        # assertion
        expected_reads = {}
        expected_kmerSize = 0
        expected_minGeneMerCoverage = 1
        expected_minEdgeCoverage = 1
        expected_nodes = {}
        expected_edges = {}
        self.assertEqual(actual_reads, expected_reads)
        self.assertEqual(actual_kmerSize, expected_kmerSize)
        self.assertEqual(actual_minGeneMerCoverage, expected_minGeneMerCoverage)
        self.assertEqual(actual_minEdgeCoverage, expected_minEdgeCoverage)
        self.assertEqual(actual_nodes, expected_nodes)
        self.assertEqual(actual_edges, expected_edges)

    def test___init_non_empty_GeneMerGraph(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene6"]
        graph = GeneMerGraph({"read1": genes,
                            "read2": genes2},
                            3)
        # execution
        actual_reads = graph.get_reads()
        actual_kmerSize = graph.get_kmerSize()
        actual_minGeneMerCoverage = graph.get_minNodeCoverage()
        actual_minEdgeCoverage = graph.get_minEdgeCoverage()
        actual_number_nodes = graph.get_total_number_of_nodes()
        actual_number_edges = graph.get_total_number_of_edges()
        actual_node_1_coverage = graph.get_nodes()[list(graph.get_nodes().keys())[0]].get_node_coverage()
        actual_node_2_coverage = graph.get_nodes()[list(graph.get_nodes().keys())[1]].get_node_coverage()
        actual_node_3_coverage = graph.get_nodes()[list(graph.get_nodes().keys())[2]].get_node_coverage()
        # assertion
        expected_reads = {"read1": genes,
                        "read2": genes2}
        expected_kmerSize = 3
        expected_minGeneMerCoverage = 1
        expected_minEdgeCoverage = 1
        expected_number_nodes = 3
        expected_number_edges = 4
        expected_node_1_coverage = 2
        expected_node_2_coverage = 1
        expected_node_3_coverage = 1
        self.assertEqual(actual_reads, expected_reads)
        self.assertEqual(actual_kmerSize, expected_kmerSize)
        self.assertEqual(actual_minGeneMerCoverage, expected_minGeneMerCoverage)
        self.assertEqual(actual_minEdgeCoverage, expected_minEdgeCoverage)
        self.assertEqual(actual_number_edges, expected_number_edges)
        self.assertEqual(actual_number_nodes, expected_number_nodes)
        self.assertEqual(actual_node_1_coverage, expected_node_1_coverage)
        self.assertEqual(actual_node_2_coverage, expected_node_2_coverage)
        self.assertEqual(actual_node_3_coverage, expected_node_3_coverage)

    def test___init_non_empty_duplicate_nodes_GeneMerGraph(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene1", "-gene2", "+gene3", "+gene8"] # 5 nodes, 10 edges
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene6", "+gene1", "-gene2", "+gene3"] # 3 nodes, 8 edges
        graph = GeneMerGraph({"read1": genes,
                            "read2": genes2},
                            3)
        # execution
        actual_reads = graph.get_reads()
        actual_kmerSize = graph.get_kmerSize()
        actual_minGeneMerCoverage = graph.get_minNodeCoverage()
        actual_minEdgeCoverage = graph.get_minEdgeCoverage()
        actual_number_nodes = graph.get_total_number_of_nodes()
        actual_number_edges = graph.get_total_number_of_edges()
        actual_node_1_coverage = graph.get_nodes()[list(graph.get_nodes().keys())[0]].get_node_coverage()
        # assertion
        expected_reads = {"read1": genes,
                        "read2": genes2}
        expected_kmerSize = 3
        expected_minGeneMerCoverage = 1
        expected_minEdgeCoverage = 1
        expected_number_nodes = 8
        expected_number_edges = 18
        expected_node_1_coverage = 4
        expected_edge_coverage = 1
        self.assertEqual(actual_reads, expected_reads)
        self.assertEqual(actual_kmerSize, expected_kmerSize)
        self.assertEqual(actual_minGeneMerCoverage, expected_minGeneMerCoverage)
        self.assertEqual(actual_minEdgeCoverage, expected_minEdgeCoverage)
        self.assertEqual(actual_number_edges, expected_number_edges)
        self.assertEqual(actual_number_nodes, expected_number_nodes)
        self.assertEqual(actual_node_1_coverage, expected_node_1_coverage)
        self.assertTrue(all(node.get_node_coverage() == 1 for node in list(graph.get_nodes().values())[1:]))
        [print(edge.get_edge_coverage()) for edge in list(graph.get_edges().values())]
        self.assertTrue(all(edge.get_edge_coverage() == expected_edge_coverage for edge in list(graph.get_edges().values())))

    def test__init_two_genemers_one_read(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"] # 2 nodes, 1 edge
        graph = GeneMerGraph({"read1": genes},
                            3)
        # execution
        actual_reads = graph.get_reads()
        actual_kmerSize = graph.get_kmerSize()
        actual_minGeneMerCoverage = graph.get_minNodeCoverage()
        actual_minEdgeCoverage = graph.get_minEdgeCoverage()
        actual_number_nodes = graph.get_total_number_of_nodes()
        actual_number_edges = graph.get_total_number_of_edges()
        actual_node_1_coverage = graph.get_nodes()[list(graph.get_nodes().keys())[0]].get_node_coverage()
        # assertion
        expected_reads = {"read1": genes}
        expected_kmerSize = 3
        expected_minGeneMerCoverage = 1
        expected_minEdgeCoverage = 1
        expected_number_nodes = 2
        expected_number_edges = 2
        expected_node_1_coverage = 1
        expected_edge_coverage = 1
        self.assertEqual(actual_reads, expected_reads)
        self.assertEqual(actual_kmerSize, expected_kmerSize)
        self.assertEqual(actual_minGeneMerCoverage, expected_minGeneMerCoverage)
        self.assertEqual(actual_minEdgeCoverage, expected_minEdgeCoverage)
        self.assertEqual(actual_number_edges, expected_number_edges)
        self.assertEqual(actual_number_nodes, expected_number_nodes)
        self.assertEqual(actual_node_1_coverage, expected_node_1_coverage)
        self.assertTrue(all(node.get_node_coverage() == 1 for node in list(graph.get_nodes().values())[1:]))
        self.assertTrue(all(edge.get_edge_coverage() == expected_edge_coverage for edge in list(graph.get_edges().values())))

    def test___all_nodes(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "-gene3", "+gene2", "-gene1"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                            3)
        for g in geneMers:
            graph.add_node(g,
                        [read1])
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
                        3)
        # execution
        actual_returned_node = graph.add_node(geneMer,
                                            [read1])
        actual_returned_node.increment_node_coverage()
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
                        3)
        node1 = graph.add_node(geneMer1,
                            [read1])
        node1.increment_node_coverage()
        genes2 = ["+gene4", "-gene3", "+gene2"]
        read2 = Read("read2",
                    genes2)
        geneMer2 = [x for x in read2.get_geneMers(3)][0]
        # execution
        actual_returned_node2 = graph.add_node(geneMer2,
                                            [read2])
        actual_returned_node2.increment_node_coverage()
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
                        3)
        graph.add_node(geneMer,
                    [read1])
        # execution
        actual_returned_node = graph.add_node(geneMer,
                                        [read1])
        actual_returned_node.increment_node_coverage()
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
                            3)
        node1 = graph.add_node(geneMer1,
                    [read1])
        node1.increment_node_coverage()
        genes2 = ["+gene4", "-gene3", "+gene2"]
        read2 = Read("read2",
                    genes2)
        geneMer2 = [x for x in read2.get_geneMers(3)][0]
        node2 = graph.add_node(geneMer2,
                    [read2])
        node2.increment_node_coverage()
        # execution
        actual_node1 = graph.get_node(geneMer1)
        actual_node1 = graph.get_node(geneMer1)
        actual_node1_read_list = [x for x in actual_node1.get_reads()]
        actual_node1_hash = actual_node1.__hash__()
        actual_node1_coverage = actual_node1.get_node_coverage()
        actual_node2 = graph.get_node(geneMer2)
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
                        3)
        graph.add_node(geneMer1,
                    [read1])
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
                            3)
        for g in geneMers:
            graph.add_node(g,
                    [read1])
        selectedGenes = ["gene2", "gene6"]
        for g in range(len(selectedGenes)):
            # execution
            actual_selectedNodes = [n for n in graph.get_nodes_containing(selectedGenes[g])]
            actual_selectedGenes = [[g.get_name() for g in n.get_canonical_geneMer()] for n in actual_selectedNodes]
            actual_geneMerCount = len(actual_selectedGenes)
            # assertion
            expected_geneMerCount = 3
            self.assertTrue(all(selectedGenes[g] in ng for ng in actual_selectedGenes))
            self.assertEqual(actual_geneMerCount, expected_geneMerCount)

    def test___get_nodes_containing_all(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "-gene3", "+gene2", "-gene1"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                            3)
        for g in geneMers:
            graph.add_node(g,
                    [read1])
        selectedGenes = ["gene1", "gene2", "gene3", "gene4", "gene5", "gene6", "gene3", "gene2", "gene1"]
        expected_counts = [1, 3, 5, 3, 3, 3, 5, 3, 1]
        for g in range(len(selectedGenes)):
            # execution
            actual_selectedNodes = [n for n in graph.get_nodes_containing(selectedGenes[g])]
            actual_selectedGenes = [[g.get_name() for g in n.get_canonical_geneMer()] for n in actual_selectedNodes]
            actual_geneMerCount = len(actual_selectedGenes)
            # assertion
            expected_geneMerCount = expected_counts[g]
            self.assertTrue(all(selectedGenes[g] in ng for ng in actual_selectedGenes))
            self.assertEqual(actual_geneMerCount, expected_geneMerCount)

    def test___get_nodes_containing_gene_not_in_graph(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "-gene3", "+gene2", "-gene1"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                            3)
        for g in geneMers:
            graph.add_node(g,
                    [read1])
        selectedGenes = "gene10"
        # execution
        actual_selectedNodes = [n for n in graph.get_nodes_containing(selectedGenes)]
        actual_selectedGenes = [[g.get_name() for g in n.get_canonical_geneMer()] for n in actual_selectedNodes]
        # assertion
        expected_selectedGenes = []
        self.assertEqual(actual_selectedGenes, expected_selectedGenes)

    def test___get_nodes_containing_strand(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "-gene3", "+gene2", "-gene1"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                            3)
        for g in geneMers:
            graph.add_node(g,
                    [read1])
        selectedGenes = ["gene2", "+gene6"]
        # assertion
        self.assertRaises(AssertionError, graph.get_nodes_containing, selectedGenes)

    def test___get_nodes_containing_strand(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "-gene3", "+gene2", "-gene1"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                            3)
        for g in geneMers:
            graph.add_node(g,
                    [read1])
        selectedGenes = ["+gene6"]
        for g in range(len(selectedGenes)):
            # assertion
            self.assertRaises(AssertionError, graph.get_nodes_containing, selectedGenes[g])

    def test___create_edges_positive_to_positive(self):

        class fakeNode:
            def __init__(self,
                        direction):
                self.direction = direction
            def get_direction(self):
                return self.direction

        # setup
        graph = GeneMerGraph({},
                            3)
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
                            3)
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
                            3)
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
                            3)
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
                            3)
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
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
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
                            3)
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
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        # execution
        actual_edge = graph.add_edge_to_edges(sourceToTargetEdge)
        actual_reverse_edge = graph.add_edge_to_edges(reverseTargetToSourceEdge)
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
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
                            3)
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
        sourceToTargetEdge.increment_edge_coverage()
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        # execution
        actual_edge = graph.add_edge_to_edges(sourceToTargetEdge)
        actual_reverse_edge = graph.add_edge_to_edges(reverseTargetToSourceEdge)
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
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
                            3)
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
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
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
                            3)
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
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        # execution
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(sourceToTargetEdge,
                                                                                reverseTargetToSourceEdge)
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
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
                            3)
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
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        # execution
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(sourceToTargetEdge,
                                                                                reverseTargetToSourceEdge)
        targetToThirdEdge, reverseThirdToTargetEdge = graph.add_edges_to_graph(targetToThirdEdge,
                                                                                reverseThirdToTargetEdge)
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        targetToThirdEdge.increment_edge_coverage()
        reverseThirdToTargetEdge.increment_edge_coverage()
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
                            3)
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
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        graph.add_edge_to_edges(targetToThirdEdge)
        graph.add_edge_to_edges(reverseThirdToTargetEdge)
        targetToThirdEdge.increment_edge_coverage()
        reverseThirdToTargetEdge.increment_edge_coverage()
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
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        targetToThirdEdge, reverseThirdToTargetEdge = graph.add_edges_to_graph(rc_targetToThirdEdge,
                                                                            rc_reverseThirdToTargetEdge)
        targetToThirdEdge.increment_edge_coverage()
        reverseThirdToTargetEdge.increment_edge_coverage()
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
                            3)
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
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        graph.add_edge_to_edges(targetToThirdEdge)
        graph.add_edge_to_edges(reverseThirdToTargetEdge)
        targetToThirdEdge.increment_edge_coverage()
        reverseThirdToTargetEdge.increment_edge_coverage()
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
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        targetToThirdEdge, reverseThirdToTargetEdge = graph.add_edges_to_graph(rc_targetToThirdEdge,
                                                                            rc_reverseThirdToTargetEdge)
        targetToThirdEdge.increment_edge_coverage()
        reverseThirdToTargetEdge.increment_edge_coverage()
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
        #### UNIT TEST NEEDED ####

        class fakeNode:
            def __init__(self,
                        id,
                        direction):
                self.id = id
                self.direction = direction
            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({},
                            3)
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
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        graph.add_edge_to_edges(targetToThirdEdge)
        graph.add_edge_to_edges(reverseThirdToTargetEdge)
        targetToThirdEdge.increment_edge_coverage()
        reverseThirdToTargetEdge.increment_edge_coverage()
        # execution
        rc_targetToThirdEdge = targetToThirdEdge
        rc_targetToThirdEdge.set_sourceNodeDirection(targetToThirdEdge.get_sourceNodeDirection() * -1)
        rc_targetToThirdEdge.set_targetNodeDirection(targetToThirdEdge.get_targetNodeDirection() * -1)
        rc_reverseThirdToTargetEdge = reverseThirdToTargetEdge
        rc_reverseThirdToTargetEdge.set_sourceNodeDirection(reverseThirdToTargetEdge.get_sourceNodeDirection() * -1)
        reverseThirdToTargetEdge.set_targetNodeDirection(reverseThirdToTargetEdge.get_targetNodeDirection() * -1)
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(sourceToTargetEdge,
                                                                                reverseTargetToSourceEdge)
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        targetToThirdEdge, reverseThirdToTargetEdge = graph.add_edges_to_graph(rc_targetToThirdEdge,
                                                                            rc_reverseThirdToTargetEdge)
        targetToThirdEdge.increment_edge_coverage()
        reverseThirdToTargetEdge.increment_edge_coverage()
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
                            3)
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
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        graph.add_edge_to_edges(targetToThirdEdge)
        graph.add_edge_to_edges(reverseThirdToTargetEdge)
        targetToThirdEdge.increment_edge_coverage()
        reverseThirdToTargetEdge.increment_edge_coverage()
        # execution
        rc_sourceToTargetEdge = sourceToTargetEdge
        rc_sourceToTargetEdge.set_sourceNodeDirection(sourceToTargetEdge.get_sourceNodeDirection() * -1)
        rc_sourceToTargetEdge.set_targetNodeDirection(sourceToTargetEdge.get_targetNodeDirection() * -1)
        rc_reverseTargetToSourceEdge = reverseTargetToSourceEdge
        rc_reverseTargetToSourceEdge.set_sourceNodeDirection(reverseTargetToSourceEdge.get_sourceNodeDirection() * -1)
        rc_reverseTargetToSourceEdge.set_targetNodeDirection(reverseTargetToSourceEdge.get_targetNodeDirection() * -1)
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(rc_sourceToTargetEdge,
                                                                                rc_reverseTargetToSourceEdge)
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        targetToThirdEdge, reverseThirdToTargetEdge = graph.add_edges_to_graph(targetToThirdEdge,
                                                                            reverseThirdToTargetEdge)
        targetToThirdEdge.increment_edge_coverage()
        reverseThirdToTargetEdge.increment_edge_coverage()
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

    def test___add_edge_to_node_forward_empty(self):
        # setup
        class FakeEdge:
            def __init__(self,
                        hash,
                        sourceNodeDirection):
                self.hash = hash
                self.sourceNodeDirection = sourceNodeDirection
            def __hash__(self):
                return self.hash
            def get_sourceNodeDirection(self):
                return self.sourceNodeDirection
        mock_edge = FakeEdge(12345,
                            +1)
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes)
        node = [Node(x) for x in read1.get_geneMers(3)][0]
        graph = GeneMerGraph({},
                            3)
        # execution
        actual_updated_node = graph.add_edge_to_node(node,
                                                    mock_edge)
        actual_forward_edges = actual_updated_node.get_forward_edge_hashes()
        actual_backward_edges = actual_updated_node.get_backward_edge_hashes()
        # assertion
        expected_forward_edges = [12345]
        expected_backward_edges = []
        self.assertEqual(actual_updated_node, node)
        self.assertEqual(actual_forward_edges, expected_forward_edges)
        self.assertEqual(actual_backward_edges, expected_backward_edges)

    def test___add_edge_to_node_forward_not_empty(self):
        # setup
        class FakeEdge:
            def __init__(self,
                        hash,
                        sourceNodeDirection):
                self.hash = hash
                self.sourceNodeDirection = sourceNodeDirection
            def __hash__(self):
                return self.hash
            def get_sourceNodeDirection(self):
                return self.sourceNodeDirection
        mock_edge = FakeEdge(12345,
                            +1)
        mock_edge2 = FakeEdge(56789,
                            +1)
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes)
        node = [Node(x) for x in read1.get_geneMers(3)][0]
        graph = GeneMerGraph({},
                            3)
        graph.add_edge_to_node(node,
                            mock_edge)
        # execution
        actual_updated_node = graph.add_edge_to_node(node,
                                                    mock_edge2)
        actual_forward_edges = actual_updated_node.get_forward_edge_hashes()
        actual_backward_edges = actual_updated_node.get_backward_edge_hashes()
        # assertion
        expected_forward_edges = [12345, 56789]
        expected_backward_edges = []
        self.assertEqual(actual_updated_node, node)
        self.assertEqual(actual_forward_edges, expected_forward_edges)
        self.assertEqual(actual_backward_edges, expected_backward_edges)

    def test___add_edge_to_node_backward_empty(self):
        # setup
        class FakeEdge:
            def __init__(self,
                        hash,
                        sourceNodeDirection):
                self.hash = hash
                self.sourceNodeDirection = sourceNodeDirection
            def __hash__(self):
                return self.hash
            def get_sourceNodeDirection(self):
                return self.sourceNodeDirection
        mock_edge = FakeEdge(12345,
                            -1)
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes)
        node = [Node(x) for x in read1.get_geneMers(3)][0]
        graph = GeneMerGraph({},
                            3)
        # execution
        actual_updated_node = graph.add_edge_to_node(node,
                                                    mock_edge)
        actual_forward_edges = actual_updated_node.get_forward_edge_hashes()
        actual_backward_edges = actual_updated_node.get_backward_edge_hashes()
        # assertion
        expected_forward_edges = []
        expected_backward_edges = [12345]
        self.assertEqual(actual_updated_node, node)
        self.assertEqual(actual_forward_edges, expected_forward_edges)
        self.assertEqual(actual_backward_edges, expected_backward_edges)

    def test___add_edge_to_node_backward_not_empty(self):
        # setup
        class FakeEdge:
            def __init__(self,
                        hash,
                        sourceNodeDirection):
                self.hash = hash
                self.sourceNodeDirection = sourceNodeDirection
            def __hash__(self):
                return self.hash
            def get_sourceNodeDirection(self):
                return self.sourceNodeDirection
        mock_edge = FakeEdge(12345,
                            -1)
        mock_edge2 = FakeEdge(56789,
                            -1)
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes)
        node = [Node(x) for x in read1.get_geneMers(3)][0]
        graph = GeneMerGraph({},
                            3)
        graph.add_edge_to_node(node,
                            mock_edge)
        # execution
        actual_updated_node = graph.add_edge_to_node(node,
                                                    mock_edge2)
        actual_forward_edges = actual_updated_node.get_forward_edge_hashes()
        actual_backward_edges = actual_updated_node.get_backward_edge_hashes()
        # assertion
        expected_forward_edges = []
        expected_backward_edges = [12345, 56789]
        self.assertEqual(actual_updated_node, node)
        self.assertEqual(actual_forward_edges, expected_forward_edges)
        self.assertEqual(actual_backward_edges, expected_backward_edges)

    def test___get_degree_1_edge_2_nodes(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        graph = GeneMerGraph({},
                            3)
        graph.add_edge(sourceGeneMer,
                    targetGeneMer)
        sourceNode = graph.get_node(sourceGeneMer)
        targetNode = graph.get_node(targetGeneMer)
        # execution
        actual_source_degrees = graph.get_degree(sourceNode)
        actual_target_degrees = graph.get_degree(targetNode)
        # assertion
        expected_source_degrees = 1
        expected_target_degrees = 1
        self.assertEqual(actual_source_degrees, expected_source_degrees)
        self.assertEqual(actual_target_degrees, expected_target_degrees)

    def test___get_degree_2_edge_3_nodes(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        thirdGeneMer = geneMers[2]
        graph = GeneMerGraph({},
                            3)
        graph.add_edge(sourceGeneMer,
                    targetGeneMer)
        graph.add_edge(targetGeneMer,
                    thirdGeneMer)
        sourceNode = graph.get_node(sourceGeneMer)
        targetNode = graph.get_node(targetGeneMer)
        thirdNode = graph.get_node(thirdGeneMer)
        # execution
        actual_source_degrees = graph.get_degree(sourceNode)
        actual_target_degrees = graph.get_degree(targetNode)
        actual_third_degrees = graph.get_degree(thirdNode)
        # assertion
        expected_source_degrees = 1
        expected_target_degrees = 2
        expected_third_degrees = 1
        self.assertEqual(actual_source_degrees, expected_source_degrees)
        self.assertEqual(actual_target_degrees, expected_target_degrees)
        self.assertEqual(actual_third_degrees, expected_third_degrees)

    def test___get_degree_3_edge_4_nodes(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        thirdGeneMer = geneMers[2]
        fourthGeneMer = geneMers[3]
        graph = GeneMerGraph({},
                            3)
        graph.add_edge(sourceGeneMer,
                    targetGeneMer)
        graph.add_edge(targetGeneMer,
                    thirdGeneMer)
        graph.add_edge(thirdGeneMer,
                    fourthGeneMer)
        sourceNode = graph.get_node(sourceGeneMer)
        targetNode = graph.get_node(targetGeneMer)
        thirdNode = graph.get_node(thirdGeneMer)
        fourthNode = graph.get_node(fourthGeneMer)
        # execution
        actual_source_degrees = graph.get_degree(sourceNode)
        actual_target_degrees = graph.get_degree(targetNode)
        actual_third_degrees = graph.get_degree(thirdNode)
        actual_fourth_degrees = graph.get_degree(fourthNode)
        # assertion
        expected_source_degrees = 1
        expected_target_degrees = 2
        expected_third_degrees = 2
        expected_fourth_degrees = 1
        self.assertEqual(actual_source_degrees, expected_source_degrees)
        self.assertEqual(actual_target_degrees, expected_target_degrees)
        self.assertEqual(actual_third_degrees, expected_third_degrees)
        self.assertEqual(actual_fourth_degrees, expected_fourth_degrees)

    def test___get_degree_3_edge_3_nodes_first(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene6"]
        read1 = Read("read1",
                    genes)
        read2 = Read("read1",
                    genes2)
        geneMers = [x for x in read1.get_geneMers(3)]
        geneMers2 = [x for x in read2.get_geneMers(3)]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        targetGeneMer2 = geneMers2[1]
        graph = GeneMerGraph({},
                            3)
        graph.add_edge(sourceGeneMer,
                    targetGeneMer)
        graph.add_edge(sourceGeneMer,
                    targetGeneMer2)
        sourceNode = graph.get_node(sourceGeneMer)
        targetNode = graph.get_node(targetGeneMer)
        targetNode2 = graph.get_node(targetGeneMer2)
        # execution
        actual_source_degrees = graph.get_degree(sourceNode)
        actual_target_degrees = graph.get_degree(targetNode)
        actual_target2_degrees = graph.get_degree(targetNode2)
        # assertion
        expected_source_degrees = 2
        expected_target_degrees = 1
        expected_target2_degrees = 1
        self.assertEqual(actual_source_degrees, expected_source_degrees)
        self.assertEqual(actual_target_degrees, expected_target_degrees)
        self.assertEqual(actual_target2_degrees, expected_target2_degrees)

    def test___get_degree_3_edge_3_nodes_last(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        genes2 = ["+gene0", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1",
                    genes)
        read2 = Read("read1",
                    genes2)
        geneMers = [x for x in read1.get_geneMers(3)]
        geneMers2 = [x for x in read2.get_geneMers(3)]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        sourceGeneMer2 = geneMers2[0]
        graph = GeneMerGraph({},
                            3)
        graph.add_edge(sourceGeneMer,
                    targetGeneMer)
        graph.add_edge(sourceGeneMer2,
                    targetGeneMer)
        sourceNode = graph.get_node(sourceGeneMer)
        targetNode = graph.get_node(targetGeneMer)
        sourceNode2 = graph.get_node(sourceGeneMer2)
        # execution
        actual_source_degrees = graph.get_degree(sourceNode)
        actual_target_degrees = graph.get_degree(targetNode)
        actual_source2_degrees = graph.get_degree(sourceNode2)
        # assertion
        expected_source_degrees = 1
        expected_target_degrees = 2
        expected_source2_degrees = 1
        self.assertEqual(actual_source_degrees, expected_source_degrees)
        self.assertEqual(actual_target_degrees, expected_target_degrees)
        self.assertEqual(actual_source2_degrees, expected_source2_degrees)

    def test___get_degree_5_edge_5_nodes_middle(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene0", "-gene2", "+gene3", "-gene4", "+gene6"]
        read1 = Read("read1",
                    genes)
        read2 = Read("read1",
                    genes2)
        geneMers = [x for x in read1.get_geneMers(3)]
        geneMers2 = [x for x in read2.get_geneMers(3)]
        firstGeneMer = geneMers[0]
        middleGeneMer = geneMers[1]
        thirdGeneMer = geneMers[2]
        firstGeneMer2 = geneMers2[0]
        middleGeneMer2 = geneMers2[1]
        thirdGeneMer2 = geneMers2[2]
        graph = GeneMerGraph({},
                            3)
        graph.add_edge(firstGeneMer,
                    middleGeneMer)
        graph.add_edge(middleGeneMer,
                    thirdGeneMer)
        graph.add_edge(firstGeneMer2,
                    middleGeneMer2)
        graph.add_edge(middleGeneMer2,
                    thirdGeneMer2)
        firstNode = graph.get_node(firstGeneMer)
        middleNode = graph.get_node(middleGeneMer)
        thirdNode = graph.get_node(thirdGeneMer)
        firstNode2 = graph.get_node(firstGeneMer2)
        thirdNode2 = graph.get_node(thirdGeneMer2)
        # execution
        actual_first_degrees = graph.get_degree(firstNode)
        actual_middle_degrees = graph.get_degree(middleNode)
        actual_third_degrees = graph.get_degree(thirdNode)
        actual_first_degrees2 = graph.get_degree(firstNode2)
        actual_third_degrees2 = graph.get_degree(thirdNode2)
        # assertion
        expected_first_degrees = 1
        expected_middle_degrees = 4
        expected_third_degrees = 1
        expected_first_degrees2 = 1
        expected_third_degrees2 = 1
        self.assertEqual(actual_first_degrees, expected_first_degrees)
        self.assertEqual(actual_middle_degrees, expected_middle_degrees)
        self.assertEqual(actual_third_degrees, expected_third_degrees)
        self.assertEqual(actual_first_degrees2, expected_first_degrees2)
        self.assertEqual(actual_third_degrees2, expected_third_degrees2)

    def test___get_degree_1_edge_2_node_duplicate(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        graph = GeneMerGraph({},
                            3)
        graph.add_edge(sourceGeneMer,
                    targetGeneMer)
        graph.add_edge(sourceGeneMer,
                    targetGeneMer)
        sourceNode = graph.get_node(sourceGeneMer)
        targetNode = graph.get_node(targetGeneMer)
        # execution
        actual_source_degrees = graph.get_degree(sourceNode)
        actual_target_degrees = graph.get_degree(targetNode)
        # assertion
        expected_source_degrees = 1
        expected_target_degrees = 1
        self.assertEqual(actual_source_degrees, expected_source_degrees)
        self.assertEqual(actual_target_degrees, expected_target_degrees)

    def test___remove_existing_edge(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        graph = GeneMerGraph({},
                            3)
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(sourceGeneMer,
                                                                    targetGeneMer)
        sourceNode = sourceToTargetEdge.get_sourceNode()
        targetNode = reverseTargetToSourceEdge.get_sourceNode()
        sourceNodeEdgeHashes = sourceNode.get_forward_edge_hashes() + sourceNode.get_backward_edge_hashes()
        targetNodeEdgeHashes = targetNode.get_forward_edge_hashes() + targetNode.get_backward_edge_hashes()
        # sanity checks
        self.assertNotEqual(sourceNodeEdgeHashes, [])
        self.assertNotEqual(targetNodeEdgeHashes, [])
        self.assertNotEqual(graph.get_edges(), {})
        # execution
        for edgeHash in sourceNodeEdgeHashes:
            graph.remove_edge(edgeHash)
        for edgeHash in targetNodeEdgeHashes:
            graph.remove_edge(edgeHash)
        actual_sourceNodeEdgeHashes = sourceNode.get_forward_edge_hashes() + sourceNode.get_backward_edge_hashes()
        actual_targetNodeEdgeHashes = targetNode.get_forward_edge_hashes() + targetNode.get_backward_edge_hashes()
        actual_edges = graph.get_edges()
        # assertion
        expected_sourceNodeEdgeHashes = []
        expected_targetNodeEdgeHashes = []
        expected_edges = {}
        self.assertEqual(actual_sourceNodeEdgeHashes, expected_sourceNodeEdgeHashes)
        self.assertEqual(actual_targetNodeEdgeHashes, expected_targetNodeEdgeHashes)
        self.assertEqual(actual_edges, expected_edges)

    def test___remove_non_existing_edge(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        graph = GeneMerGraph({},
                            3)
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(sourceGeneMer,
                                                                    targetGeneMer)
        sourceNode = sourceToTargetEdge.get_sourceNode()
        targetNode = reverseTargetToSourceEdge.get_sourceNode()
        sourceNodeEdgeHashes = sourceNode.get_forward_edge_hashes() + sourceNode.get_backward_edge_hashes()
        targetNodeEdgeHashes = targetNode.get_forward_edge_hashes() + targetNode.get_backward_edge_hashes()
        # sanity checks
        self.assertNotEqual(sourceNodeEdgeHashes, [])
        self.assertNotEqual(targetNodeEdgeHashes, [])
        self.assertNotEqual(graph.get_edges(), {})
        # execution
        self.assertRaises(AssertionError, graph.remove_edge, 12345)

    def test___assign_Id_to_nodes(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        graph = GeneMerGraph({},
                            3)
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(sourceGeneMer,
                                                                    targetGeneMer)
        # execution
        graph.assign_Id_to_nodes()
        # assertion
        expected_sourceNode_Id = 0
        expected_targetNode_Id = 1
        self.assertEqual(sourceToTargetEdge.get_sourceNode().get_nodeId(), expected_sourceNode_Id)
        self.assertEqual(reverseTargetToSourceEdge.get_sourceNode().get_nodeId(), expected_targetNode_Id)

    def test___write_gml_to_file(self):
        # setup
        import os
        graph = GeneMerGraph({},
                            3)
        # execution
        graph.write_gml_to_file("tests/test_graph",
                            ["graph\t[", "\t\tnode\t[]", "\t\tedge\t[]", "]"])
        # assertion
        self.assertTrue(os.path.exists("tests/test_graph.gml"))
        with open("tests/test_graph.gml", "r") as i:
            actual_content = i.read().splitlines()
        self.assertEqual(actual_content, ["graph\t[", "\t\tnode\t[]", "\t\tedge\t[]", "]"])
        os.remove("tests/test_graph.gml")

    def test___write_gml_to_file_in_subdirectory(self):
        # setup
        import os
        graph = GeneMerGraph({},
                            3)
        # execution
        graph.write_gml_to_file("tests/tests/test_graph",
                                ["graph\t[", "\t\tnode\t[]", "\t\tedge\t[]", "]"])
        # assertion
        self.assertTrue(os.path.exists("tests/tests/test_graph.gml"))
        with open("tests/tests/test_graph.gml", "r") as i:
            actual_content = i.read().splitlines()
        self.assertEqual(actual_content, ["graph\t[", "\t\tnode\t[]", "\t\tedge\t[]", "]"])
        import shutil
        shutil.rmtree("tests/tests")

    def test___generate_gml(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        graph = GeneMerGraph({},
                            3)
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(sourceGeneMer,
                                                                    targetGeneMer)
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        sourceToTargetEdge.get_sourceNode().increment_node_coverage()
        reverseTargetToSourceEdge.get_sourceNode().increment_node_coverage()
        # execution
        actual_writtenGraph = graph.generate_gml("tests/test_graph",
                                                3,
                                                1,
                                                1)
        # assertion
        expected_writtenGraph = [['graph\t[',
                                '\tnode\t[\n\t\tid\t0\n\t\tlabel\t"-gene3~~~+gene2~~~-gene1"\n\t\tcoverage\t1\n\t\treads\t""\n\t]',
                                '\tedge\t[\n\t\tsource\t0\n\t\ttarget\t1\n\t\tweight\t1\n\t]',
                                '\tnode\t[\n\t\tid\t1\n\t\tlabel\t"+gene4~~~-gene3~~~+gene2"\n\t\tcoverage\t1\n\t\treads\t""\n\t]',
                                '\tedge\t[\n\t\tsource\t1\n\t\ttarget\t0\n\t\tweight\t1\n\t]',
                                ']'],
                                ['graph\t[',
                                '\tnode\t[\n\t\tid\t0\n\t\tlabel\t"+gene1~~~-gene2~~~+gene3"\n\t\tcoverage\t1\n\t\treads\t""\n\t]',
                                '\tedge\t[\n\t\tsource\t0\n\t\ttarget\t1\n\t\tweight\t1\n\t]',
                                '\tnode\t[\n\t\tid\t1\n\t\tlabel\t"-gene2~~~+gene3~~~-gene4"\n\t\tcoverage\t1\n\t\treads\t""\n\t]',
                                '\tedge\t[\n\t\tsource\t1\n\t\ttarget\t0\n\t\tweight\t1\n\t]',
                                ']'],
                                ['graph\t[',
                                '\tnode\t[\n\t\tid\t0\n\t\tlabel\t"-gene3~~~+gene2~~~-gene1"\n\t\tcoverage\t1\n\t\treads\t""\n\t]',
                                '\tedge\t[\n\t\tsource\t0\n\t\ttarget\t1\n\t\tweight\t1\n\t]',
                                '\tnode\t[\n\t\tid\t1\n\t\tlabel\t"-gene2~~~+gene3~~~-gene4"\n\t\tcoverage\t1\n\t\treads\t""\n\t]',
                                '\tedge\t[\n\t\tsource\t1\n\t\ttarget\t0\n\t\tweight\t1\n\t]',
                                ']'],
                                ['graph\t[',
                                '\tnode\t[\n\t\tid\t0\n\t\tlabel\t"+gene1~~~-gene2~~~+gene3"\n\t\tcoverage\t1\n\t\treads\t""\n\t]',
                                '\tedge\t[\n\t\tsource\t0\n\t\ttarget\t1\n\t\tweight\t1\n\t]',
                                '\tnode\t[\n\t\tid\t1\n\t\tlabel\t"+gene4~~~-gene3~~~+gene2"\n\t\tcoverage\t1\n\t\treads\t""\n\t]',
                                '\tedge\t[\n\t\tsource\t1\n\t\ttarget\t0\n\t\tweight\t1\n\t]',
                                ']']]
        import os
        self.assertTrue(os.path.exists("tests/test_graph.3.1.1.gml"))
        self.assertTrue(any(actual_writtenGraph == e for e in expected_writtenGraph))
        os.remove("tests/test_graph.3.1.1.gml")

    def test___get_gene_mer_label(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes)
        sourceGeneMer = [x for x in read1.get_geneMers(3)][0]
        sourceNode = Node(sourceGeneMer)
        graph = GeneMerGraph({},
                            3)
        # execution
        actual_geneMerString = graph.get_gene_mer_label(sourceNode)
        # assertion
        expected_geneMerString = ["+gene1~~~-gene2~~~+gene3", "-gene3~~~+gene2~~~-gene1"]
        self.assertTrue(any(actual_geneMerString == e for e in expected_geneMerString))

    def test___filter_graph(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2},
                            3)
        # execution
        graph.filter_graph(2,
                        2)
        actual_writtenGraph = graph.generate_gml("tests/test_graph",
                                                3,
                                                2,
                                                2)
        actual_numberOfNodes = graph.get_total_number_of_nodes()
        actual_numberOfEdges = graph.get_total_number_of_edges()
        # assertion
        expected_numberOfNodes = 6
        expected_numberOfEdges = 10
        import os
        self.assertEqual(actual_numberOfNodes, expected_numberOfNodes)
        self.assertEqual(actual_numberOfEdges, expected_numberOfEdges)
        os.remove("tests/test_graph.3.2.2.gml")

    def test___filter_all_graph(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2},
                            3)
        # execution
        graph.filter_graph(10,
                        10)
        actual_writtenGraph = graph.generate_gml("tests/test_graph",
                                                3,
                                                10,
                                                10)
        actual_numberOfNodes = graph.get_total_number_of_nodes()
        actual_numberOfEdges = graph.get_total_number_of_edges()
        # assertion
        expected_writtenGraph = ['graph\t[', ']']
        expected_numberOfNodes = 0
        expected_numberOfEdges = 0
        import os
        self.assertEqual(actual_writtenGraph, expected_writtenGraph)
        self.assertEqual(actual_numberOfNodes, expected_numberOfNodes)
        self.assertEqual(actual_numberOfEdges, expected_numberOfEdges)
        os.remove("tests/test_graph.3.10.10.gml")

    def test___add_node_to_read(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes)
        sourceGeneMer = [x for x in read1.get_geneMers(3)][0]
        sourceNode = Node(sourceGeneMer)
        graph = GeneMerGraph({},
                            3)
        # execution
        actual_listOfNodes = graph.add_node_to_read(sourceNode,
                                                    "read1")
        actual_readNodeDict = graph.get_readNodes()
        # assertion
        expected_listOfNodes = [sourceNode.__hash__()]
        expected_readNodeDict = {"read1": [sourceNode.__hash__()]}
        self.assertEqual(actual_listOfNodes, expected_listOfNodes)
        self.assertEqual(actual_readNodeDict, expected_readNodeDict)

    def test___get_nodes_containing_read_filtered_graph(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        graph = GeneMerGraph({"read1": genes1},
                            3)
        graph.filter_graph(2,
                        2)
        # execution
        actual_listOfNodes = graph.get_nodes_containing_read("read1")
        # assertion
        self.assertEqual(len(actual_listOfNodes), 2)

    def test___remove_node_from_reads_one_copy(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5"]
        read1 = Read("read1",
                    genes)
        graph = GeneMerGraph({},
                            3)
        sourceNodes = []
        for s in read1.get_geneMers(3):
            sourceNode = graph.add_node(s, [read1])
            graph.add_node_to_read(sourceNode,
                                read1.get_readId())
            sourceNodes.append(sourceNode)
        # execution
        test_node = sourceNodes[1]
        graph.remove_node_from_reads(test_node)
        # assertion
        expectedReadNodes = {"read1": [
                    sourceNodes[0].__hash__(),
                    sourceNodes[2].__hash__()
                                    ]}
        self.assertEqual(graph.get_readNodes(), expectedReadNodes)

    def test___remove_node_from_reads_more_than_one_copy(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1",
                    genes)
        graph = GeneMerGraph({},
                            3)
        sourceNodes = []
        for s in read1.get_geneMers(3):
            sourceNode = graph.add_node(s, [read1])
            graph.add_node_to_read(sourceNode,
                                read1.get_readId())
            sourceNodes.append(sourceNode)
        # execution
        test_node = sourceNodes[1]
        graph.remove_node_from_reads(test_node)
        # assertion
        expectedReadNodes = {"read1": [
                    sourceNodes[0].__hash__(),
                    sourceNodes[2].__hash__(),
                    sourceNodes[3].__hash__(),
                    sourceNodes[4].__hash__()
                                    ]}
        self.assertEqual(graph.get_readNodes(), expectedReadNodes)

    def test___get_existing_forward_node_from_node(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g,
                                [read1])
            nodes.append(node)
        mock_forward_edge = Edge(nodes[0],
                                nodes[1],
                                1,
                                1)
        mock_rc_forward_edge = Edge(nodes[1],
                                    nodes[0],
                                    -1,
                                    -1)
        graph.add_edges_to_graph(mock_forward_edge,
                                mock_rc_forward_edge)
        graph.add_edge_to_node(nodes[0],
                            mock_forward_edge)
        graph.add_edge_to_node(nodes[1],
                            mock_rc_forward_edge)
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = graph.get_forward_node_from_node(nodes[0])
        actual_targetNode = actual_targetNode
        actual_targetDirection = actual_targetDirection
        # assertion
        expected_targetNode = nodes[1]
        expected_targetDirection = 1
        expected_extend = True
        self.assertEqual(actual_targetNode, expected_targetNode)
        self.assertEqual(actual_targetDirection, expected_targetDirection)
        self.assertEqual(actual_extend, expected_extend)

    def test___get_existing_forward_node_from_node_in_middle(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g,
                                [read1])
            nodes.append(node)
        for n in range(len(nodes) - 1):
            mock_forward_edge = Edge(nodes[n],
                                    nodes[n+1],
                                    1,
                                    1)
            mock_rc_forward_edge = Edge(nodes[n+1],
                                        nodes[n],
                                        -1,
                                        -1)
            graph.add_edges_to_graph(mock_forward_edge,
                                    mock_rc_forward_edge)
            graph.add_edge_to_node(nodes[n],
                                mock_forward_edge)
            graph.add_edge_to_node(nodes[n+1],
                                mock_rc_forward_edge)
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = graph.get_forward_node_from_node(nodes[0])
        actual_targetNode = actual_targetNode
        actual_targetDirection = actual_targetDirection
        # assertion
        expected_targetNode = nodes[1]
        expected_targetDirection = 1
        expected_extend = True
        self.assertEqual(actual_targetNode, expected_targetNode)
        self.assertEqual(actual_targetDirection, expected_targetDirection)
        self.assertEqual(actual_extend, expected_extend)

    def test___get_branched_forward_node_from_node(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene7", "-gene8"]
        read1 = Read("read1",
                    genes1)
        read2 =  Read("read2",
                    genes2)
        geneMers1 = [x for x in read1.get_geneMers(3)]
        geneMers2 = [x for x in read2.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes1 = []
        for g in geneMers1:
            node = graph.add_node(g,
                                [read1])
            nodes1.append(node)
        nodes2 = []
        for g in geneMers2:
            node = graph.add_node(g,
                                [read2])
            nodes2.append(node)
        for n in [nodes1, nodes2]:
            mock_forward_edge = Edge(n[1],
                                    n[2],
                                    1,
                                    1)
            mock_rc_forward_edge = Edge(n[2],
                                        n[1],
                                        -1,
                                        -1)
            graph.add_edges_to_graph(mock_forward_edge,
                                    mock_rc_forward_edge)
            graph.add_edge_to_node(n[1],
                                mock_forward_edge)
            graph.add_edge_to_node(n[2],
                                mock_rc_forward_edge)
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = graph.get_forward_node_from_node(nodes1[0])
        actual_targetNode = actual_targetNode
        actual_targetDirection = actual_targetDirection
        # assertion
        expected_targetNode = None
        expected_targetDirection = None
        expected_extend = False
        self.assertEqual(actual_targetNode, expected_targetNode)
        self.assertEqual(actual_targetDirection, expected_targetDirection)
        self.assertEqual(actual_extend, expected_extend)

    def test___get_non_existing_forward_node_from_node(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g,
                                [read1])
            nodes.append(node)
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = graph.get_forward_node_from_node(nodes[0])
        # assertion
        expected_targetNode = None
        expected_targetDirection = None
        expected_extend = False
        self.assertEqual(actual_targetNode, expected_targetNode)
        self.assertEqual(actual_targetDirection, expected_targetDirection)
        self.assertEqual(actual_extend, expected_extend)

    def test___get_existing_backward_node_from_node(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g,
                                [read1])
            nodes.append(node)
        mock_backward_edge = Edge(nodes[0],
                                nodes[1],
                                -1,
                                -1)
        mock_rc_backward_edge = Edge(nodes[1],
                                    nodes[0],
                                    1,
                                    1)
        graph.add_edges_to_graph(mock_backward_edge,
                                mock_rc_backward_edge)
        graph.add_edge_to_node(nodes[0],
                            mock_backward_edge)
        graph.add_edge_to_node(nodes[1],
                            mock_rc_backward_edge)
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = graph.get_backward_node_from_node(nodes[0])
        actual_targetNode = actual_targetNode
        actual_targetDirection = actual_targetDirection
        # assertion
        expected_targetNode = nodes[1]
        expected_targetDirection = -1
        expected_extend = True
        self.assertEqual(actual_targetNode, expected_targetNode)
        self.assertEqual(actual_targetDirection, expected_targetDirection)
        self.assertEqual(actual_extend, expected_extend)

    def test___get_non_existing_backward_node_from_node(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g,
                                [read1])
            nodes.append(node)
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = graph.get_backward_node_from_node(nodes[0])
        # assertion
        expected_targetNode = None
        expected_targetDirection = None
        expected_extend = False
        self.assertEqual(actual_targetNode, expected_targetNode)
        self.assertEqual(actual_targetDirection, expected_targetDirection)
        self.assertEqual(actual_extend, expected_extend)

    def test___get_existing_backward_node_from_node_in_middle(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g,
                                [read1])
            nodes.append(node)
        for n in range(len(nodes) - 1):
            mock_backward_edge = Edge(nodes[n],
                                    nodes[n+1],
                                    -1,
                                    -1)
            mock_rc_backward_edge = Edge(nodes[n+1],
                                        nodes[n],
                                        1,
                                        1)
            graph.add_edges_to_graph(mock_backward_edge,
                                    mock_rc_backward_edge)
            graph.add_edge_to_node(nodes[n],
                                mock_backward_edge)
            graph.add_edge_to_node(nodes[n+1],
                                mock_rc_backward_edge)
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = graph.get_backward_node_from_node(nodes[0])
        actual_targetNode = actual_targetNode
        actual_targetDirection = actual_targetDirection
        # assertion
        expected_targetNode = nodes[1]
        expected_targetDirection = -1
        expected_extend = True
        self.assertEqual(actual_targetNode, expected_targetNode)
        self.assertEqual(actual_targetDirection, expected_targetDirection)
        self.assertEqual(actual_extend, expected_extend)

    def test___get_branched_backward_node_from_node(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene7", "-gene8"]
        read1 = Read("read1",
                    genes1)
        read2 =  Read("read2",
                    genes2)
        geneMers1 = [x for x in read1.get_geneMers(3)]
        geneMers2 = [x for x in read2.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes1 = []
        for g in geneMers1:
            node = graph.add_node(g,
                                [read1])
            nodes1.append(node)
        nodes2 = []
        for g in geneMers2:
            node = graph.add_node(g,
                                [read2])
            nodes2.append(node)
        for n in [nodes1, nodes2]:
            mock_backward_edge = Edge(n[1],
                                    n[2],
                                    -1,
                                    -1)
            mock_rc_backward_edge = Edge(n[2],
                                        n[1],
                                        1,
                                        1)
            graph.add_edges_to_graph(mock_backward_edge,
                                    mock_rc_backward_edge)
            graph.add_edge_to_node(n[1],
                                mock_backward_edge)
            graph.add_edge_to_node(n[2],
                                mock_rc_backward_edge)
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = graph.get_forward_node_from_node(nodes1[0])
        actual_targetNode = actual_targetNode
        actual_targetDirection = actual_targetDirection
        # assertion
        expected_targetNode = None
        expected_targetDirection = None
        expected_extend = False
        self.assertEqual(actual_targetNode, expected_targetNode)
        self.assertEqual(actual_targetDirection, expected_targetDirection)
        self.assertEqual(actual_extend, expected_extend)

    def test___get_forward_path_from_node_linear(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene7", "-gene8"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g,
                                [read1])
            nodes.append(node)
        for n in range(len(nodes) - 1):
            mock_forward_edge = Edge(nodes[n],
                                    nodes[n+1],
                                    1,
                                    1)
            mock_rc_forward_edge = Edge(nodes[n+1],
                                        nodes[n],
                                        -1,
                                        -1)
            graph.add_edges_to_graph(mock_forward_edge,
                                    mock_rc_forward_edge)
            graph.add_edge_to_node(nodes[n],
                                mock_forward_edge)
            graph.add_edge_to_node(nodes[n+1],
                                mock_rc_forward_edge)
        # execution
        actual_forward_path_from_node = graph.get_forward_path_from_node(nodes[1])
        actual_path_length = len(actual_forward_path_from_node)
        # assertion
        expected_forward_path_from_node = [n.__hash__() for n in nodes[1:]]
        expected_path_length = 5
        self.assertEqual(actual_forward_path_from_node, expected_forward_path_from_node)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___get_forward_path_from_node_branched(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene7"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene8", "+gene8"]
        read1 = Read("read1",
                    genes1)
        read2 =  Read("read2",
                    genes2)
        geneMers1 = [x for x in read1.get_geneMers(3)]
        geneMers2 = [x for x in read2.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes1 = []
        for g in geneMers1:
            node = graph.add_node(g,
                                [read1])
            nodes1.append(node)
        nodes2 = []
        for g in geneMers2:
            node = graph.add_node(g,
                                [read2])
            nodes2.append(node)
        for nodes in [nodes1, nodes2]:
            for n in range(len(nodes) - 1):
                mock_forward_edge = Edge(nodes[n],
                                        nodes[n+1],
                                        1,
                                        1)
                mock_rc_forward_edge = Edge(nodes[n+1],
                                            nodes[n],
                                            -1,
                                            -1)
                graph.add_edges_to_graph(mock_forward_edge,
                                        mock_rc_forward_edge)
                graph.add_edge_to_node(nodes[n],
                                    mock_forward_edge)
                graph.add_edge_to_node(nodes[n+1],
                                    mock_rc_forward_edge)
        # execution
        actual_forward_path_from_node = graph.get_forward_path_from_node(nodes1[0])
        actual_path_length = len(actual_forward_path_from_node)
        # assertion
        expected_forward_path_from_node = [n.__hash__() for n in nodes1[:2]]
        expected_path_length = 2
        self.assertEqual(actual_forward_path_from_node, expected_forward_path_from_node)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___get_forward_path_from_node_branched_want_branched(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene7"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene8", "+gene8"]
        read1 = Read("read1",
                    genes1)
        read2 =  Read("read2",
                    genes2)
        geneMers1 = [x for x in read1.get_geneMers(3)]
        geneMers2 = [x for x in read2.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes1 = []
        for g in geneMers1:
            node = graph.add_node(g,
                                [read1])
            nodes1.append(node)
        nodes2 = []
        for g in geneMers2:
            node = graph.add_node(g,
                                [read2])
            nodes2.append(node)
        for nodes in [nodes1, nodes2]:
            for n in range(len(nodes) - 1):
                mock_forward_edge = Edge(nodes[n],
                                        nodes[n+1],
                                        1,
                                        1)
                mock_rc_forward_edge = Edge(nodes[n+1],
                                            nodes[n],
                                            -1,
                                            -1)
                graph.add_edges_to_graph(mock_forward_edge,
                                        mock_rc_forward_edge)
                graph.add_edge_to_node(nodes[n],
                                    mock_forward_edge)
                graph.add_edge_to_node(nodes[n+1],
                                    mock_rc_forward_edge)
        # execution
        actual_forward_path_from_node = graph.get_forward_path_from_node(nodes1[0],
                                                                        True)
        actual_path_length = len(actual_forward_path_from_node)
        # assertion
        expected_forward_path_from_node = [n.__hash__() for n in nodes1[:3]]
        expected_path_length = 3
        self.assertEqual(actual_forward_path_from_node, expected_forward_path_from_node)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___get_backward_path_from_node_linear(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene7", "-gene8"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g,
                                [read1])
            nodes.append(node)
        for n in range(len(nodes) - 1):
            mock_backward_edge = Edge(nodes[n],
                                    nodes[n+1],
                                    -1,
                                    -1)
            mock_rc_backward_edge = Edge(nodes[n+1],
                                        nodes[n],
                                        1,
                                        1)
            graph.add_edges_to_graph(mock_backward_edge,
                                    mock_rc_backward_edge)
            graph.add_edge_to_node(nodes[n],
                                mock_backward_edge)
            graph.add_edge_to_node(nodes[n+1],
                                mock_rc_backward_edge)
        # execution
        actual_backward_path_from_node = graph.get_backward_path_from_node(nodes[1])
        actual_path_length = len(actual_backward_path_from_node)
        # assertion
        expected_backward_path_from_node = list(reversed([n.__hash__() for n in nodes[1:]]))
        expected_path_length = 5
        self.assertEqual(actual_backward_path_from_node, expected_backward_path_from_node)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___get_backward_path_from_node_branched_want_branched(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene7"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene8", "+gene8"]
        read1 = Read("read1",
                    genes1)
        read2 =  Read("read2",
                    genes2)
        geneMers1 = [x for x in read1.get_geneMers(3)]
        geneMers2 = [x for x in read2.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes1 = []
        for g in geneMers1:
            node = graph.add_node(g,
                                [read1])
            nodes1.append(node)
        nodes2 = []
        for g in geneMers2:
            node = graph.add_node(g,
                                [read2])
            nodes2.append(node)
        for nodes in [nodes1, nodes2]:
            for n in range(len(nodes) - 1):
                mock_backward_edge = Edge(nodes[n],
                                        nodes[n+1],
                                        -1,
                                        -1)
                mock_rc_backward_edge = Edge(nodes[n+1],
                                            nodes[n],
                                            1,
                                            1)
                graph.add_edges_to_graph(mock_backward_edge,
                                        mock_rc_backward_edge)
                graph.add_edge_to_node(nodes[n],
                                    mock_backward_edge)
                graph.add_edge_to_node(nodes[n+1],
                                    mock_rc_backward_edge)
        # execution
        actual_backward_path_from_node = graph.get_backward_path_from_node(nodes1[0],
                                                                        True)
        actual_path_length = len(actual_backward_path_from_node)
        # assertion
        expected_backward_path_from_node = list(reversed([n.__hash__() for n in nodes1[:3]]))
        expected_path_length = 3
        self.assertEqual(actual_backward_path_from_node, expected_backward_path_from_node)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___get_linear_path_for_node(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        read1 = Read("read1",
                    genes1)
        read2 = Read("read2",
                    genes2)
        geneMers1 = [g for g in read1.get_geneMers(3)]
        geneMers2 = [g for g in read2.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodeHashes1 = []
        for g in range(len(geneMers1) - 1):
            sourceNode = graph.add_node(geneMers1[g],
                                    [read1])
            nodeHashes1.append(sourceNode.__hash__())
            sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers1[g],
                                                                        geneMers1[g+1])
        nodeHashes2 = []
        for g in range(len(geneMers2) - 1):
            sourceNode = graph.add_node(geneMers2[g],
                                    [read1])
            nodeHashes2.append(sourceNode.__hash__())
            sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers2[g],
                                                                        geneMers2[g+1])
        expected_paths = [nodeHashes1[0:2],
                        nodeHashes1[3:8],
                        nodeHashes2[3:8]]
        nodes_to_test = [nodeHashes1[1], nodeHashes1[3], nodeHashes2[3]]
        for n in range(len(nodes_to_test)):
            # execution
            actual_path = graph.get_linear_path_for_node(graph.get_node_by_hash(nodes_to_test[n]))
            # assertion
            self.assertTrue((actual_path == expected_paths[n] or actual_path == list(reversed(expected_paths[n]))))

    def test___get_linear_path_for_single_node(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "+gene5", "-gene6", "+gene7"]
        genes2 = ["+gene1", "-gene2", "+gene3", "+gene5", "-gene8", "+gene8"]
        read1 = Read("read1",
                    genes1)
        read2 = Read("read2",
                    genes2)
        geneMers1 = [g for g in read1.get_geneMers(3)]
        geneMers2 = [g for g in read2.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes1 = []
        for g in range(len(geneMers1) - 1):
            sourceNode = graph.add_node(geneMers1[g],
                                    [read1])
            nodes1.append(sourceNode)
            sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers1[g],
                                                                        geneMers1[g+1])
        for g in range(len(geneMers2) - 1):
            sourceNode = graph.add_node(geneMers2[g],
                                    [read2])
            sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers2[g],
                                                                        geneMers2[g+1])
        # execution
        actual_path = graph.get_linear_path_for_node(nodes1[0])
        actual_path_length = len(actual_path)
        # assertion
        expected_path = [nodes1[0].__hash__()]
        expected_path_length = 1
        self.assertEqual(actual_path, expected_path)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___insert_between_nodes_on_read(self):
        # setup
        # setup
        genes1 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3},
                            3)
        mock_readNodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 2, 1, 2]
        # execution
        actual_modifiedReadNodes = graph.insert_between_nodes_on_read(mock_readNodes,
                                                                        1,
                                                                        2,
                                                                        "*")
        # assertion
        expected_modifiedReadNodes = [1, "*", 2, 3, 4, 5, 6, 7, 8, 9, 10, 2, "*", 1, "*",2]
        self.assertEqual(expected_modifiedReadNodes, actual_modifiedReadNodes)

    def test___insert_after_node_at_end_of_read(self):
        # setup
        # setup
        genes1 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3},
                            3)
        mock_readNodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 2, 1, 1]
        # execution
        actual_modifiedReadNodes = graph.insert_between_nodes_on_read(mock_readNodes,
                                                                        1,
                                                                        2,
                                                                        "*")
        # assertion
        expected_modifiedReadNodes = [1, "*", 2, 3, 4, 5, 6, 7, 8, 9, 10, 2, "*", 1, 1, "*"]
        self.assertEqual(expected_modifiedReadNodes, actual_modifiedReadNodes)

    def test___insert_before_node_at_start_of_read(self):
        # setup
        # setup
        genes1 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3},
                            3)
        mock_readNodes = [2, 3, 4, 5, 6, 7, 8, 9, 10, 2, 1, 2]
        # execution
        actual_modifiedReadNodes = graph.insert_between_nodes_on_read(mock_readNodes,
                                                                        1,
                                                                        2,
                                                                        "*")
        # assertion
        expected_modifiedReadNodes = ["*", 2, 3, 4, 5, 6, 7, 8, 9, 10, 2, "*", 1, "*", 2]
        self.assertEqual(expected_modifiedReadNodes, actual_modifiedReadNodes)

    def test___insert_before_node_single_node(self):
        # setup
        genes1 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3},
                            3)
        mock_readNodes = [2]
        # execution
        actual_modifiedReadNodes = graph.insert_between_nodes_on_read(mock_readNodes,
                                                                        1,
                                                                        2,
                                                                        "*")
        # assertion
        expected_modifiedReadNodes = [2]
        self.assertEqual(expected_modifiedReadNodes, actual_modifiedReadNodes)

    def test___modify_readNode_list(self):
        # setup
        genes1 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3},
                            3)
        mock_readNodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 2, 1, 2]
        mock_replacementDict = {1: 11, "*": 12, 2: 13}
        # execution
        actual_modifiedReadNodes = graph.modify_readNode_list(mock_readNodes,
                                                            mock_replacementDict)
        # assertion
        expected_modifiedReadNodes = [1, "*", 2, 3, 4, 5, 6, 7, 8, 9, 10, 2, "*", 1, "*", 2]
        self.assertEqual(expected_modifiedReadNodes, actual_modifiedReadNodes)

    def test___replace_nodes_on_read(self):
        # setup
        genes1 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["-gene6", "+gene10", "+gene9", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3"]
        read1 = Read("read1",
                    genes1)
        read2 = Read("read2",
                    genes2)
        geneMers1 = [g for g in read1.get_geneMers(3)]
        geneMers2 = [g for g in read2.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes1 = []
        for g in range(len(geneMers1) - 1):
            sourceNode = graph.add_node(geneMers1[g],
                                    [read1])
            nodes1.append(sourceNode)
            graph.add_node_to_read(sourceNode,
                                "read1")
            sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers1[g],
                                                                        geneMers1[g+1])
        nodes2 = []
        for g in range(len(geneMers2) - 1):
            sourceNode = graph.add_node(geneMers2[g],
                                    [read2])
            graph.add_node_to_read(sourceNode,
                                "read2")
            nodes2.append(sourceNode)
            sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers2[g],
                                                                        geneMers2[g+1])
        replacementDict = {nodes2[0].__hash__(): nodes1[0].__hash__(),
                        nodes2[1].__hash__(): nodes1[1].__hash__(),
                        "*":  nodes1[2].__hash__(),
                        nodes2[2].__hash__():  nodes1[3].__hash__(),
                        nodes2[3].__hash__():  nodes1[4].__hash__()}
        # execution
        actual_modifiedReadNodes = graph.correct_read_nodes("read2",
                                                            replacementDict)
        actual_modified_length = len(actual_modifiedReadNodes)
        # assertion
        expected_modified_length = len(nodes2) + 1
        expected_modifiedReadNodes = graph.get_readNodes()["read1"][:6]
        self.assertEqual(actual_modified_length, expected_modified_length)
        self.assertEqual(actual_modifiedReadNodes, expected_modifiedReadNodes)

    def test__get_nodes_with_degree(self):
        # setup
        genes1 = ["-gene0", "+gene1", "-gene0", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene7", "-gene8", "-gene10", "+gene11", "-gene12", "+gene13"]
        genes2 = ["-gene", "+gene1","-gene0", "-gene2", "+gene3", "+gene5", "-gene6", "+gene7", "-gene8", "+gene9", "-gene10", "+gene11", "-gene12"]
        genes3 = ["-gene", "+gene1","-gene0", "-gene2", "+gene5", "-gene6", "+gene7", "-gene8", "+gene9", "-gene10", "+gene11", "-gene12"]
        read1 = Read("read1",
                    genes1)
        read2 = Read("read2",
                    genes2)
        read3 = Read("read3",
                    genes3)
        geneMers1 = [g for g in read1.get_geneMers(3)]
        geneMers2 = [g for g in read2.get_geneMers(3)]
        geneMers3 = [g for g in read3.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes1 = []
        for g in range(len(geneMers1)):
            sourceNode = graph.add_node(geneMers1[g],
                                    [read1])
            nodes1.append(sourceNode)
            graph.add_node_to_read(sourceNode,
                                "read1")
            if not g == len(geneMers1) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers1[g],
                                                                            geneMers1[g+1])
        nodes2 = []
        for g in range(len(geneMers2) - 1):
            sourceNode = graph.add_node(geneMers2[g],
                                    [read2])
            graph.add_node_to_read(sourceNode,
                                "read2")
            nodes2.append(sourceNode)
            if not g == len(geneMers2) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers2[g],
                                                                            geneMers2[g+1])
        nodes3 = []
        for g in range(len(geneMers3) - 1):
            sourceNode = graph.add_node(geneMers3[g],
                                    [read3])
            graph.add_node_to_read(sourceNode,
                                "read3")
            nodes3.append(sourceNode)
            if not g == len(geneMers3) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers3[g],
                                                                            geneMers3[g+1])
        # execution
        actual_two_degree_nodes = graph.get_nodes_with_degree(2)
        actual_three_degree_nodes = graph.get_nodes_with_degree(3)
        actual_four_degree_nodes =  graph.get_nodes_with_degree(4)
        # assertion
        expected_two_degree_nodes = nodes1[3:6] + nodes2[3:5] + nodes3[2:4] + nodes2[7:10] + nodes1[8:10]
        expected_three_degree_nodes = [nodes1[2], nodes1[7], nodes1[10]]
        expected_four_degree_nodes = [nodes1[1], nodes1[6]]
        self.assertEqual(len(actual_two_degree_nodes), 12)
        self.assertTrue(all(n in actual_two_degree_nodes for n in expected_two_degree_nodes))
        self.assertEqual(len(actual_three_degree_nodes), 3)
        self.assertTrue(all(n in actual_three_degree_nodes for n in expected_three_degree_nodes))
        self.assertEqual(len(actual_four_degree_nodes), 2)
        self.assertTrue(all(n in actual_four_degree_nodes for n in expected_four_degree_nodes))

    def test___make_replacement_dict_odd_one(self):
        # setup
        graph = GeneMerGraph({},
                            3)
        old_path = [1, 3]
        new_path = [1, 2, 3]
        # execution
        actual_replacementDict = graph.make_replacement_dict(old_path,
                                                            new_path)
        # assertion
        expected_replacementDict = {1: 1, "*": 2, 3: 3}
        self.assertEqual(actual_replacementDict, expected_replacementDict)

    def test___make_replacement_dict_odd_three(self):
        # setup
        graph = GeneMerGraph({},
                            3)
        old_path = [1, 2, 3, 4]
        new_path = [1, 5, 6, 7, 4]
        # execution
        actual_replacementDict = graph.make_replacement_dict(old_path,
                                                            new_path)
        # assertion
        expected_replacementDict = {1: 1, 2: 5, "*": 6, 3: 7, 4: 4}
        self.assertEqual(actual_replacementDict, expected_replacementDict)

    def test___make_replacement_dict_odd_three(self):
        # setup
        graph = GeneMerGraph({},
                            3)
        old_path = [1, 2, 3, 4, 5, 6]
        new_path = [1, 2, 7, 8, 9, 5, 6]
        # execution
        actual_replacementDict = graph.make_replacement_dict(old_path,
                                                            new_path)
        # assertion
        expected_replacementDict = {1: 1, 2: 2, 3: 7, "*": 8, 4: 9, 5: 5, 6: 6}
        self.assertEqual(actual_replacementDict, expected_replacementDict)

    def test___make_replacement_dict_even(self):
        # setup
        graph = GeneMerGraph({},
                            3)
        old_path = [1, 2, 3, 4]
        new_path = [1, 5, 6, 4]
        # assertion
        self.assertRaises(AssertionError, graph.make_replacement_dict, old_path, new_path)

    def test___make_replacement_dict_empty_old_path(self):
        # setup
        graph = GeneMerGraph({},
                            3)
        old_path = []
        new_path = ['x']
        # execution
        result = graph.make_replacement_dict(old_path, new_path)
        # assertion
        assert result == {'*': 'x'}

    def test___make_replacement_dict_empty_new_path(self):
        # setup
        graph = GeneMerGraph({},
                            3)
        old_path = ['a']
        new_path = []
        # assertion
        self.assertRaises(AssertionError, graph.make_replacement_dict, old_path, new_path)


    def test___make_replacement_dict_new_path_even(self):
        # setup
        graph = GeneMerGraph({},
                            3)
        old_path = ['a', 'b']
        new_path = ['x', 'y']
        # assertion
        self.assertRaises(AssertionError, graph.make_replacement_dict, old_path, new_path)

    def test___make_replacement_dict_old_path_longer_than_new_path(self):
        # setup
        graph = GeneMerGraph({},
                            3)
        old_path = ['a', 'b', 'c']
        new_path = ['x', 'y', 'z']
        self.assertRaises(AssertionError, graph.make_replacement_dict, old_path, new_path)

    def test___make_replacement_dict_old_path_has_star(self):
        # setup
        graph = GeneMerGraph({},
                            3)
        old_path = ['a', '*', 'b']
        new_path = ['x', 'y', 'z', 'w']
        # assertion
        self.assertRaises(AssertionError, graph.make_replacement_dict, old_path, new_path)

    def test___remove_short_linear_paths(self):
        # setup
        genes1 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes1)
        read2 = Read("read2",
                    genes2)
        read3 = Read("read3",
                    genes3)
        geneMers1 = [g for g in read1.get_geneMers(3)]
        geneMers2 = [g for g in read2.get_geneMers(3)]
        geneMers3 = [g for g in read3.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes1 = []
        for g in range(len(geneMers1)):
            sourceNode = graph.add_node(geneMers1[g],
                                    [read1])
            sourceNode.increment_node_coverage()
            nodes1.append(sourceNode)
            graph.add_node_to_read(sourceNode,
                                "read1")
            if not g == len(geneMers1) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers1[g],
                                                                            geneMers1[g+1])
        sourceNode = graph.add_node(geneMers1[-1],
                                    [read1])
        sourceNode.increment_node_coverage()
        nodes1.append(sourceNode)
        nodes2 = []
        for g in range(len(geneMers2) - 1):
            sourceNode = graph.add_node(geneMers2[g],
                                    [read2])
            sourceNode.increment_node_coverage()
            graph.add_node_to_read(sourceNode,
                                "read2")
            nodes2.append(sourceNode)
            if not g == len(geneMers2) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers2[g],
                                                                            geneMers2[g+1])
        sourceNode = graph.add_node(geneMers2[-1],
                                    [read2])
        sourceNode.increment_node_coverage()
        nodes2.append(sourceNode)
        nodes3 = []
        for g in range(len(geneMers3) - 1):
            sourceNode = graph.add_node(geneMers3[g],
                                    [read3])
            sourceNode.increment_node_coverage()
            graph.add_node_to_read(sourceNode,
                                "read3")
            nodes3.append(sourceNode)
            if not g == len(geneMers3) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers3[g],
                                                                            geneMers3[g+1])
        sourceNode = graph.add_node(geneMers3[-1],
                                    [read3])
        sourceNode.increment_node_coverage()
        nodes3.append(sourceNode)
        # execution
        actual_removed_nodeHashes = graph.remove_short_linear_paths(4)
        # assertion
        expected_number_removed_nodeHashes = 8
        expected_removed_nodeHashes = [nodes1[0].__hash__(),
                                    nodes1[1].__hash__(),
                                    nodes1[2].__hash__(),
                                    nodes2[0].__hash__(),
                                    nodes2[1].__hash__(),
                                    nodes2[2].__hash__(),
                                    nodes3[0].__hash__(),
                                    nodes3[1].__hash__()]
        self.assertEqual(len(actual_removed_nodeHashes), expected_number_removed_nodeHashes)
        self.assertTrue(all(h in actual_removed_nodeHashes for h in expected_removed_nodeHashes))
        for nodeHash in expected_removed_nodeHashes:
            self.assertTrue(nodeHash not in graph.get_nodes())

    def test___remove_short_linear_paths_longer_than_min(self):
        # setup
        genes1 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes1)
        read2 = Read("read2",
                    genes2)
        read3 = Read("read3",
                    genes3)
        geneMers1 = [g for g in read1.get_geneMers(3)]
        geneMers2 = [g for g in read2.get_geneMers(3)]
        geneMers3 = [g for g in read3.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes1 = []
        for g in range(len(geneMers1)):
            sourceNode = graph.add_node(geneMers1[g],
                                    [read1])
            sourceNode.increment_node_coverage()
            nodes1.append(sourceNode)
            graph.add_node_to_read(sourceNode,
                                "read1")
            if not g == len(geneMers1) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers1[g],
                                                                            geneMers1[g+1])
        sourceNode = graph.add_node(geneMers1[-1],
                                    [read1])
        sourceNode.increment_node_coverage()
        nodes1.append(sourceNode)
        nodes2 = []
        for g in range(len(geneMers2) - 1):
            sourceNode = graph.add_node(geneMers2[g],
                                    [read2])
            sourceNode.increment_node_coverage()
            graph.add_node_to_read(sourceNode,
                                "read2")
            nodes2.append(sourceNode)
            if not g == len(geneMers2) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers2[g],
                                                                            geneMers2[g+1])
        sourceNode = graph.add_node(geneMers2[-1],
                                    [read2])
        sourceNode.increment_node_coverage()
        nodes2.append(sourceNode)
        nodes3 = []
        for g in range(len(geneMers3) - 1):
            sourceNode = graph.add_node(geneMers3[g],
                                    [read3])
            sourceNode.increment_node_coverage()
            graph.add_node_to_read(sourceNode,
                                "read3")
            nodes3.append(sourceNode)
            if not g == len(geneMers3) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers3[g],
                                                                            geneMers3[g+1])
        sourceNode = graph.add_node(geneMers3[-1],
                                    [read3])
        sourceNode.increment_node_coverage()
        nodes3.append(sourceNode)
        # execution
        actual_removed_nodeHashes = graph.remove_short_linear_paths(3)
        # assertion
        expected_number_removed_nodeHashes = 2
        expected_removed_nodeHashes = [nodes3[0].__hash__(),
                                    nodes3[1].__hash__()]
        self.assertEqual(len(actual_removed_nodeHashes), expected_number_removed_nodeHashes)
        self.assertTrue(all(h in actual_removed_nodeHashes for h in expected_removed_nodeHashes))
        for nodeHash in expected_removed_nodeHashes:
            self.assertTrue(nodeHash not in graph.get_nodes())

    def test___remove_short_linear_paths_linear_path_of_length_one(self):
        # setup
        genes1 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5", "-gene12"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes1)
        read2 = Read("read2",
                    genes2)
        read3 = Read("read3",
                    genes3)
        geneMers1 = [g for g in read1.get_geneMers(3)]
        geneMers2 = [g for g in read2.get_geneMers(3)]
        geneMers3 = [g for g in read3.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3)
        nodes1 = []
        for g in range(len(geneMers1)):
            sourceNode = graph.add_node(geneMers1[g],
                                    [read1])
            sourceNode.increment_node_coverage()
            nodes1.append(sourceNode)
            graph.add_node_to_read(sourceNode,
                                "read1")
            if not g == len(geneMers1) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers1[g],
                                                                            geneMers1[g+1])
        sourceNode = graph.add_node(geneMers1[-1],
                                    [read1])
        sourceNode.increment_node_coverage()
        nodes1.append(sourceNode)
        nodes2 = []
        for g in range(len(geneMers2) - 1):
            sourceNode = graph.add_node(geneMers2[g],
                                    [read2])
            sourceNode.increment_node_coverage()
            graph.add_node_to_read(sourceNode,
                                "read2")
            nodes2.append(sourceNode)
            if not g == len(geneMers2) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers2[g],
                                                                            geneMers2[g+1])
        sourceNode = graph.add_node(geneMers2[-1],
                                    [read2])
        sourceNode.increment_node_coverage()
        nodes2.append(sourceNode)
        nodes3 = []
        for g in range(len(geneMers3) - 1):
            sourceNode = graph.add_node(geneMers3[g],
                                    [read3])
            sourceNode.increment_node_coverage()
            graph.add_node_to_read(sourceNode,
                                "read3")
            nodes3.append(sourceNode)
            if not g == len(geneMers3) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(geneMers3[g],
                                                                            geneMers3[g+1])
        sourceNode = graph.add_node(geneMers3[-1],
                                    [read3])
        sourceNode.increment_node_coverage()
        nodes3.append(sourceNode)
        # execution
        actual_removed_nodeHashes = graph.remove_short_linear_paths(2)
        # assertion
        expected_number_removed_nodeHashes = 1
        expected_removed_nodeHashes = [nodes2[4].__hash__()]
        self.assertEqual(len(actual_removed_nodeHashes), expected_number_removed_nodeHashes)
        self.assertTrue(all(h in actual_removed_nodeHashes for h in expected_removed_nodeHashes))
        for nodeHash in expected_removed_nodeHashes:
            self.assertTrue(nodeHash not in graph.get_nodes())

    def test___pop_bubbles(self):
        # setup
        genes1 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        graph = GeneMerGraph({"read1": genes1,
                            "read2": genes2},
                            5)
        # execution
        graph.generate_gml("test", 5, 1, 1)
        graph.pop_bubbles()
        # assertion
        self.assertEqual(graph.get_readNodes()["read1"], graph.get_readNodes()["read2"])