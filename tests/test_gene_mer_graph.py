import unittest
from collections import Counter

from amira.__main__ import parse_fastq
from amira.construct_edge import Edge
from amira.construct_graph import GeneMerGraph
from amira.construct_node import Node
from amira.construct_read import Read
from amira.path_finding_operations import construct_suffix_tree


class TestGeneMerGraphConstructor(unittest.TestCase):
    def test___init_empty_GeneMerGraph(self):
        # setup
        graph = GeneMerGraph({}, 0)
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
        graph = GeneMerGraph({"read1": genes, "read2": genes2}, 3)
        # execution
        actual_reads = graph.get_reads()
        actual_kmerSize = graph.get_kmerSize()
        actual_minGeneMerCoverage = graph.get_minNodeCoverage()
        actual_minEdgeCoverage = graph.get_minEdgeCoverage()
        actual_number_nodes = graph.get_total_number_of_nodes()
        actual_number_edges = graph.get_total_number_of_edges()
        actual_node_1_coverage = graph.get_nodes()[
            list(graph.get_nodes().keys())[0]
        ].get_node_coverage()
        actual_node_2_coverage = graph.get_nodes()[
            list(graph.get_nodes().keys())[1]
        ].get_node_coverage()
        actual_node_3_coverage = graph.get_nodes()[
            list(graph.get_nodes().keys())[2]
        ].get_node_coverage()
        # assertion
        expected_reads = {"read1": genes, "read2": genes2}
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
        genes = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene1",
            "-gene2",
            "+gene3",
            "+gene8",
        ]  # 5 nodes, 10 edges
        genes2 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene6",
            "+gene1",
            "-gene2",
            "+gene3",
        ]  # 3 nodes, 8 edges
        graph = GeneMerGraph({"read1": genes, "read2": genes2}, 3)
        # execution
        actual_reads = graph.get_reads()
        actual_kmerSize = graph.get_kmerSize()
        actual_minGeneMerCoverage = graph.get_minNodeCoverage()
        actual_minEdgeCoverage = graph.get_minEdgeCoverage()
        actual_number_nodes = graph.get_total_number_of_nodes()
        actual_number_edges = graph.get_total_number_of_edges()
        actual_node_1_coverage = graph.get_nodes()[
            list(graph.get_nodes().keys())[0]
        ].get_node_coverage()
        # assertion
        expected_reads = {"read1": genes, "read2": genes2}
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
        self.assertTrue(
            all(node.get_node_coverage() == 1 for node in list(graph.get_nodes().values())[1:])
        )
        self.assertTrue(
            all(
                edge.get_edge_coverage() == expected_edge_coverage
                for edge in list(graph.get_edges().values())
            )
        )

    def test__init_two_genemers_one_read(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]  # 2 nodes, 1 edge
        graph = GeneMerGraph({"read1": genes}, 3)
        # execution
        actual_reads = graph.get_reads()
        actual_kmerSize = graph.get_kmerSize()
        actual_minGeneMerCoverage = graph.get_minNodeCoverage()
        actual_minEdgeCoverage = graph.get_minEdgeCoverage()
        actual_number_nodes = graph.get_total_number_of_nodes()
        actual_number_edges = graph.get_total_number_of_edges()
        actual_node_1_coverage = graph.get_nodes()[
            list(graph.get_nodes().keys())[0]
        ].get_node_coverage()
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
        self.assertTrue(
            all(node.get_node_coverage() == 1 for node in list(graph.get_nodes().values())[1:])
        )
        self.assertTrue(
            all(
                edge.get_edge_coverage() == expected_edge_coverage
                for edge in list(graph.get_edges().values())
            )
        )

    def test___all_nodes(self):
        # setup
        genes = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "-gene3",
            "+gene2",
            "-gene1",
        ]
        read1 = Read("read1", genes)
        geneMers, geneMerPositions = read1.get_geneMers(3)
        graph = GeneMerGraph({}, 3)
        for g in geneMers:
            graph.add_node(g, [read1])
        # execution
        actual_all_nodes = [x for x in graph.all_nodes()]
        actual_number_of_nodes = len(actual_all_nodes)
        # assertion
        expected_number_of_nodes = 6
        self.assertEqual(actual_number_of_nodes, expected_number_of_nodes)

    def test___add_node_to_empty_graph(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1", genes)
        geneMers, geneMerPositions = read1.get_geneMers(3)
        geneMer = geneMers[0]
        graph = GeneMerGraph({}, 3)
        # execution
        actual_returned_node = graph.add_node(geneMer, [read1.get_readId()])
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
        read1 = Read("read1", genes1)
        geneMers1, geneMerPositions1 = read1.get_geneMers(3)
        graph = GeneMerGraph([], 3)
        node1 = graph.add_node(geneMers1[0], [read1.get_readId()])
        node1.increment_node_coverage()
        genes2 = ["+gene4", "-gene3", "+gene2"]
        read2 = Read("read2", genes2)
        geneMers2, geneMerPositions2 = read2.get_geneMers(3)
        geneMer2 = geneMers2[0]
        # execution
        actual_returned_node2 = graph.add_node(geneMer2, [read2.get_readId()])
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
        read1 = Read("read1", genes)
        geneMers, geneMerPositions = read1.get_geneMers(3)
        geneMer = geneMers[0]
        graph = GeneMerGraph([], 3)
        graph.add_node(geneMer, [read1.get_readId()])
        # execution
        actual_returned_node = graph.add_node(geneMer, [read1.get_readId()])
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
        read1 = Read("read1", genes1)
        geneMers1, geneMerPositions1 = read1.get_geneMers(3)
        geneMer1 = geneMers1[0]
        graph = GeneMerGraph([], 3)
        node1 = graph.add_node(geneMer1, [read1.get_readId()])
        node1.increment_node_coverage()
        genes2 = ["+gene4", "-gene3", "+gene2"]
        read2 = Read("read2", genes2)
        geneMers2, geneMerPositions2 = read2.get_geneMers(3)
        geneMer2 = geneMers2[0]
        node2 = graph.add_node(geneMer2, [read2.get_readId()])
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
        read1 = Read("read1", genes1)
        geneMer1 = [x for x in read1.get_geneMers(3)[0]][0]
        graph = GeneMerGraph([], 3)
        graph.add_node(geneMer1, [read1])
        genes2 = ["+gene4", "-gene3", "+gene2"]
        read2 = Read("read2", genes2)
        geneMer2 = [x for x in read2.get_geneMers(3)[0]][0]
        # assertion
        self.assertRaises(AssertionError, graph.get_node, geneMer2)

    def test___get_nodes_containing_subset(self):
        # setup
        genes = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "-gene3",
            "+gene2",
            "-gene1",
        ]
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        for g in geneMers:
            graph.add_node(g, [read1])
        selectedGenes = ["gene2", "gene6"]
        for g in range(len(selectedGenes)):
            # execution
            actual_selectedNodes = [n for n in graph.get_nodes_containing(selectedGenes[g])]
            actual_selectedGenes = [
                [g.get_name() for g in n.get_canonical_geneMer()] for n in actual_selectedNodes
            ]
            actual_geneMerCount = len(actual_selectedGenes)
            # assertion
            expected_geneMerCount = 3
            self.assertTrue(all(selectedGenes[g] in ng for ng in actual_selectedGenes))
            self.assertEqual(actual_geneMerCount, expected_geneMerCount)

    def test___get_nodes_containing_all(self):
        # setup
        genes = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "-gene3",
            "+gene2",
            "-gene1",
        ]
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        for g in geneMers:
            graph.add_node(g, [read1])
        selectedGenes = [
            "gene1",
            "gene2",
            "gene3",
            "gene4",
            "gene5",
            "gene6",
            "gene3",
            "gene2",
            "gene1",
        ]
        expected_counts = [1, 3, 5, 3, 3, 3, 5, 3, 1]
        for g in range(len(selectedGenes)):
            # execution
            actual_selectedNodes = [n for n in graph.get_nodes_containing(selectedGenes[g])]
            actual_selectedGenes = [
                [g.get_name() for g in n.get_canonical_geneMer()] for n in actual_selectedNodes
            ]
            actual_geneMerCount = len(actual_selectedGenes)
            # assertion
            expected_geneMerCount = expected_counts[g]
            self.assertTrue(all(selectedGenes[g] in ng for ng in actual_selectedGenes))
            self.assertEqual(actual_geneMerCount, expected_geneMerCount)

    def test___get_nodes_containing_gene_not_in_graph(self):
        # setup
        genes = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "-gene3",
            "+gene2",
            "-gene1",
        ]
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        for g in geneMers:
            graph.add_node(g, [read1])
        selectedGenes = "gene10"
        # execution
        actual_selectedNodes = [n for n in graph.get_nodes_containing(selectedGenes)]
        actual_selectedGenes = [
            [g.get_name() for g in n.get_canonical_geneMer()] for n in actual_selectedNodes
        ]
        # assertion
        expected_selectedGenes = []
        self.assertEqual(actual_selectedGenes, expected_selectedGenes)

    def test___get_nodes_containing_strand_list(self):
        # setup
        genes = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "-gene3",
            "+gene2",
            "-gene1",
        ]
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        for g in geneMers:
            graph.add_node(g, [read1])
        selectedGenes = ["gene2", "+gene6"]
        # assertion
        self.assertRaises(AssertionError, graph.get_nodes_containing, selectedGenes)

    def test___get_nodes_containing_strand_string(self):
        # setup
        genes = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "-gene3",
            "+gene2",
            "-gene1",
        ]
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        for g in geneMers:
            graph.add_node(g, [read1])
        selectedGenes = ["+gene6"]
        for g in range(len(selectedGenes)):
            # assertion
            self.assertRaises(AssertionError, graph.get_nodes_containing, selectedGenes[g])

    def test___create_edges_positive_to_positive(self):
        class fakeNode:
            def __init__(self, direction):
                self.direction = direction

            def get_direction(self):
                return self.direction

        # setup
        graph = GeneMerGraph({}, 3)
        mock_sourceNode = fakeNode(1)
        mock_targetNode = fakeNode(1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        # execution
        actual_sourceToTargetEdge, actual_reverseTargetToSourceEdge = graph.create_edges(
            mock_sourceNode, mock_targetNode, mock_sourceDirection, mock_targetDirection
        )
        actual_sourceToTargetEdge_sourceNode = actual_sourceToTargetEdge.get_sourceNode()
        actual_sourceToTargetEdge_targetNode = actual_sourceToTargetEdge.get_targetNode()
        actual_sourceToTargetEdge_sourceDirection = (
            actual_sourceToTargetEdge.get_sourceNodeDirection()
        )
        actual_sourceToTargetEdge_targetDirection = (
            actual_sourceToTargetEdge.get_targetNodeDirection()
        )
        actual_sourceToTargetEdge_coverage = actual_sourceToTargetEdge.get_edge_coverage()
        expected_sourceToTargetEdge_sourceNode = mock_sourceNode
        expected_sourceToTargetEdge_targetNode = mock_targetNode
        expected_sourceToTargetEdge_sourceDirection = mock_sourceDirection
        expected_sourceToTargetEdge_targetDirection = mock_targetDirection
        expected_sourceToTargetEdge_coverage = 0

        actual_reverseTargetToSourceEdge_sourceNode = (
            actual_reverseTargetToSourceEdge.get_sourceNode()
        )
        actual_reverseTargetToSourceEdge_targetNode = (
            actual_reverseTargetToSourceEdge.get_targetNode()
        )
        actual_reverseTargetToSourceEdge_sourceDirection = (
            actual_reverseTargetToSourceEdge.get_sourceNodeDirection()
        )
        actual_reverseTargetToSourceEdge_targetDirection = (
            actual_reverseTargetToSourceEdge.get_targetNodeDirection()
        )
        actual_reverseTargetToSourceEdge_coverage = (
            actual_reverseTargetToSourceEdge.get_edge_coverage()
        )
        expected_reverseTargetToSourceEdge_sourceNode = mock_targetNode
        expected_reverseTargetToSourceEdge_targetNode = mock_sourceNode
        expected_reverseTargetToSourceEdge_sourceDirection = mock_targetDirection * -1
        expected_reverseTargetToSourceEdge_targetDirection = mock_sourceDirection * -1
        expected_reverseTargetToSourceEdge_coverage = 0

        actual_sourceToTargetEdge_hash = actual_sourceToTargetEdge.__hash__()
        actual_reverseTargetToSourceEdge_hash = actual_reverseTargetToSourceEdge.__hash__()
        # assertion
        self.assertEqual(
            actual_sourceToTargetEdge_sourceNode, expected_sourceToTargetEdge_sourceNode
        )
        self.assertEqual(
            actual_sourceToTargetEdge_targetNode, expected_sourceToTargetEdge_targetNode
        )
        self.assertEqual(
            actual_sourceToTargetEdge_sourceDirection, expected_sourceToTargetEdge_sourceDirection
        )
        self.assertEqual(
            actual_sourceToTargetEdge_targetDirection, expected_sourceToTargetEdge_targetDirection
        )
        self.assertEqual(actual_sourceToTargetEdge_coverage, expected_sourceToTargetEdge_coverage)

        self.assertEqual(
            actual_reverseTargetToSourceEdge_sourceNode,
            expected_reverseTargetToSourceEdge_sourceNode,
        )
        self.assertEqual(
            actual_reverseTargetToSourceEdge_targetNode,
            expected_reverseTargetToSourceEdge_targetNode,
        )
        self.assertEqual(
            actual_reverseTargetToSourceEdge_sourceDirection,
            expected_reverseTargetToSourceEdge_sourceDirection,
        )
        self.assertEqual(
            actual_reverseTargetToSourceEdge_targetDirection,
            expected_reverseTargetToSourceEdge_targetDirection,
        )
        self.assertEqual(
            actual_reverseTargetToSourceEdge_coverage, expected_reverseTargetToSourceEdge_coverage
        )

        self.assertNotEqual(actual_sourceToTargetEdge_hash, actual_reverseTargetToSourceEdge_hash)

    def test___create_edges_negative_to_negative(self):
        class fakeNode:
            def __init__(self, direction):
                self.direction = direction

            def get_direction(self):
                return self.direction

        # setup
        graph = GeneMerGraph({}, 3)
        mock_sourceNode = fakeNode(-1)
        mock_targetNode = fakeNode(-1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        # execution
        actual_sourceToTargetEdge, actual_reverseTargetToSourceEdge = graph.create_edges(
            mock_sourceNode, mock_targetNode, mock_sourceDirection, mock_targetDirection
        )
        actual_sourceToTargetEdge_sourceNode = actual_sourceToTargetEdge.get_sourceNode()
        actual_sourceToTargetEdge_targetNode = actual_sourceToTargetEdge.get_targetNode()
        actual_sourceToTargetEdge_sourceDirection = (
            actual_sourceToTargetEdge.get_sourceNodeDirection()
        )
        actual_sourceToTargetEdge_targetDirection = (
            actual_sourceToTargetEdge.get_targetNodeDirection()
        )
        actual_sourceToTargetEdge_coverage = actual_sourceToTargetEdge.get_edge_coverage()
        expected_sourceToTargetEdge_sourceNode = mock_sourceNode
        expected_sourceToTargetEdge_targetNode = mock_targetNode
        expected_sourceToTargetEdge_sourceDirection = -1
        expected_sourceToTargetEdge_targetDirection = -1
        expected_sourceToTargetEdge_coverage = 0

        actual_reverseTargetToSourceEdge_sourceNode = (
            actual_reverseTargetToSourceEdge.get_sourceNode()
        )
        actual_reverseTargetToSourceEdge_targetNode = (
            actual_reverseTargetToSourceEdge.get_targetNode()
        )
        actual_reverseTargetToSourceEdge_sourceDirection = (
            actual_reverseTargetToSourceEdge.get_sourceNodeDirection()
        )
        actual_reverseTargetToSourceEdge_targetDirection = (
            actual_reverseTargetToSourceEdge.get_targetNodeDirection()
        )
        actual_reverseTargetToSourceEdge_coverage = (
            actual_reverseTargetToSourceEdge.get_edge_coverage()
        )
        expected_reverseTargetToSourceEdge_sourceNode = mock_targetNode
        expected_reverseTargetToSourceEdge_targetNode = mock_sourceNode
        expected_reverseTargetToSourceEdge_sourceDirection = 1
        expected_reverseTargetToSourceEdge_targetDirection = 1
        expected_reverseTargetToSourceEdge_coverage = 0

        actual_sourceToTargetEdge_hash = actual_sourceToTargetEdge.__hash__()
        actual_reverseTargetToSourceEdge_hash = actual_reverseTargetToSourceEdge.__hash__()
        # assertion
        self.assertEqual(
            actual_sourceToTargetEdge_sourceNode, expected_sourceToTargetEdge_sourceNode
        )
        self.assertEqual(
            actual_sourceToTargetEdge_targetNode, expected_sourceToTargetEdge_targetNode
        )
        self.assertEqual(
            actual_sourceToTargetEdge_sourceDirection, expected_sourceToTargetEdge_sourceDirection
        )
        self.assertEqual(
            actual_sourceToTargetEdge_targetDirection, expected_sourceToTargetEdge_targetDirection
        )
        self.assertEqual(actual_sourceToTargetEdge_coverage, expected_sourceToTargetEdge_coverage)

        self.assertEqual(
            actual_reverseTargetToSourceEdge_sourceNode,
            expected_reverseTargetToSourceEdge_sourceNode,
        )
        self.assertEqual(
            actual_reverseTargetToSourceEdge_targetNode,
            expected_reverseTargetToSourceEdge_targetNode,
        )
        self.assertEqual(
            actual_reverseTargetToSourceEdge_sourceDirection,
            expected_reverseTargetToSourceEdge_sourceDirection,
        )
        self.assertEqual(
            actual_reverseTargetToSourceEdge_targetDirection,
            expected_reverseTargetToSourceEdge_targetDirection,
        )
        self.assertEqual(
            actual_reverseTargetToSourceEdge_coverage, expected_reverseTargetToSourceEdge_coverage
        )

        self.assertNotEqual(actual_sourceToTargetEdge_hash, actual_reverseTargetToSourceEdge_hash)

    def test___create_edges_positive_to_negative(self):
        class fakeNode:
            def __init__(self, direction):
                self.direction = direction

            def get_direction(self):
                return self.direction

        # setup
        graph = GeneMerGraph({}, 3)
        mock_sourceNode = fakeNode(1)
        mock_targetNode = fakeNode(-1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        # execution
        actual_sourceToTargetEdge, actual_reverseTargetToSourceEdge = graph.create_edges(
            mock_sourceNode, mock_targetNode, mock_sourceDirection, mock_targetDirection
        )
        actual_sourceToTargetEdge_sourceNode = actual_sourceToTargetEdge.get_sourceNode()
        actual_sourceToTargetEdge_targetNode = actual_sourceToTargetEdge.get_targetNode()
        actual_sourceToTargetEdge_sourceDirection = (
            actual_sourceToTargetEdge.get_sourceNodeDirection()
        )
        actual_sourceToTargetEdge_targetDirection = (
            actual_sourceToTargetEdge.get_targetNodeDirection()
        )
        actual_sourceToTargetEdge_coverage = actual_sourceToTargetEdge.get_edge_coverage()
        expected_sourceToTargetEdge_sourceNode = mock_sourceNode
        expected_sourceToTargetEdge_targetNode = mock_targetNode
        expected_sourceToTargetEdge_sourceDirection = 1
        expected_sourceToTargetEdge_targetDirection = -1
        expected_sourceToTargetEdge_coverage = 0

        actual_reverseTargetToSourceEdge_sourceNode = (
            actual_reverseTargetToSourceEdge.get_sourceNode()
        )
        actual_reverseTargetToSourceEdge_targetNode = (
            actual_reverseTargetToSourceEdge.get_targetNode()
        )
        actual_reverseTargetToSourceEdge_sourceDirection = (
            actual_reverseTargetToSourceEdge.get_sourceNodeDirection()
        )
        actual_reverseTargetToSourceEdge_targetDirection = (
            actual_reverseTargetToSourceEdge.get_targetNodeDirection()
        )
        actual_reverseTargetToSourceEdge_coverage = (
            actual_reverseTargetToSourceEdge.get_edge_coverage()
        )
        expected_reverseTargetToSourceEdge_sourceNode = mock_targetNode
        expected_reverseTargetToSourceEdge_targetNode = mock_sourceNode
        expected_reverseTargetToSourceEdge_sourceDirection = mock_targetDirection * -1
        expected_reverseTargetToSourceEdge_targetDirection = mock_sourceDirection * -1
        expected_reverseTargetToSourceEdge_coverage = 0

        actual_sourceToTargetEdge_hash = actual_sourceToTargetEdge.__hash__()
        actual_reverseTargetToSourceEdge_hash = actual_reverseTargetToSourceEdge.__hash__()
        # assertion
        self.assertEqual(
            actual_sourceToTargetEdge_sourceNode, expected_sourceToTargetEdge_sourceNode
        )
        self.assertEqual(
            actual_sourceToTargetEdge_targetNode, expected_sourceToTargetEdge_targetNode
        )
        self.assertEqual(
            actual_sourceToTargetEdge_sourceDirection, expected_sourceToTargetEdge_sourceDirection
        )
        self.assertEqual(
            actual_sourceToTargetEdge_targetDirection, expected_sourceToTargetEdge_targetDirection
        )
        self.assertEqual(actual_sourceToTargetEdge_coverage, expected_sourceToTargetEdge_coverage)

        self.assertEqual(
            actual_reverseTargetToSourceEdge_sourceNode,
            expected_reverseTargetToSourceEdge_sourceNode,
        )
        self.assertEqual(
            actual_reverseTargetToSourceEdge_targetNode,
            expected_reverseTargetToSourceEdge_targetNode,
        )
        self.assertEqual(
            actual_reverseTargetToSourceEdge_sourceDirection,
            expected_reverseTargetToSourceEdge_sourceDirection,
        )
        self.assertEqual(
            actual_reverseTargetToSourceEdge_targetDirection,
            expected_reverseTargetToSourceEdge_targetDirection,
        )
        self.assertEqual(
            actual_reverseTargetToSourceEdge_coverage, expected_reverseTargetToSourceEdge_coverage
        )

        self.assertNotEqual(actual_sourceToTargetEdge_hash, actual_reverseTargetToSourceEdge_hash)

    def test___create_edges_negative_to_positive(self):
        class fakeNode:
            def __init__(self, direction):
                self.direction = direction

            def get_direction(self):
                return self.direction

        # setup
        graph = GeneMerGraph({}, 3)
        mock_sourceNode = fakeNode(1)
        mock_targetNode = fakeNode(-1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        # execution
        actual_sourceToTargetEdge, actual_reverseTargetToSourceEdge = graph.create_edges(
            mock_sourceNode, mock_targetNode, mock_sourceDirection, mock_targetDirection
        )
        actual_sourceToTargetEdge_sourceNode = actual_sourceToTargetEdge.get_sourceNode()
        actual_sourceToTargetEdge_targetNode = actual_sourceToTargetEdge.get_targetNode()
        actual_sourceToTargetEdge_sourceDirection = (
            actual_sourceToTargetEdge.get_sourceNodeDirection()
        )
        actual_sourceToTargetEdge_targetDirection = (
            actual_sourceToTargetEdge.get_targetNodeDirection()
        )
        actual_sourceToTargetEdge_coverage = actual_sourceToTargetEdge.get_edge_coverage()
        expected_sourceToTargetEdge_sourceNode = mock_sourceNode
        expected_sourceToTargetEdge_targetNode = mock_targetNode
        expected_sourceToTargetEdge_sourceDirection = 1
        expected_sourceToTargetEdge_targetDirection = -1
        expected_sourceToTargetEdge_coverage = 0

        actual_reverseTargetToSourceEdge_sourceNode = (
            actual_reverseTargetToSourceEdge.get_sourceNode()
        )
        actual_reverseTargetToSourceEdge_targetNode = (
            actual_reverseTargetToSourceEdge.get_targetNode()
        )
        actual_reverseTargetToSourceEdge_sourceDirection = (
            actual_reverseTargetToSourceEdge.get_sourceNodeDirection()
        )
        actual_reverseTargetToSourceEdge_targetDirection = (
            actual_reverseTargetToSourceEdge.get_targetNodeDirection()
        )
        actual_reverseTargetToSourceEdge_coverage = (
            actual_reverseTargetToSourceEdge.get_edge_coverage()
        )
        expected_reverseTargetToSourceEdge_sourceNode = mock_targetNode
        expected_reverseTargetToSourceEdge_targetNode = mock_sourceNode
        expected_reverseTargetToSourceEdge_sourceDirection = mock_targetDirection * -1
        expected_reverseTargetToSourceEdge_targetDirection = mock_sourceDirection * -1
        expected_reverseTargetToSourceEdge_coverage = 0

        actual_sourceToTargetEdge_hash = actual_sourceToTargetEdge.__hash__()
        actual_reverseTargetToSourceEdge_hash = actual_reverseTargetToSourceEdge.__hash__()
        # assertion
        self.assertEqual(
            actual_sourceToTargetEdge_sourceNode, expected_sourceToTargetEdge_sourceNode
        )
        self.assertEqual(
            actual_sourceToTargetEdge_targetNode, expected_sourceToTargetEdge_targetNode
        )
        self.assertEqual(
            actual_sourceToTargetEdge_sourceDirection, expected_sourceToTargetEdge_sourceDirection
        )
        self.assertEqual(
            actual_sourceToTargetEdge_targetDirection, expected_sourceToTargetEdge_targetDirection
        )
        self.assertEqual(actual_sourceToTargetEdge_coverage, expected_sourceToTargetEdge_coverage)

        self.assertEqual(
            actual_reverseTargetToSourceEdge_sourceNode,
            expected_reverseTargetToSourceEdge_sourceNode,
        )
        self.assertEqual(
            actual_reverseTargetToSourceEdge_targetNode,
            expected_reverseTargetToSourceEdge_targetNode,
        )
        self.assertEqual(
            actual_reverseTargetToSourceEdge_sourceDirection,
            expected_reverseTargetToSourceEdge_sourceDirection,
        )
        self.assertEqual(
            actual_reverseTargetToSourceEdge_targetDirection,
            expected_reverseTargetToSourceEdge_targetDirection,
        )
        self.assertEqual(
            actual_reverseTargetToSourceEdge_coverage, expected_reverseTargetToSourceEdge_coverage
        )

        self.assertNotEqual(actual_sourceToTargetEdge_hash, actual_reverseTargetToSourceEdge_hash)

    def test___add_missing_edge_to_edges(self):
        # setup

        class fakeNode:
            def __init__(self, direction):
                self.direction = direction

            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({}, 3)
        mock_sourceNode = fakeNode(-1)
        mock_targetNode = fakeNode(1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(
            mock_sourceNode, mock_targetNode, mock_sourceDirection, mock_targetDirection
        )
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
            def __init__(self, direction):
                self.direction = direction

            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({}, 3)
        mock_sourceNode = fakeNode(-1)
        mock_targetNode = fakeNode(1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(
            mock_sourceNode, mock_targetNode, mock_sourceDirection, mock_targetDirection
        )
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
            def __init__(self, direction):
                self.direction = direction

            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({}, 3)
        mock_sourceNode = fakeNode(-1)
        mock_targetNode = fakeNode(1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(
            mock_sourceNode, mock_targetNode, mock_sourceDirection, mock_targetDirection
        )
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
            def __init__(self, direction):
                self.direction = direction

            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({}, 3)
        mock_sourceNode = fakeNode(-1)
        mock_targetNode = fakeNode(1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(
            mock_sourceNode, mock_targetNode, mock_sourceDirection, mock_targetDirection
        )
        # execution
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(
            sourceToTargetEdge, reverseTargetToSourceEdge
        )
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
        self.assertTrue(
            all(
                actual_graph_edges[h].get_edge_coverage() == expected_edge_coverage
                for h in expected_edge_hashes
            )
        )

    def test___add_duplicate_edges_to_graph(self):
        # setup

        class fakeNode:
            def __init__(self, direction):
                self.direction = direction

            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({}, 3)
        mock_sourceNode = fakeNode(-1)
        mock_targetNode = fakeNode(1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(
            mock_sourceNode, mock_targetNode, mock_sourceDirection, mock_targetDirection
        )
        graph.add_edge_to_edges(sourceToTargetEdge)
        graph.add_edge_to_edges(reverseTargetToSourceEdge)
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        # execution
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(
            sourceToTargetEdge, reverseTargetToSourceEdge
        )
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
        self.assertTrue(
            all(
                actual_graph_edges[h].get_edge_coverage() == expected_edge_coverage
                for h in expected_edge_hashes
            )
        )

    def test___add_one_new_one_duplicate_edges_to_graph(self):
        # setup

        class fakeNode:
            def __init__(self, id, direction):
                self.id = id
                self.direction = direction

            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({}, 3)
        mock_sourceNode = fakeNode(0, -1)
        mock_targetNode = fakeNode(1, 1)
        mock_thirdNode = fakeNode(1, -1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        mock_thirdDirection = mock_thirdNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(
            mock_sourceNode, mock_targetNode, mock_sourceDirection, mock_targetDirection
        )
        targetToThirdEdge, reverseThirdToTargetEdge = graph.create_edges(
            mock_targetNode, mock_thirdNode, mock_targetDirection, mock_thirdDirection
        )
        graph.add_edge_to_edges(sourceToTargetEdge)
        graph.add_edge_to_edges(reverseTargetToSourceEdge)
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        # execution
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(
            sourceToTargetEdge, reverseTargetToSourceEdge
        )
        targetToThirdEdge, reverseThirdToTargetEdge = graph.add_edges_to_graph(
            targetToThirdEdge, reverseThirdToTargetEdge
        )
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        targetToThirdEdge.increment_edge_coverage()
        reverseThirdToTargetEdge.increment_edge_coverage()
        actual_graph_edges = graph.get_edges()
        actual_number_of_edges = len(actual_graph_edges)
        # assertion
        expected_number_of_edges = 4
        expected_edge_hashes = [
            sourceToTargetEdge.__hash__(),
            reverseTargetToSourceEdge.__hash__(),
            targetToThirdEdge.__hash__(),
            reverseThirdToTargetEdge.__hash__(),
        ]
        expected_source_target_edge_coverage = 2
        expected_target_third_edge_coverage = 1
        self.assertEqual(actual_number_of_edges, expected_number_of_edges)
        self.assertTrue(all(h in list(actual_graph_edges.keys()) for h in expected_edge_hashes))
        self.assertEqual(
            sourceToTargetEdge.get_edge_coverage(), expected_source_target_edge_coverage
        )
        self.assertEqual(
            reverseTargetToSourceEdge.get_edge_coverage(), expected_source_target_edge_coverage
        )
        self.assertEqual(targetToThirdEdge.get_edge_coverage(), expected_target_third_edge_coverage)
        self.assertEqual(
            reverseThirdToTargetEdge.get_edge_coverage(), expected_target_third_edge_coverage
        )

    def test___add_two_duplicate_reverse_edges_to_graph_all_positive(self):
        # setup

        class fakeNode:
            def __init__(self, id, direction):
                self.id = id
                self.direction = direction

            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({}, 3)
        mock_sourceNode = fakeNode(0, 1)
        mock_targetNode = fakeNode(1, 1)
        mock_thirdNode = fakeNode(1, 1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        mock_thirdDirection = mock_thirdNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(
            mock_sourceNode, mock_targetNode, mock_sourceDirection, mock_targetDirection
        )
        targetToThirdEdge, reverseThirdToTargetEdge = graph.create_edges(
            mock_targetNode, mock_thirdNode, mock_targetDirection, mock_thirdDirection
        )
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
        rc_sourceToTargetEdge.set_sourceNodeDirection(
            sourceToTargetEdge.get_sourceNodeDirection() * -1
        )
        rc_sourceToTargetEdge.set_targetNodeDirection(
            sourceToTargetEdge.get_targetNodeDirection() * -1
        )
        rc_reverseTargetToSourceEdge = reverseTargetToSourceEdge
        rc_reverseTargetToSourceEdge.set_sourceNodeDirection(
            reverseTargetToSourceEdge.get_sourceNodeDirection() * -1
        )
        rc_reverseTargetToSourceEdge.set_targetNodeDirection(
            reverseTargetToSourceEdge.get_targetNodeDirection() * -1
        )
        rc_targetToThirdEdge = targetToThirdEdge
        rc_targetToThirdEdge.set_sourceNodeDirection(
            targetToThirdEdge.get_sourceNodeDirection() * -1
        )
        rc_targetToThirdEdge.set_targetNodeDirection(
            targetToThirdEdge.get_targetNodeDirection() * -1
        )
        rc_reverseThirdToTargetEdge = reverseThirdToTargetEdge
        rc_reverseThirdToTargetEdge.set_sourceNodeDirection(
            reverseThirdToTargetEdge.get_sourceNodeDirection() * -1
        )
        reverseThirdToTargetEdge.set_targetNodeDirection(
            reverseThirdToTargetEdge.get_targetNodeDirection() * -1
        )
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(
            rc_sourceToTargetEdge, rc_reverseTargetToSourceEdge
        )
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        targetToThirdEdge, reverseThirdToTargetEdge = graph.add_edges_to_graph(
            rc_targetToThirdEdge, rc_reverseThirdToTargetEdge
        )
        targetToThirdEdge.increment_edge_coverage()
        reverseThirdToTargetEdge.increment_edge_coverage()
        actual_graph_edges = graph.get_edges()
        actual_number_of_edges = len(actual_graph_edges)
        # assertion
        expected_number_of_edges = 4
        expected_edge_hashes = [
            sourceToTargetEdge.__hash__(),
            reverseTargetToSourceEdge.__hash__(),
            targetToThirdEdge.__hash__(),
            reverseThirdToTargetEdge.__hash__(),
        ]
        expected_source_target_edge_coverage = 2
        expected_target_third_edge_coverage = 2
        self.assertEqual(actual_number_of_edges, expected_number_of_edges)
        self.assertTrue(all(h in list(actual_graph_edges.keys()) for h in expected_edge_hashes))
        self.assertEqual(
            sourceToTargetEdge.get_edge_coverage(), expected_source_target_edge_coverage
        )
        self.assertEqual(
            reverseTargetToSourceEdge.get_edge_coverage(), expected_source_target_edge_coverage
        )
        self.assertEqual(targetToThirdEdge.get_edge_coverage(), expected_target_third_edge_coverage)
        self.assertEqual(
            reverseThirdToTargetEdge.get_edge_coverage(), expected_target_third_edge_coverage
        )

    def test___add_two_duplicate_reverse_edges_to_graph_all_negative(self):
        # setup

        class fakeNode:
            def __init__(self, id, direction):
                self.id = id
                self.direction = direction

            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({}, 3)
        mock_sourceNode = fakeNode(0, -1)
        mock_targetNode = fakeNode(1, -1)
        mock_thirdNode = fakeNode(1, -1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        mock_thirdDirection = mock_thirdNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(
            mock_sourceNode, mock_targetNode, mock_sourceDirection, mock_targetDirection
        )
        targetToThirdEdge, reverseThirdToTargetEdge = graph.create_edges(
            mock_targetNode, mock_thirdNode, mock_targetDirection, mock_thirdDirection
        )
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
        rc_sourceToTargetEdge.set_sourceNodeDirection(
            sourceToTargetEdge.get_sourceNodeDirection() * -1
        )
        rc_sourceToTargetEdge.set_targetNodeDirection(
            sourceToTargetEdge.get_targetNodeDirection() * -1
        )
        rc_reverseTargetToSourceEdge = reverseTargetToSourceEdge
        rc_reverseTargetToSourceEdge.set_sourceNodeDirection(
            reverseTargetToSourceEdge.get_sourceNodeDirection() * -1
        )
        rc_reverseTargetToSourceEdge.set_targetNodeDirection(
            reverseTargetToSourceEdge.get_targetNodeDirection() * -1
        )
        rc_targetToThirdEdge = targetToThirdEdge
        rc_targetToThirdEdge.set_sourceNodeDirection(
            targetToThirdEdge.get_sourceNodeDirection() * -1
        )
        rc_targetToThirdEdge.set_targetNodeDirection(
            targetToThirdEdge.get_targetNodeDirection() * -1
        )
        rc_reverseThirdToTargetEdge = reverseThirdToTargetEdge
        rc_reverseThirdToTargetEdge.set_sourceNodeDirection(
            reverseThirdToTargetEdge.get_sourceNodeDirection() * -1
        )
        reverseThirdToTargetEdge.set_targetNodeDirection(
            reverseThirdToTargetEdge.get_targetNodeDirection() * -1
        )
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(
            rc_sourceToTargetEdge, rc_reverseTargetToSourceEdge
        )
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        targetToThirdEdge, reverseThirdToTargetEdge = graph.add_edges_to_graph(
            rc_targetToThirdEdge, rc_reverseThirdToTargetEdge
        )
        targetToThirdEdge.increment_edge_coverage()
        reverseThirdToTargetEdge.increment_edge_coverage()
        actual_graph_edges = graph.get_edges()
        actual_number_of_edges = len(actual_graph_edges)
        # assertion
        expected_number_of_edges = 4
        expected_edge_hashes = [
            sourceToTargetEdge.__hash__(),
            reverseTargetToSourceEdge.__hash__(),
            targetToThirdEdge.__hash__(),
            reverseThirdToTargetEdge.__hash__(),
        ]
        expected_source_target_edge_coverage = 2
        expected_target_third_edge_coverage = 2
        self.assertEqual(actual_number_of_edges, expected_number_of_edges)
        self.assertTrue(all(h in list(actual_graph_edges.keys()) for h in expected_edge_hashes))
        self.assertEqual(
            sourceToTargetEdge.get_edge_coverage(), expected_source_target_edge_coverage
        )
        self.assertEqual(
            reverseTargetToSourceEdge.get_edge_coverage(), expected_source_target_edge_coverage
        )
        self.assertEqual(targetToThirdEdge.get_edge_coverage(), expected_target_third_edge_coverage)
        self.assertEqual(
            reverseThirdToTargetEdge.get_edge_coverage(), expected_target_third_edge_coverage
        )

    def test___add_two_duplicate_reverse_edges_to_graph_one_positive_one_negative(self):
        # setup

        class fakeNode:
            def __init__(self, id, direction):
                self.id = id
                self.direction = direction

            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({}, 3)
        mock_sourceNode = fakeNode(0, 1)
        mock_targetNode = fakeNode(1, 1)
        mock_thirdNode = fakeNode(1, -1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        mock_thirdDirection = mock_thirdNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(
            mock_sourceNode, mock_targetNode, mock_sourceDirection, mock_targetDirection
        )
        targetToThirdEdge, reverseThirdToTargetEdge = graph.create_edges(
            mock_targetNode, mock_thirdNode, mock_targetDirection, mock_thirdDirection
        )
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
        rc_targetToThirdEdge.set_sourceNodeDirection(
            targetToThirdEdge.get_sourceNodeDirection() * -1
        )
        rc_targetToThirdEdge.set_targetNodeDirection(
            targetToThirdEdge.get_targetNodeDirection() * -1
        )
        rc_reverseThirdToTargetEdge = reverseThirdToTargetEdge
        rc_reverseThirdToTargetEdge.set_sourceNodeDirection(
            reverseThirdToTargetEdge.get_sourceNodeDirection() * -1
        )
        reverseThirdToTargetEdge.set_targetNodeDirection(
            reverseThirdToTargetEdge.get_targetNodeDirection() * -1
        )
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(
            sourceToTargetEdge, reverseTargetToSourceEdge
        )
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        targetToThirdEdge, reverseThirdToTargetEdge = graph.add_edges_to_graph(
            rc_targetToThirdEdge, rc_reverseThirdToTargetEdge
        )
        targetToThirdEdge.increment_edge_coverage()
        reverseThirdToTargetEdge.increment_edge_coverage()
        actual_graph_edges = graph.get_edges()
        actual_number_of_edges = len(actual_graph_edges)
        # assertion
        expected_number_of_edges = 4
        expected_edge_hashes = [
            sourceToTargetEdge.__hash__(),
            reverseTargetToSourceEdge.__hash__(),
            targetToThirdEdge.__hash__(),
            reverseThirdToTargetEdge.__hash__(),
        ]
        expected_source_target_edge_coverage = 2
        expected_target_third_edge_coverage = 2
        self.assertEqual(actual_number_of_edges, expected_number_of_edges)
        self.assertTrue(all(h in list(actual_graph_edges.keys()) for h in expected_edge_hashes))
        self.assertEqual(
            sourceToTargetEdge.get_edge_coverage(), expected_source_target_edge_coverage
        )
        self.assertEqual(
            reverseTargetToSourceEdge.get_edge_coverage(), expected_source_target_edge_coverage
        )
        self.assertEqual(targetToThirdEdge.get_edge_coverage(), expected_target_third_edge_coverage)
        self.assertEqual(
            reverseThirdToTargetEdge.get_edge_coverage(), expected_target_third_edge_coverage
        )

    def test___add_two_duplicate_reverse_edges_to_graph_one_negative_one_positive(self):
        # setup

        class fakeNode:
            def __init__(self, id, direction):
                self.id = id
                self.direction = direction

            def get_direction(self):
                return self.direction

        graph = GeneMerGraph({}, 3)
        mock_sourceNode = fakeNode(0, -1)
        mock_targetNode = fakeNode(1, 1)
        mock_thirdNode = fakeNode(1, 1)
        mock_sourceDirection = mock_sourceNode.get_direction()
        mock_targetDirection = mock_targetNode.get_direction()
        mock_thirdDirection = mock_thirdNode.get_direction()
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.create_edges(
            mock_sourceNode, mock_targetNode, mock_sourceDirection, mock_targetDirection
        )
        targetToThirdEdge, reverseThirdToTargetEdge = graph.create_edges(
            mock_targetNode, mock_thirdNode, mock_targetDirection, mock_thirdDirection
        )
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
        rc_sourceToTargetEdge.set_sourceNodeDirection(
            sourceToTargetEdge.get_sourceNodeDirection() * -1
        )
        rc_sourceToTargetEdge.set_targetNodeDirection(
            sourceToTargetEdge.get_targetNodeDirection() * -1
        )
        rc_reverseTargetToSourceEdge = reverseTargetToSourceEdge
        rc_reverseTargetToSourceEdge.set_sourceNodeDirection(
            reverseTargetToSourceEdge.get_sourceNodeDirection() * -1
        )
        rc_reverseTargetToSourceEdge.set_targetNodeDirection(
            reverseTargetToSourceEdge.get_targetNodeDirection() * -1
        )
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edges_to_graph(
            rc_sourceToTargetEdge, rc_reverseTargetToSourceEdge
        )
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        targetToThirdEdge, reverseThirdToTargetEdge = graph.add_edges_to_graph(
            targetToThirdEdge, reverseThirdToTargetEdge
        )
        targetToThirdEdge.increment_edge_coverage()
        reverseThirdToTargetEdge.increment_edge_coverage()
        actual_graph_edges = graph.get_edges()
        actual_number_of_edges = len(actual_graph_edges)
        # assertion
        expected_number_of_edges = 4
        expected_edge_hashes = [
            sourceToTargetEdge.__hash__(),
            reverseTargetToSourceEdge.__hash__(),
            targetToThirdEdge.__hash__(),
            reverseThirdToTargetEdge.__hash__(),
        ]
        expected_source_target_edge_coverage = 2
        expected_target_third_edge_coverage = 2
        self.assertEqual(actual_number_of_edges, expected_number_of_edges)
        self.assertTrue(all(h in list(actual_graph_edges.keys()) for h in expected_edge_hashes))
        self.assertEqual(
            sourceToTargetEdge.get_edge_coverage(), expected_source_target_edge_coverage
        )
        self.assertEqual(
            reverseTargetToSourceEdge.get_edge_coverage(), expected_source_target_edge_coverage
        )
        self.assertEqual(targetToThirdEdge.get_edge_coverage(), expected_target_third_edge_coverage)
        self.assertEqual(
            reverseThirdToTargetEdge.get_edge_coverage(), expected_target_third_edge_coverage
        )

    def test___add_edge_to_node_forward_empty(self):
        # setup
        class FakeEdge:
            def __init__(self, hash, sourceNodeDirection):
                self.hash = hash
                self.sourceNodeDirection = sourceNodeDirection

            def __hash__(self):
                return self.hash

            def get_sourceNodeDirection(self):
                return self.sourceNodeDirection

        mock_edge = FakeEdge(12345, +1)
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1", genes)
        node = [Node(x) for x in read1.get_geneMers(3)[0]][0]
        graph = GeneMerGraph({}, 3)
        # execution
        actual_updated_node = graph.add_edge_to_node(node, mock_edge)
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
            def __init__(self, hash, sourceNodeDirection):
                self.hash = hash
                self.sourceNodeDirection = sourceNodeDirection

            def __hash__(self):
                return self.hash

            def get_sourceNodeDirection(self):
                return self.sourceNodeDirection

        mock_edge = FakeEdge(12345, +1)
        mock_edge2 = FakeEdge(56789, +1)
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1", genes)
        node = [Node(x) for x in read1.get_geneMers(3)[0]][0]
        graph = GeneMerGraph({}, 3)
        graph.add_edge_to_node(node, mock_edge)
        # execution
        actual_updated_node = graph.add_edge_to_node(node, mock_edge2)
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
            def __init__(self, hash, sourceNodeDirection):
                self.hash = hash
                self.sourceNodeDirection = sourceNodeDirection

            def __hash__(self):
                return self.hash

            def get_sourceNodeDirection(self):
                return self.sourceNodeDirection

        mock_edge = FakeEdge(12345, -1)
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1", genes)
        node = [Node(x) for x in read1.get_geneMers(3)[0]][0]
        graph = GeneMerGraph({}, 3)
        # execution
        actual_updated_node = graph.add_edge_to_node(node, mock_edge)
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
            def __init__(self, hash, sourceNodeDirection):
                self.hash = hash
                self.sourceNodeDirection = sourceNodeDirection

            def __hash__(self):
                return self.hash

            def get_sourceNodeDirection(self):
                return self.sourceNodeDirection

        mock_edge = FakeEdge(12345, -1)
        mock_edge2 = FakeEdge(56789, -1)
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1", genes)
        node = [Node(x) for x in read1.get_geneMers(3)[0]][0]
        graph = GeneMerGraph({}, 3)
        graph.add_edge_to_node(node, mock_edge)
        # execution
        actual_updated_node = graph.add_edge_to_node(node, mock_edge2)
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
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        graph = GeneMerGraph({}, 3)
        graph.add_edge(sourceGeneMer, targetGeneMer)
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
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        thirdGeneMer = geneMers[2]
        graph = GeneMerGraph({}, 3)
        graph.add_edge(sourceGeneMer, targetGeneMer)
        graph.add_edge(targetGeneMer, thirdGeneMer)
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
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        thirdGeneMer = geneMers[2]
        fourthGeneMer = geneMers[3]
        graph = GeneMerGraph({}, 3)
        graph.add_edge(sourceGeneMer, targetGeneMer)
        graph.add_edge(targetGeneMer, thirdGeneMer)
        graph.add_edge(thirdGeneMer, fourthGeneMer)
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
        read1 = Read("read1", genes)
        read2 = Read("read1", genes2)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        geneMers2 = [x for x in read2.get_geneMers(3)[0]]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        targetGeneMer2 = geneMers2[1]
        graph = GeneMerGraph({}, 3)
        graph.add_edge(sourceGeneMer, targetGeneMer)
        graph.add_edge(sourceGeneMer, targetGeneMer2)
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
        read1 = Read("read1", genes)
        read2 = Read("read1", genes2)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        geneMers2 = [x for x in read2.get_geneMers(3)[0]]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        sourceGeneMer2 = geneMers2[0]
        graph = GeneMerGraph({}, 3)
        graph.add_edge(sourceGeneMer, targetGeneMer)
        graph.add_edge(sourceGeneMer2, targetGeneMer)
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
        read1 = Read("read1", genes)
        read2 = Read("read1", genes2)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        geneMers2 = [x for x in read2.get_geneMers(3)[0]]
        firstGeneMer = geneMers[0]
        middleGeneMer = geneMers[1]
        thirdGeneMer = geneMers[2]
        firstGeneMer2 = geneMers2[0]
        middleGeneMer2 = geneMers2[1]
        thirdGeneMer2 = geneMers2[2]
        graph = GeneMerGraph({}, 3)
        graph.add_edge(firstGeneMer, middleGeneMer)
        graph.add_edge(middleGeneMer, thirdGeneMer)
        graph.add_edge(firstGeneMer2, middleGeneMer2)
        graph.add_edge(middleGeneMer2, thirdGeneMer2)
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
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        graph = GeneMerGraph({}, 3)
        graph.add_edge(sourceGeneMer, targetGeneMer)
        graph.add_edge(sourceGeneMer, targetGeneMer)
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
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        graph = GeneMerGraph({}, 3)
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(sourceGeneMer, targetGeneMer)
        sourceNode = sourceToTargetEdge.get_sourceNode()
        targetNode = reverseTargetToSourceEdge.get_sourceNode()
        sourceNodeEdgeHashes = (
            sourceNode.get_forward_edge_hashes() + sourceNode.get_backward_edge_hashes()
        )
        targetNodeEdgeHashes = (
            targetNode.get_forward_edge_hashes() + targetNode.get_backward_edge_hashes()
        )
        # sanity checks
        self.assertNotEqual(sourceNodeEdgeHashes, [])
        self.assertNotEqual(targetNodeEdgeHashes, [])
        self.assertNotEqual(graph.get_edges(), {})
        # execution
        for edgeHash in sourceNodeEdgeHashes:
            graph.remove_edge(edgeHash)
        for edgeHash in targetNodeEdgeHashes:
            graph.remove_edge(edgeHash)
        actual_sourceNodeEdgeHashes = (
            sourceNode.get_forward_edge_hashes() + sourceNode.get_backward_edge_hashes()
        )
        actual_targetNodeEdgeHashes = (
            targetNode.get_forward_edge_hashes() + targetNode.get_backward_edge_hashes()
        )
        actual_edges = graph.get_edges()
        # assertion
        expected_sourceNodeEdgeHashes = []
        expected_targetNodeEdgeHashes = []
        expected_edges = {}
        self.assertEqual(actual_sourceNodeEdgeHashes, expected_sourceNodeEdgeHashes)
        self.assertEqual(actual_targetNodeEdgeHashes, expected_targetNodeEdgeHashes)
        self.assertEqual(actual_edges, expected_edges)

    def test___assign_Id_to_nodes(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        graph = GeneMerGraph({}, 3)
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(sourceGeneMer, targetGeneMer)
        # execution
        graph.assign_Id_to_nodes()
        # assertion
        expected_sourceNode_Id = 0
        expected_targetNode_Id = 1
        self.assertEqual(sourceToTargetEdge.get_sourceNode().get_node_Id(), expected_sourceNode_Id)
        self.assertEqual(
            reverseTargetToSourceEdge.get_sourceNode().get_node_Id(), expected_targetNode_Id
        )

    def test___write_gml_to_file(self):
        # setup
        import os

        graph = GeneMerGraph({}, 3)
        # execution
        graph.write_gml_to_file(
            "tests/test_graph", ["graph\t[", "\t\tnode\t[]", "\t\tedge\t[]", "]"]
        )
        # assertion
        self.assertTrue(os.path.exists("tests/test_graph.gml"))
        with open("tests/test_graph.gml", "r") as i:
            actual_content = i.read().splitlines()
        self.assertEqual(actual_content, ["graph\t[", "\t\tnode\t[]", "\t\tedge\t[]", "]"])
        os.remove("tests/test_graph.gml")

    def test___write_gml_to_file_in_subdirectory(self):
        # setup
        import os

        graph = GeneMerGraph({}, 3)
        # execution
        graph.write_gml_to_file(
            "tests/tests/test_graph", ["graph\t[", "\t\tnode\t[]", "\t\tedge\t[]", "]"]
        )
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
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        sourceGeneMer = geneMers[0]
        targetGeneMer = geneMers[1]
        graph = GeneMerGraph({}, 3)
        sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(sourceGeneMer, targetGeneMer)
        sourceToTargetEdge.increment_edge_coverage()
        reverseTargetToSourceEdge.increment_edge_coverage()
        sourceToTargetEdge.get_sourceNode().increment_node_coverage()
        reverseTargetToSourceEdge.get_sourceNode().increment_node_coverage()
        # execution
        actual_writtenGraph = graph.generate_gml("tests/test_graph", 3, 1, 1)
        # assertion
        node_prefix = "\tnode\t[\n\t\tid\t"
        edge_prefix = "\tedge\t[\n\t\tsource\t"
        gap = "\n\t\t"
        source_bw = "-gene3~~~+gene2~~~-gene1"
        source_fw = "+gene1~~~-gene2~~~+gene3"
        target_bw = "+gene4~~~-gene3~~~+gene2"
        target_fw = "-gene2~~~+gene3~~~-gene4"
        sd_str = "source_direction"
        td_str = "target_direction"
        expected_writtenGraph = [
            [
                "graph\t[",
                "multigraph 1",
                f'{node_prefix}0{gap}label\t"{source_bw}"{gap}coverage\t1{gap}reads\t""\n\t]',
                f"{edge_prefix}0{gap}target\t1{gap}{sd_str}\t-1{gap}{td_str}\t1{gap}weight\t1\n\t]",
                f'{node_prefix}1{gap}label\t"{target_bw}"{gap}coverage\t1{gap}reads\t""\n\t]',
                f"{edge_prefix}1{gap}target\t0{gap}{sd_str}\t1{gap}{td_str}\t-1{gap}weight\t1\n\t]",
                "]",
            ],
            [
                "graph\t[",
                "multigraph 1",
                f'{node_prefix}0{gap}label\t"{source_fw}"{gap}coverage\t1{gap}reads\t""\n\t]',
                f"{edge_prefix}0{gap}target\t1{gap}{sd_str}\t1{gap}{td_str}\t1{gap}weight\t1\n\t]",
                f'{node_prefix}1{gap}label\t"{target_fw}"{gap}coverage\t1{gap}reads\t""\n\t]',
                f"{edge_prefix}1{gap}target\t0{gap}{sd_str}\t-1{gap}{td_str}\t1{gap}weight\t1\n\t]",
                "]",
            ],
            [
                "graph\t[",
                "multigraph 1",
                f'{node_prefix}0{gap}label\t"{source_bw}"{gap}coverage\t1{gap}reads\t""\n\t]',
                f"{edge_prefix}0{gap}target\t1{gap}{sd_str}\t-1{gap}{td_str}\t1{gap}weight\t1\n\t]",
                f'{node_prefix}1{gap}label\t"{target_fw}"{gap}coverage\t1{gap}reads\t""\n\t]',
                f"{edge_prefix}1{gap}target\t0{gap}{sd_str}\t-1{gap}{td_str}\t1{gap}weight\t1\n\t]",
                "]",
            ],
            [
                "graph\t[",
                "multigraph 1",
                f'{node_prefix}0{gap}label\t"{source_fw}"{gap}coverage\t1{gap}reads\t""\n\t]',
                f"{edge_prefix}0{gap}target\t1{gap}{sd_str}\t1{gap}{td_str}\t-1{gap}weight\t1\n\t]",
                f'{node_prefix}1{gap}label\t"{target_bw}"{gap}coverage\t1{gap}reads\t""\n\t]',
                f"{edge_prefix}1{gap}target\t0{gap}{sd_str}\t1{gap}{td_str}\t-1{gap}weight\t1\n\t]",
                "]",
            ],
        ]
        import os

        self.assertTrue(os.path.exists("tests/test_graph.3.1.1.gml"))
        self.assertTrue(any(actual_writtenGraph == e for e in expected_writtenGraph))
        os.remove("tests/test_graph.3.1.1.gml")

    def test___get_gene_mer_label(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1", genes)
        sourceGeneMer = [x for x in read1.get_geneMers(3)[0]][0]
        sourceNode = Node(sourceGeneMer)
        graph = GeneMerGraph({}, 3)
        # execution
        actual_geneMerString = graph.get_gene_mer_label(sourceNode)
        # assertion
        expected_geneMerString = ["+gene1~~~-gene2~~~+gene3", "-gene3~~~+gene2~~~-gene1"]
        self.assertTrue(any(actual_geneMerString == e for e in expected_geneMerString))

    def test___filter_graph(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene10",
            "+gene9",
            "-gene6",
            "+gene3",
            "-gene7",
            "+gene5",
            "-gene6",
            "+gene3",
            "-gene7",
            "-gene6",
            "+gene3",
            "-gene7",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene3",
            "-gene4",
            "+gene5",
        ]
        genes2 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene9",
            "-gene6",
            "+gene7",
            "+gene3",
            "-gene4",
            "+gene5",
        ]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2}, 3)
        # execution
        graph.filter_graph(2, 2)
        graph.generate_gml("tests/test_graph", 3, 2, 2)
        actual_numberOfNodes = graph.get_total_number_of_nodes()
        actual_numberOfEdges = graph.get_total_number_of_edges()
        # assertion
        expected_numberOfNodes = 6
        expected_numberOfEdges = 10
        import os

        self.assertEqual(actual_numberOfNodes, expected_numberOfNodes)
        self.assertEqual(actual_numberOfEdges, expected_numberOfEdges)
        os.remove("tests/test_graph.3.2.2.gml")

    def test___filter_graph_k_1_cut_edge(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene10",
            "+gene9",
            "+gene3",
            "-gene7",
        ]
        genes2 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene10",
            "+gene9",
            "+gene3",
            "-gene7",
        ]
        graph = GeneMerGraph({"read1": genes1, "read2": genes1, "read3": genes2}, 1)
        # execution
        graph.filter_graph(1, 2)
        graph.generate_gml("tests/test_graph", 1, 1, 2)
        actual_numberOfNodes = graph.get_total_number_of_nodes()
        actual_numberOfEdges = graph.get_total_number_of_edges()
        # assertion
        expected_numberOfNodes = 9
        expected_numberOfEdges = 18
        import os

        self.assertEqual(actual_numberOfNodes, expected_numberOfNodes)
        self.assertEqual(actual_numberOfEdges, expected_numberOfEdges)
        os.remove("tests/test_graph.1.1.2.gml")

    def test___filter_all_graph(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene10",
            "+gene9",
            "-gene6",
            "+gene3",
            "-gene7",
            "+gene5",
            "-gene6",
            "+gene3",
            "-gene7",
            "-gene6",
            "+gene3",
            "-gene7",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene3",
            "-gene4",
            "+gene5",
        ]
        genes2 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene9",
            "-gene6",
            "+gene7",
            "+gene3",
            "-gene4",
            "+gene5",
        ]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2}, 3)
        # execution
        graph.filter_graph(10, 10)
        actual_writtenGraph = graph.generate_gml("tests/test_graph", 3, 10, 10)
        actual_numberOfNodes = graph.get_total_number_of_nodes()
        actual_numberOfEdges = graph.get_total_number_of_edges()
        # assertion
        expected_writtenGraph = ["graph\t[", "multigraph 1", "]"]
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
        read1 = Read("read1", genes)
        sourceGeneMer = [x for x in read1.get_geneMers(3)[0]][0]
        sourceNode = Node(sourceGeneMer)
        graph = GeneMerGraph({}, 3)
        # execution
        actual_listOfNodes = graph.add_node_to_read(
            sourceNode, "read1", sourceGeneMer.get_geneMerDirection()
        )
        actual_readNodeDict = graph.get_readNodes()
        # assertion
        expected_listOfNodes = [sourceNode.__hash__()]
        expected_readNodeDict = {"read1": [sourceNode.__hash__()]}
        self.assertEqual(actual_listOfNodes, expected_listOfNodes)
        self.assertEqual(actual_readNodeDict, expected_readNodeDict)

    def test___get_nodes_containing_read_filtered_graph(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene9",
            "-gene6",
            "+gene7",
            "+gene3",
            "-gene4",
            "+gene5",
        ]
        graph = GeneMerGraph({"read1": genes1}, 3)
        graph.filter_graph(2, 2)
        # execution
        actual_listOfNodes = graph.get_nodes_containing_read("read1")
        # assertion
        self.assertEqual(len(actual_listOfNodes), 2)

    def test___remove_node_from_reads_one_copy(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5"]
        read1 = Read("read1", genes)
        graph = GeneMerGraph({}, 3)
        sourceNodes = []
        for s in read1.get_geneMers(3)[0]:
            sourceNode = graph.add_node(s, [read1.get_readId()])
            graph.add_node_to_read(sourceNode, read1.get_readId(), s.get_geneMerDirection())
            sourceNodes.append(sourceNode)
        # execution
        test_node = sourceNodes[1]
        graph.remove_node_from_reads(test_node)
        # assertion
        expectedReadNodes = {"read1": [sourceNodes[0].__hash__(), None, sourceNodes[2].__hash__()]}
        self.assertEqual(graph.get_readNodes(), expectedReadNodes)

    def test___remove_node_from_reads_more_than_one_copy(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1", genes)
        graph = GeneMerGraph({}, 3)
        sourceNodes = []
        for s in read1.get_geneMers(3)[0]:
            sourceNode = graph.add_node(s, [read1.get_readId()])
            graph.add_node_to_read(sourceNode, read1.get_readId(), s.get_geneMerDirection())
            sourceNodes.append(sourceNode)
        # execution
        test_node = sourceNodes[1]
        graph.remove_node_from_reads(test_node)
        # assertion
        expectedReadNodes = {
            "read1": [
                sourceNodes[0].__hash__(),
                None,
                sourceNodes[2].__hash__(),
                sourceNodes[3].__hash__(),
                sourceNodes[4].__hash__(),
                None,
            ]
        }
        self.assertEqual(graph.get_readNodes(), expectedReadNodes)

    def test___get_existing_forward_node_from_node(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g, [read1.get_readId()])
            nodes.append(node)
        mock_forward_edge = Edge(nodes[0], nodes[1], 1, 1)
        mock_rc_forward_edge = Edge(nodes[1], nodes[0], -1, -1)
        graph.add_edges_to_graph(mock_forward_edge, mock_rc_forward_edge)
        graph.add_edge_to_node(nodes[0], mock_forward_edge)
        graph.add_edge_to_node(nodes[1], mock_rc_forward_edge)
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = graph.get_forward_node_from_node(
            nodes[0]
        )
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
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g, [read1.get_readId()])
            nodes.append(node)
        for n in range(len(nodes) - 1):
            mock_forward_edge = Edge(nodes[n], nodes[n + 1], 1, 1)
            mock_rc_forward_edge = Edge(nodes[n + 1], nodes[n], -1, -1)
            graph.add_edges_to_graph(mock_forward_edge, mock_rc_forward_edge)
            graph.add_edge_to_node(nodes[n], mock_forward_edge)
            graph.add_edge_to_node(nodes[n + 1], mock_rc_forward_edge)
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = graph.get_forward_node_from_node(
            nodes[0]
        )
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
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        geneMers1 = [x for x in read1.get_geneMers(3)[0]]
        geneMers2 = [x for x in read2.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes1 = []
        for g in geneMers1:
            node = graph.add_node(g, [read1.get_readId()])
            nodes1.append(node)
        nodes2 = []
        for g in geneMers2:
            node = graph.add_node(g, [read2.get_readId()])
            nodes2.append(node)
        for n in [nodes1, nodes2]:
            mock_forward_edge = Edge(n[1], n[2], 1, 1)
            mock_rc_forward_edge = Edge(n[2], n[1], -1, -1)
            graph.add_edges_to_graph(mock_forward_edge, mock_rc_forward_edge)
            graph.add_edge_to_node(n[1], mock_forward_edge)
            graph.add_edge_to_node(n[2], mock_rc_forward_edge)
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = graph.get_forward_node_from_node(
            nodes1[0]
        )
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
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g, [read1])
            nodes.append(node)
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = graph.get_forward_node_from_node(
            nodes[0]
        )
        # assertion
        expected_targetNode = None
        expected_targetDirection = None
        expected_extend = False
        self.assertEqual(actual_targetNode, expected_targetNode)
        self.assertEqual(actual_targetDirection, expected_targetDirection)
        self.assertEqual(actual_extend, expected_extend)

    def test___get_self_loop_forward_node_from_node(self):
        # setup
        genes = ["+gene1", "+gene1", "+gene1", "+gene1"]
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g, [read1])
            nodes.append(node)
        mock_forward_edge = Edge(nodes[0], nodes[1], 1, 1)
        mock_rc_forward_edge = Edge(nodes[1], nodes[0], -1, -1)
        graph.add_edges_to_graph(mock_forward_edge, mock_rc_forward_edge)
        graph.add_edge_to_node(nodes[0], mock_forward_edge)
        graph.add_edge_to_node(nodes[1], mock_rc_forward_edge)
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = graph.get_forward_node_from_node(
            nodes[0]
        )
        # assertion
        expected_targetNode = nodes[0]
        expected_targetDirection = 1
        expected_extend = False
        self.assertEqual(actual_targetNode, expected_targetNode)
        self.assertEqual(actual_targetDirection, expected_targetDirection)
        self.assertEqual(actual_extend, expected_extend)

    def test___get_existing_backward_node_from_node(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g, [read1])
            nodes.append(node)
        mock_backward_edge = Edge(nodes[0], nodes[1], -1, -1)
        mock_rc_backward_edge = Edge(nodes[1], nodes[0], 1, 1)
        graph.add_edges_to_graph(mock_backward_edge, mock_rc_backward_edge)
        graph.add_edge_to_node(nodes[0], mock_backward_edge)
        graph.add_edge_to_node(nodes[1], mock_rc_backward_edge)
        # execution
        (
            actual_extend,
            actual_targetNode,
            actual_targetDirection,
        ) = graph.get_backward_node_from_node(nodes[0])
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
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g, [read1])
            nodes.append(node)
        # execution
        (
            actual_extend,
            actual_targetNode,
            actual_targetDirection,
        ) = graph.get_backward_node_from_node(nodes[0])
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
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g, [read1])
            nodes.append(node)
        for n in range(len(nodes) - 1):
            mock_backward_edge = Edge(nodes[n], nodes[n + 1], -1, -1)
            mock_rc_backward_edge = Edge(nodes[n + 1], nodes[n], 1, 1)
            graph.add_edges_to_graph(mock_backward_edge, mock_rc_backward_edge)
            graph.add_edge_to_node(nodes[n], mock_backward_edge)
            graph.add_edge_to_node(nodes[n + 1], mock_rc_backward_edge)
        # execution
        (
            actual_extend,
            actual_targetNode,
            actual_targetDirection,
        ) = graph.get_backward_node_from_node(nodes[0])
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
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        geneMers1 = [x for x in read1.get_geneMers(3)[0]]
        geneMers2 = [x for x in read2.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes1 = []
        for g in geneMers1:
            node = graph.add_node(g, [read1])
            nodes1.append(node)
        nodes2 = []
        for g in geneMers2:
            node = graph.add_node(g, [read2])
            nodes2.append(node)
        for n in [nodes1, nodes2]:
            mock_backward_edge = Edge(n[1], n[2], -1, -1)
            mock_rc_backward_edge = Edge(n[2], n[1], 1, 1)
            graph.add_edges_to_graph(mock_backward_edge, mock_rc_backward_edge)
            graph.add_edge_to_node(n[1], mock_backward_edge)
            graph.add_edge_to_node(n[2], mock_rc_backward_edge)
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = graph.get_forward_node_from_node(
            nodes1[0]
        )
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
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g, [read1])
            nodes.append(node)
        for n in range(len(nodes) - 1):
            mock_forward_edge = Edge(nodes[n], nodes[n + 1], 1, 1)
            mock_rc_forward_edge = Edge(nodes[n + 1], nodes[n], -1, -1)
            graph.add_edges_to_graph(mock_forward_edge, mock_rc_forward_edge)
            graph.add_edge_to_node(nodes[n], mock_forward_edge)
            graph.add_edge_to_node(nodes[n + 1], mock_rc_forward_edge)
        # execution
        actual_forward_path_from_node = graph.get_forward_path_from_node(nodes[1], 1)
        actual_path_length = len(actual_forward_path_from_node)
        # assertion
        expected_forward_path_from_node = [n.__hash__() for n in nodes[1:]]
        expected_path_length = 5
        self.assertEqual(actual_forward_path_from_node, expected_forward_path_from_node)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___get_forward_path_from_node_circular(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene1", "-gene2", "+gene3"]
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g, [read1])
            nodes.append(node)
        for n in range(len(nodes) - 1):
            mock_forward_edge = Edge(nodes[n], nodes[n + 1], 1, 1)
            mock_rc_forward_edge = Edge(nodes[n + 1], nodes[n], -1, -1)
            graph.add_edges_to_graph(mock_forward_edge, mock_rc_forward_edge)
            graph.add_edge_to_node(nodes[n], mock_forward_edge)
            graph.add_edge_to_node(nodes[n + 1], mock_rc_forward_edge)
            mock_forward_edge.increment_edge_coverage()
            mock_rc_forward_edge.increment_edge_coverage()
        # execution
        actual_forward_path_from_node = graph.get_forward_path_from_node(nodes[0], 1, True)
        actual_path_length = len(actual_forward_path_from_node)
        # assertion
        expected_forward_path_from_node = [n.__hash__() for n in nodes]
        expected_path_length = 5
        self.assertEqual(actual_forward_path_from_node, expected_forward_path_from_node)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___get_forward_path_from_node_branched(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene7"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene8", "+gene8"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        geneMers1 = [x for x in read1.get_geneMers(3)[0]]
        geneMers2 = [x for x in read2.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes1 = []
        for g in geneMers1:
            node = graph.add_node(g, [read1])
            nodes1.append(node)
        nodes2 = []
        for g in geneMers2:
            node = graph.add_node(g, [read2])
            nodes2.append(node)
        for nodes in [nodes1, nodes2]:
            for n in range(len(nodes) - 1):
                mock_forward_edge = Edge(nodes[n], nodes[n + 1], 1, 1)
                mock_rc_forward_edge = Edge(nodes[n + 1], nodes[n], -1, -1)
                graph.add_edges_to_graph(mock_forward_edge, mock_rc_forward_edge)
                graph.add_edge_to_node(nodes[n], mock_forward_edge)
                graph.add_edge_to_node(nodes[n + 1], mock_rc_forward_edge)
        # execution
        actual_forward_path_from_node = graph.get_forward_path_from_node(nodes1[0], 1)
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
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        geneMers1 = [x for x in read1.get_geneMers(3)[0]]
        geneMers2 = [x for x in read2.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes1 = []
        for g in geneMers1:
            node = graph.add_node(g, [read1])
            nodes1.append(node)
        nodes2 = []
        for g in geneMers2:
            node = graph.add_node(g, [read2])
            nodes2.append(node)
        for nodes in [nodes1, nodes2]:
            for n in range(len(nodes) - 1):
                mock_forward_edge = Edge(nodes[n], nodes[n + 1], 1, 1)
                mock_rc_forward_edge = Edge(nodes[n + 1], nodes[n], -1, -1)
                graph.add_edges_to_graph(mock_forward_edge, mock_rc_forward_edge)
                graph.add_edge_to_node(nodes[n], mock_forward_edge)
                graph.add_edge_to_node(nodes[n + 1], mock_rc_forward_edge)
        # execution
        actual_forward_path_from_node = graph.get_forward_path_from_node(nodes1[0], 1, True)
        actual_path_length = len(actual_forward_path_from_node)
        # assertion
        expected_forward_path_from_node = [n.__hash__() for n in nodes1[:3]]
        expected_path_length = 3
        self.assertEqual(actual_forward_path_from_node, expected_forward_path_from_node)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___get_forward_path_from_middle_node_to_branched_want_branched(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene7"]
        genes2 = ["+gene0", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene8"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        geneMers1 = [x for x in read1.get_geneMers(3)[0]]
        geneMers2 = [x for x in read2.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes1 = []
        for g in geneMers1:
            node = graph.add_node(g, [read1])
            nodes1.append(node)
        nodes2 = []
        for g in geneMers2:
            node = graph.add_node(g, [read2])
            nodes2.append(node)
        for nodes in [nodes1, nodes2]:
            for n in range(len(nodes) - 1):
                mock_forward_edge = Edge(nodes[n], nodes[n + 1], 1, 1)
                mock_rc_forward_edge = Edge(nodes[n + 1], nodes[n], -1, -1)
                stedge, tsedge = graph.add_edges_to_graph(mock_forward_edge, mock_rc_forward_edge)
                graph.add_edge_to_node(nodes[n], mock_forward_edge)
                graph.add_edge_to_node(nodes[n + 1], mock_rc_forward_edge)
                stedge.increment_edge_coverage()
                tsedge.increment_edge_coverage()
        # execution
        actual_forward_path_from_node = graph.get_forward_path_from_node(nodes1[2], 1, True)
        actual_path_length = len(actual_forward_path_from_node)
        # assertion
        expected_forward_path_from_node = [n.__hash__() for n in nodes1[2:4]]
        expected_path_length = 2
        self.assertEqual(actual_forward_path_from_node, expected_forward_path_from_node)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___get_forward_path_from_branched_node_to_branched_want_branched(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
        genes2 = ["+gene0", "-gene2", "+gene3", "-gene4", "+gene5", "-gene7"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        geneMers1 = [x for x in read1.get_geneMers(3)[0]]
        geneMers2 = [x for x in read2.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes1 = []
        for g in geneMers1:
            node = graph.add_node(g, [read1])
            nodes1.append(node)
        nodes2 = []
        for g in geneMers2:
            node = graph.add_node(g, [read2])
            nodes2.append(node)
        for nodes in [nodes1, nodes2]:
            for n in range(len(nodes) - 1):
                mock_forward_edge = Edge(nodes[n], nodes[n + 1], 1, 1)
                mock_rc_forward_edge = Edge(nodes[n + 1], nodes[n], -1, -1)
                stedge, tsedge = graph.add_edges_to_graph(mock_forward_edge, mock_rc_forward_edge)
                graph.add_edge_to_node(nodes[n], mock_forward_edge)
                graph.add_edge_to_node(nodes[n + 1], mock_rc_forward_edge)
                stedge.increment_edge_coverage()
                tsedge.increment_edge_coverage()
        # execution
        actual_forward_path_from_node = graph.get_forward_path_from_node(nodes1[1], 1, True)
        actual_path_length = len(actual_forward_path_from_node)
        # assertion
        expected_forward_path_from_node = [n.__hash__() for n in nodes1[1:3]]
        expected_path_length = 2
        self.assertEqual(actual_forward_path_from_node, expected_forward_path_from_node)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___get_forward_path_from_terminal_node_to_branched_want_branched(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene7"]
        genes2 = ["+gene0", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene8"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        geneMers1 = [x for x in read1.get_geneMers(3)[0]]
        geneMers2 = [x for x in read2.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes1 = []
        for g in geneMers1:
            node = graph.add_node(g, [read1])
            nodes1.append(node)
        nodes2 = []
        for g in geneMers2:
            node = graph.add_node(g, [read2])
            nodes2.append(node)
        for nodes in [nodes1, nodes2]:
            for n in range(len(nodes) - 1):
                mock_forward_edge = Edge(nodes[n], nodes[n + 1], 1, 1)
                mock_rc_forward_edge = Edge(nodes[n + 1], nodes[n], -1, -1)
                stedge, tsedge = graph.add_edges_to_graph(mock_forward_edge, mock_rc_forward_edge)
                graph.add_edge_to_node(nodes[n], mock_forward_edge)
                graph.add_edge_to_node(nodes[n + 1], mock_rc_forward_edge)
                stedge.increment_edge_coverage()
                tsedge.increment_edge_coverage()
        # execution
        actual_forward_path_from_node = graph.get_forward_path_from_node(nodes2[0], 1, True)
        actual_path_length = len(actual_forward_path_from_node)
        # assertion
        expected_forward_path_from_node = [n.__hash__() for n in nodes2[0:2]]
        expected_path_length = 2
        self.assertEqual(actual_forward_path_from_node, expected_forward_path_from_node)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___get_backward_path_from_node_linear(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene7", "-gene8"]
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g, [read1])
            nodes.append(node)
        for n in range(len(nodes) - 1):
            mock_backward_edge = Edge(nodes[n], nodes[n + 1], -1, -1)
            mock_rc_backward_edge = Edge(nodes[n + 1], nodes[n], 1, 1)
            graph.add_edges_to_graph(mock_backward_edge, mock_rc_backward_edge)
            graph.add_edge_to_node(nodes[n], mock_backward_edge)
            graph.add_edge_to_node(nodes[n + 1], mock_rc_backward_edge)
        # execution
        actual_backward_path_from_node = graph.get_backward_path_from_node(nodes[1], -1)
        actual_path_length = len(actual_backward_path_from_node)
        # assertion
        expected_backward_path_from_node = list(reversed([n.__hash__() for n in nodes[1:]]))
        expected_path_length = 5
        self.assertEqual(actual_backward_path_from_node, expected_backward_path_from_node)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___get_backward_path_from_node_circular(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene1", "-gene2", "+gene3"]
        read1 = Read("read1", genes)
        geneMers = [x for x in read1.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g, [read1])
            nodes.append(node)
        for n in range(len(nodes) - 1):
            mock_backward_edge = Edge(nodes[n], nodes[n + 1], -1, -1)
            mock_rc_backward_edge = Edge(nodes[n + 1], nodes[n], 1, 1)
            graph.add_edges_to_graph(mock_backward_edge, mock_rc_backward_edge)
            graph.add_edge_to_node(nodes[n], mock_backward_edge)
            graph.add_edge_to_node(nodes[n + 1], mock_rc_backward_edge)
            mock_backward_edge.increment_edge_coverage()
            mock_rc_backward_edge.increment_edge_coverage()
        # execution
        actual_backward_path_from_node = graph.get_backward_path_from_node(nodes[0], -1, True)
        actual_path_length = len(actual_backward_path_from_node)
        # assertion
        expected_backward_path_from_node = list(reversed([n.__hash__() for n in nodes]))
        expected_path_length = 5
        self.assertEqual(actual_backward_path_from_node, expected_backward_path_from_node)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___get_backward_path_from_node_branched_want_branched(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene7"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene8", "+gene8"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        geneMers1 = [x for x in read1.get_geneMers(3)[0]]
        geneMers2 = [x for x in read2.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes1 = []
        for g in geneMers1:
            node = graph.add_node(g, [read1])
            nodes1.append(node)
        nodes2 = []
        for g in geneMers2:
            node = graph.add_node(g, [read2])
            nodes2.append(node)
        for nodes in [nodes1, nodes2]:
            for n in range(len(nodes) - 1):
                mock_backward_edge = Edge(nodes[n], nodes[n + 1], -1, -1)
                mock_rc_backward_edge = Edge(nodes[n + 1], nodes[n], 1, 1)
                graph.add_edges_to_graph(mock_backward_edge, mock_rc_backward_edge)
                graph.add_edge_to_node(nodes[n], mock_backward_edge)
                graph.add_edge_to_node(nodes[n + 1], mock_rc_backward_edge)
        # execution
        actual_backward_path_from_node = graph.get_backward_path_from_node(nodes1[0], -1, True)
        actual_path_length = len(actual_backward_path_from_node)
        # assertion
        expected_backward_path_from_node = list(reversed([n.__hash__() for n in nodes1[:3]]))
        expected_path_length = 3
        self.assertEqual(actual_backward_path_from_node, expected_backward_path_from_node)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___get_backward_path_from_middle_node_to_branched_want_branched(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene7"]
        genes2 = ["+gene0", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene8"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        geneMers1 = [x for x in read1.get_geneMers(3)[0]]
        geneMers2 = [x for x in read2.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes1 = []
        for g in geneMers1:
            node = graph.add_node(g, [read1])
            nodes1.append(node)
        nodes2 = []
        for g in geneMers2:
            node = graph.add_node(g, [read2])
            nodes2.append(node)
        for nodes in [nodes1, nodes2]:
            for n in range(len(nodes) - 1):
                mock_backward_edge = Edge(nodes[n], nodes[n + 1], -1, -1)
                mock_rc_backward_edge = Edge(nodes[n + 1], nodes[n], 1, 1)
                stedge, tsedge = graph.add_edges_to_graph(mock_backward_edge, mock_rc_backward_edge)
                graph.add_edge_to_node(nodes[n], mock_backward_edge)
                graph.add_edge_to_node(nodes[n + 1], mock_rc_backward_edge)
                stedge.increment_edge_coverage()
                tsedge.increment_edge_coverage()
        # execution
        actual_backward_path_from_node = graph.get_backward_path_from_node(nodes1[2], -1, True)
        actual_path_length = len(actual_backward_path_from_node)
        # assertion
        expected_backward_path_from_node = list(reversed([n.__hash__() for n in nodes1[2:4]]))
        expected_path_length = 2
        self.assertEqual(actual_backward_path_from_node, expected_backward_path_from_node)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___get_backward_path_from_branched_node_to_branched_want_branched(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
        genes2 = ["+gene0", "-gene2", "+gene3", "-gene4", "+gene5", "-gene7"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        geneMers1 = [x for x in read1.get_geneMers(3)[0]]
        geneMers2 = [x for x in read2.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes1 = []
        for g in geneMers1:
            node = graph.add_node(g, [read1])
            nodes1.append(node)
        nodes2 = []
        for g in geneMers2:
            node = graph.add_node(g, [read2])
            nodes2.append(node)
        for nodes in [nodes1, nodes2]:
            for n in range(len(nodes) - 1):
                mock_backward_edge = Edge(nodes[n], nodes[n + 1], -1, -1)
                mock_rc_backward_edge = Edge(nodes[n + 1], nodes[n], 1, 1)
                stedge, tsedge = graph.add_edges_to_graph(mock_backward_edge, mock_rc_backward_edge)
                graph.add_edge_to_node(nodes[n], mock_backward_edge)
                graph.add_edge_to_node(nodes[n + 1], mock_rc_backward_edge)
                stedge.increment_edge_coverage()
                tsedge.increment_edge_coverage()
        # execution
        actual_backward_path_from_node = graph.get_backward_path_from_node(nodes1[1], -1, True)
        actual_path_length = len(actual_backward_path_from_node)
        # assertion
        expected_backward_path_from_node = list(reversed([n.__hash__() for n in nodes1[1:3]]))
        expected_path_length = 2
        self.assertEqual(actual_backward_path_from_node, expected_backward_path_from_node)
        self.assertEqual(actual_path_length, expected_path_length)

    def test___get_linear_path_for_node(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene10",
            "+gene9",
            "-gene6",
            "+gene3",
            "-gene7",
            "+gene5",
            "-gene6",
            "+gene3",
            "-gene7",
            "-gene6",
            "+gene3",
            "-gene7",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene3",
            "-gene4",
            "+gene5",
        ]
        genes2 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene9",
            "-gene6",
            "+gene7",
            "+gene3",
            "-gene4",
            "+gene5",
        ]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        geneMers1 = [g for g in read1.get_geneMers(3)[0]]
        geneMers2 = [g for g in read2.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodeHashes1 = []
        for g in range(len(geneMers1) - 1):
            sourceNode = graph.add_node(geneMers1[g], [read1])
            nodeHashes1.append(sourceNode.__hash__())
            sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(
                geneMers1[g], geneMers1[g + 1]
            )
        nodeHashes2 = []
        for g in range(len(geneMers2) - 1):
            sourceNode = graph.add_node(geneMers2[g], [read1])
            nodeHashes2.append(sourceNode.__hash__())
            sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(
                geneMers2[g], geneMers2[g + 1]
            )
        expected_paths = [nodeHashes1[0:2], nodeHashes1[3:8], nodeHashes2[3:8]]
        nodes_to_test = [nodeHashes1[1], nodeHashes1[3], nodeHashes2[3]]
        for n in range(len(nodes_to_test)):
            # execution
            actual_path = graph.get_linear_path_for_node(graph.get_node_by_hash(nodes_to_test[n]))
            # assertion
            self.assertTrue(
                (
                    actual_path == expected_paths[n]
                    or actual_path == list(reversed(expected_paths[n]))
                )
            )

    def test___get_linear_path_for_single_node(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "+gene5", "-gene6", "+gene7"]
        genes2 = ["+gene1", "-gene2", "+gene3", "+gene5", "-gene8", "+gene8"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        geneMers1 = [g for g in read1.get_geneMers(3)[0]]
        geneMers2 = [g for g in read2.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes1 = []
        for g in range(len(geneMers1) - 1):
            sourceNode = graph.add_node(geneMers1[g], [read1])
            nodes1.append(sourceNode)
            sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(
                geneMers1[g], geneMers1[g + 1]
            )
        for g in range(len(geneMers2) - 1):
            sourceNode = graph.add_node(geneMers2[g], [read2])
            sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(
                geneMers2[g], geneMers2[g + 1]
            )
        # execution
        actual_path = graph.get_linear_path_for_node(nodes1[0])
        actual_path_length = len(actual_path)
        # assertion
        expected_path = [nodes1[0].__hash__()]
        expected_path_length = 1
        self.assertEqual(actual_path, expected_path)
        self.assertEqual(actual_path_length, expected_path_length)

    def test__get_nodes_with_degree(self):
        # setup
        genes1 = [
            "-gene0",
            "+gene1",
            "-gene0",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
            "-gene10",
            "+gene11",
            "-gene12",
            "+gene13",
        ]
        genes2 = [
            "-gene",
            "+gene1",
            "-gene0",
            "-gene2",
            "+gene3",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10",
            "+gene11",
            "-gene12",
        ]
        genes3 = [
            "-gene",
            "+gene1",
            "-gene0",
            "-gene2",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10",
            "+gene11",
            "-gene12",
        ]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        read3 = Read("read3", genes3)
        geneMers1 = [g for g in read1.get_geneMers(3)[0]]
        geneMers2 = [g for g in read2.get_geneMers(3)[0]]
        geneMers3 = [g for g in read3.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes1 = []
        for g in range(len(geneMers1)):
            sourceNode = graph.add_node(geneMers1[g], [read1])
            nodes1.append(sourceNode)
            graph.add_node_to_read(sourceNode, "read1", geneMers1[g].get_geneMerDirection())
            if not g == len(geneMers1) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(
                    geneMers1[g], geneMers1[g + 1]
                )
        nodes2 = []
        for g in range(len(geneMers2) - 1):
            sourceNode = graph.add_node(geneMers2[g], [read2])
            graph.add_node_to_read(sourceNode, "read2", geneMers1[g].get_geneMerDirection())
            nodes2.append(sourceNode)
            if not g == len(geneMers2) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(
                    geneMers2[g], geneMers2[g + 1]
                )
        nodes3 = []
        for g in range(len(geneMers3) - 1):
            sourceNode = graph.add_node(geneMers3[g], [read3])
            graph.add_node_to_read(sourceNode, "read3", geneMers3[g].get_geneMerDirection())
            nodes3.append(sourceNode)
            if not g == len(geneMers3) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(
                    geneMers3[g], geneMers3[g + 1]
                )
        # execution
        actual_two_degree_nodes = graph.get_nodes_with_degree(2)
        actual_three_degree_nodes = graph.get_nodes_with_degree(3)
        actual_four_degree_nodes = graph.get_nodes_with_degree(4)
        # assertion
        expected_two_degree_nodes = (
            nodes1[3:6] + nodes2[3:5] + nodes3[2:4] + nodes2[7:10] + nodes1[8:10]
        )
        expected_three_degree_nodes = [nodes1[2], nodes1[7], nodes1[10]]
        expected_four_degree_nodes = [nodes1[1], nodes1[6]]
        self.assertEqual(len(actual_two_degree_nodes), 12)
        self.assertTrue(all(n in actual_two_degree_nodes for n in expected_two_degree_nodes))
        self.assertEqual(len(actual_three_degree_nodes), 3)
        self.assertTrue(all(n in actual_three_degree_nodes for n in expected_three_degree_nodes))
        self.assertEqual(len(actual_four_degree_nodes), 2)
        self.assertTrue(all(n in actual_four_degree_nodes for n in expected_four_degree_nodes))

    def test___remove_short_linear_paths(self):
        # setup
        genes1 = [
            "-gene6",
            "+gene10",
            "+gene9",
            "-gene6",
            "+gene3",
            "-gene7",
            "+gene5",
            "-gene6",
            "+gene3",
            "-gene7",
            "-gene6",
            "+gene3",
            "-gene7",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene3",
            "-gene4",
            "+gene5",
        ]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        read3 = Read("read3", genes3)
        geneMers1 = [g for g in read1.get_geneMers(3)[0]]
        geneMers2 = [g for g in read2.get_geneMers(3)[0]]
        geneMers3 = [g for g in read3.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes1 = []
        for g in range(len(geneMers1)):
            sourceNode = graph.add_node(geneMers1[g], [read1.get_readId()])
            sourceNode.increment_node_coverage()
            nodes1.append(sourceNode)
            graph.add_node_to_read(sourceNode, "read1", geneMers1[g].get_geneMerDirection())
            if not g == len(geneMers1) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(
                    geneMers1[g], geneMers1[g + 1]
                )
        sourceNode = graph.add_node(geneMers1[-1], [read1.get_readId()])
        sourceNode.increment_node_coverage()
        nodes1.append(sourceNode)
        nodes2 = []
        for g in range(len(geneMers2) - 1):
            sourceNode = graph.add_node(geneMers2[g], [read2.get_readId()])
            sourceNode.increment_node_coverage()
            graph.add_node_to_read(sourceNode, "read2", geneMers2[g].get_geneMerDirection())
            nodes2.append(sourceNode)
            if not g == len(geneMers2) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(
                    geneMers2[g], geneMers2[g + 1]
                )
        sourceNode = graph.add_node(geneMers2[-1], [read2.get_readId()])
        sourceNode.increment_node_coverage()
        nodes2.append(sourceNode)
        nodes3 = []
        for g in range(len(geneMers3) - 1):
            sourceNode = graph.add_node(geneMers3[g], [read3.get_readId()])
            sourceNode.increment_node_coverage()
            graph.add_node_to_read(sourceNode, "read3", geneMers3[g].get_geneMerDirection())
            nodes3.append(sourceNode)
            if not g == len(geneMers3) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(
                    geneMers3[g], geneMers3[g + 1]
                )
        sourceNode = graph.add_node(geneMers3[-1], [read3.get_readId()])
        sourceNode.increment_node_coverage()
        nodes3.append(sourceNode)
        # execution
        actual_removed_nodeHashes = graph.remove_short_linear_paths(4)
        # assertion
        expected_number_removed_nodeHashes = 8
        expected_removed_nodeHashes = [
            nodes1[0].__hash__(),
            nodes1[1].__hash__(),
            nodes1[2].__hash__(),
            nodes2[0].__hash__(),
            nodes2[1].__hash__(),
            nodes2[2].__hash__(),
            nodes3[0].__hash__(),
            nodes3[1].__hash__(),
        ]
        self.assertEqual(len(actual_removed_nodeHashes), expected_number_removed_nodeHashes)
        self.assertTrue(all(h in actual_removed_nodeHashes for h in expected_removed_nodeHashes))
        for nodeHash in expected_removed_nodeHashes:
            self.assertTrue(nodeHash not in graph.get_nodes())

    def test___remove_short_linear_paths_longer_than_min(self):
        # setup
        genes1 = [
            "-gene6",
            "+gene10",
            "+gene9",
            "-gene6",
            "+gene3",
            "-gene7",
            "+gene5",
            "-gene6",
            "+gene3",
            "-gene7",
            "-gene6",
            "+gene3",
            "-gene7",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene3",
            "-gene4",
            "+gene5",
        ]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        read3 = Read("read3", genes3)
        geneMers1 = [g for g in read1.get_geneMers(3)[0]]
        geneMers2 = [g for g in read2.get_geneMers(3)[0]]
        geneMers3 = [g for g in read3.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes1 = []
        for g in range(len(geneMers1)):
            sourceNode = graph.add_node(geneMers1[g], [read1.get_readId()])
            sourceNode.increment_node_coverage()
            nodes1.append(sourceNode)
            graph.add_node_to_read(sourceNode, "read1", geneMers1[g].get_geneMerDirection())
            if not g == len(geneMers1) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(
                    geneMers1[g], geneMers1[g + 1]
                )
        sourceNode = graph.add_node(geneMers1[-1], [read1.get_readId()])
        sourceNode.increment_node_coverage()
        nodes1.append(sourceNode)
        nodes2 = []
        for g in range(len(geneMers2) - 1):
            sourceNode = graph.add_node(geneMers2[g], [read2.get_readId()])
            sourceNode.increment_node_coverage()
            graph.add_node_to_read(sourceNode, "read2", geneMers2[g].get_geneMerDirection())
            nodes2.append(sourceNode)
            if not g == len(geneMers2) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(
                    geneMers2[g], geneMers2[g + 1]
                )
        sourceNode = graph.add_node(geneMers2[-1], [read2.get_readId()])
        sourceNode.increment_node_coverage()
        nodes2.append(sourceNode)
        nodes3 = []
        for g in range(len(geneMers3) - 1):
            sourceNode = graph.add_node(geneMers3[g], [read3.get_readId()])
            sourceNode.increment_node_coverage()
            graph.add_node_to_read(sourceNode, "read3", geneMers3[g].get_geneMerDirection())
            nodes3.append(sourceNode)
            if not g == len(geneMers3) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(
                    geneMers3[g], geneMers3[g + 1]
                )
        sourceNode = graph.add_node(geneMers3[-1], [read3.get_readId()])
        sourceNode.increment_node_coverage()
        nodes3.append(sourceNode)
        # execution
        actual_removed_nodeHashes = graph.remove_short_linear_paths(3)
        # assertion
        expected_number_removed_nodeHashes = 2
        expected_removed_nodeHashes = [nodes3[0].__hash__(), nodes3[1].__hash__()]
        self.assertEqual(len(actual_removed_nodeHashes), expected_number_removed_nodeHashes)
        self.assertTrue(all(h in actual_removed_nodeHashes for h in expected_removed_nodeHashes))
        for nodeHash in expected_removed_nodeHashes:
            self.assertTrue(nodeHash not in graph.get_nodes())

    def test___remove_short_linear_paths_linear_path_of_length_one(self):
        # setup
        genes1 = [
            "-gene6",
            "+gene10",
            "+gene9",
            "-gene6",
            "+gene3",
            "-gene7",
            "+gene5",
            "-gene6",
            "+gene3",
            "-gene7",
            "-gene6",
            "+gene3",
            "-gene7",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene3",
            "-gene4",
            "+gene5",
        ]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5", "-gene12"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        read3 = Read("read3", genes3)
        geneMers1 = [g for g in read1.get_geneMers(3)[0]]
        geneMers2 = [g for g in read2.get_geneMers(3)[0]]
        geneMers3 = [g for g in read3.get_geneMers(3)[0]]
        graph = GeneMerGraph({}, 3)
        nodes1 = []
        for g in range(len(geneMers1)):
            sourceNode = graph.add_node(geneMers1[g], [read1.get_readId()])
            sourceNode.increment_node_coverage()
            nodes1.append(sourceNode)
            graph.add_node_to_read(sourceNode, "read1", geneMers1[g].get_geneMerDirection())
            if not g == len(geneMers1) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(
                    geneMers1[g], geneMers1[g + 1]
                )
        sourceNode = graph.add_node(geneMers1[-1], [read1.get_readId()])
        sourceNode.increment_node_coverage()
        nodes1.append(sourceNode)
        nodes2 = []
        for g in range(len(geneMers2) - 1):
            sourceNode = graph.add_node(geneMers2[g], [read2.get_readId()])
            sourceNode.increment_node_coverage()
            graph.add_node_to_read(sourceNode, "read2", geneMers2[g].get_geneMerDirection())
            nodes2.append(sourceNode)
            if not g == len(geneMers2) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(
                    geneMers2[g], geneMers2[g + 1]
                )
        sourceNode = graph.add_node(geneMers2[-1], [read2.get_readId()])
        sourceNode.increment_node_coverage()
        nodes2.append(sourceNode)
        nodes3 = []
        for g in range(len(geneMers3) - 1):
            sourceNode = graph.add_node(geneMers3[g], [read3.get_readId()])
            sourceNode.increment_node_coverage()
            graph.add_node_to_read(sourceNode, "read3", geneMers3[g].get_geneMerDirection())
            nodes3.append(sourceNode)
            if not g == len(geneMers3) - 1:
                sourceToTargetEdge, reverseTargetToSourceEdge = graph.add_edge(
                    geneMers3[g], geneMers3[g + 1]
                )
        sourceNode = graph.add_node(geneMers3[-1], [read3.get_readId()])
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

    def test___all_paths_for_subgraph_junctions(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "-gene6",
            "+gene7",
            "+gene9",
            "-gene10",
            "+gene16",
            "-gene17",
            "+gene18",
            "-gene19",
            "+gene20",
        ]
        genes2 = [
            "+gene11",
            "-gene12",
            "+gene3",
            "-gene4",
            "-gene6",
            "+gene13",
            "+gene14",
            "-gene15",
            "+gene16",
            "-gene17",
            "+gene18",
            "-gene21",
            "+gene22",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes2, "read3": genes1, "read4": genes2}, 3
        )
        for geneOfInterest in set(["gene4", "gene17"]):
            # get the graph nodes containing this gene
            nodeOfInterest = graph.get_nodes_containing(geneOfInterest)
            # get the node hashes containing this gene
            nodeHashesOfInterest = [n.__hash__() for n in nodeOfInterest]
            # get the nodes that contain this gene but are at the end of a path
            anchor_nodes, anchor_junctions = graph.get_anchors_of_interest(nodeHashesOfInterest)
            # execution
            actual_paths = graph.all_paths_for_subgraph(nodeHashesOfInterest, anchor_nodes)
            # assertion
            self.assertEqual(len(actual_paths), 6)
            self.assertTrue(all(len(actual_paths[p]) == 1 for p in actual_paths))
            self.assertTrue(all(len(actual_paths[p][0]) == 3 for p in actual_paths))

    def test___all_paths_for_subgraph_linear(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "-gene6",
            "+gene7",
            "+gene9",
            "-gene10",
            "+gene16",
            "-gene17",
            "+gene18",
            "-gene19",
            "+gene20",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes1, "read4": genes1}, 3
        )
        for geneOfInterest in set(["gene7"]):
            # get the graph nodes containing this gene
            nodeOfInterest = graph.get_nodes_containing(geneOfInterest)
            # get the node hashes containing this gene
            nodeHashesOfInterest = [n.__hash__() for n in nodeOfInterest]
            # get the nodes that contain this gene but are at the end of a path
            anchor_nodes, anchor_junctions = graph.get_anchors_of_interest(nodeHashesOfInterest)
            # execution
            actual_paths = graph.all_paths_for_subgraph(nodeHashesOfInterest, anchor_nodes)
            # assertion
            self.assertEqual(len(actual_paths), 1)
            self.assertTrue(all(len(actual_paths[p]) == 1 for p in actual_paths))
            self.assertTrue(all(len(actual_paths[p][0]) == 3 for p in actual_paths))

    def test___all_paths_for_subgraph_linear_duplicate(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "-gene6",
            "+gene7",
            "+gene9",
            "-gene7",
            "+gene16",
            "-gene17",
            "+gene18",
            "-gene19",
            "+gene20",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes1, "read4": genes1}, 3
        )
        for geneOfInterest in set(["gene7"]):
            # get the graph nodes containing this gene
            nodeOfInterest = graph.get_nodes_containing(geneOfInterest)
            # get the node hashes containing this gene
            nodeHashesOfInterest = [n.__hash__() for n in nodeOfInterest]
            # get the nodes that contain this gene but are at the end of a path
            anchor_nodes, anchor_junctions = graph.get_anchors_of_interest(nodeHashesOfInterest)
            # execution
            actual_paths = graph.all_paths_for_subgraph(nodeHashesOfInterest, anchor_nodes)
            # assertion
            self.assertEqual(len(actual_paths), 1)
            self.assertTrue(all(len(actual_paths[p]) == 1 for p in actual_paths))
            self.assertTrue(all(len(actual_paths[p][0]) == 5 for p in actual_paths))

    def test___find_read_boundaries(self):
        # setup
        genes1 = [
            "-gene0",
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10",
            "+gene11",
        ]
        genes2 = [
            "-gene0",
            "+gene1",
            "-gene2",
            "-gene4",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10",
            "+gene11",
        ]
        genes3 = [
            "-gene0",
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10",
            "+gene11",
        ]
        genes4 = [
            "-gene0",
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "-gene8",
            "+gene9",
            "-gene10",
            "+gene11",
        ]
        graph = GeneMerGraph(
            {
                "read1": genes1,
                "read2": genes2,
                "read3": genes3,
                "read4": genes1,
                "read5": genes1,
                "read6": genes1,
                "read7": genes1,
                "read8": genes1,
                "read9": genes1,
                "read10": genes1,
                "read11": genes1,
                "read12": genes3,
                "read13": genes3,
                "read14": genes3,
                "read15": genes3,
                "read16": genes3,
                "read17": genes3,
                "read18": genes3,
                "read19": genes3,
                "read20": genes3,
                "read21": genes3,
                "read22": genes4,
                "read23": genes1,
            },
            3,
        )
        # execution
        actual_start, actual_end = graph.find_read_boundaries(
            [None, 1, None, 2, 3, 4, 5, None, None, 6, 7, 8, 9, None, None, None]
        )
        # assertion
        self.assertEqual(actual_start, 1)
        self.assertEqual(actual_end, 12)

    def test_empty_insert_dict(self):
        graph = GeneMerGraph({}, 3)
        base_list = [(1, 1), (2, -1)]
        insert_dict = {}
        expected = [base_list]
        result = graph.insert_elements(base_list, insert_dict)
        self.assertEqual(result, expected)

    def test_single_insert(self):
        graph = GeneMerGraph({}, 3)
        base_list = [(1, 1), (2, -1)]
        insert_dict = {(0, 1): [[(1, 1), (3, 1), (2, -1)]]}
        expected = [[(1, 1), (3, 1), (2, -1)]]
        result = graph.insert_elements(base_list, insert_dict)
        self.assertEqual(result, expected)

    def test_multiple_inserts_single_path(self):
        graph = GeneMerGraph({}, 3)
        base_list = [(1, 1), (2, -1)]
        insert_dict = {(0, 1): [[(1, 1), (3, 1), (2, -1)], [(1, 1), (4, -1), (5, 1), (2, -1)]]}
        expected = [[(1, 1), (3, 1), (2, -1)], [(1, 1), (4, -1), (5, 1), (2, -1)]]
        result = graph.insert_elements(base_list, insert_dict)
        self.assertEqual(result, expected)

    def test_multiple_paths(self):
        graph = GeneMerGraph({}, 3)
        base_list = [(1, 1), (2, -1), (3, 1)]
        insert_dict = {
            (0, 1): [[(1, 1), (6, 1), (2, -1)], [(1, 1), (4, -1), (5, 1), (2, -1)]],
            (1, 2): [[(2, -1), (4, -1), (3, 1)], [(2, -1), (5, -1), (6, 1), (3, 1)]],
        }
        expected = sorted(
            [
                [(1, 1), (6, 1), (2, -1), (4, -1), (3, 1)],
                [(1, 1), (6, 1), (2, -1), (5, -1), (6, 1), (3, 1)],
                [(1, 1), (4, -1), (5, 1), (2, -1), (4, -1), (3, 1)],
                [(1, 1), (4, -1), (5, 1), (2, -1), (5, -1), (6, 1), (3, 1)],
            ]
        )
        result = graph.insert_elements(base_list, insert_dict)
        self.assertEqual(sorted(result), expected)

    def test___get_genes_in_unitig_length_one(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3"]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes1, "read4": genes1}, 3
        )
        node = [k for k, v in graph.get_nodes().items()]
        # execution
        actual_genes = graph.get_genes_in_unitig(node)
        # assertion
        self.assertTrue(
            actual_genes == ["+gene1", "-gene2", "+gene3"]
            or actual_genes == ["-gene3", "+gene2", "-gene1"]
        )

    def test___get_genes_in_unitig_length_greater_than_one(self):
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "-gene6",
            "+gene7",
            "+gene9",
            "-gene10",
            "+gene16",
            "-gene17",
            "+gene18",
            "-gene19",
            "+gene20",
        ]
        genes2 = [
            "+gene11",
            "-gene12",
            "+gene3",
            "-gene4",
            "-gene6",
            "+gene13",
            "+gene14",
            "-gene15",
            "+gene16",
            "-gene17",
            "+gene18",
            "-gene21",
            "+gene22",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes2, "read3": genes1, "read4": genes2}, 3
        )
        nodes = [n.__hash__() for n in graph.get_nodes_containing("gene15")]
        # execution
        actual_genes = graph.get_genes_in_unitig(nodes)
        # assertion
        self.assertTrue(
            actual_genes == ["+gene13", "+gene14", "-gene15", "+gene16", "-gene17"]
            or actual_genes == ["+gene17", "-gene16", "+gene15", "-gene14", "-gene13"]
        )

    def test___get_genes_in_unitig_length_zero(self):
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "-gene6",
            "+gene7",
            "+gene9",
            "-gene10",
            "+gene16",
            "-gene17",
            "+gene18",
            "-gene19",
            "+gene20",
        ]
        genes2 = [
            "+gene11",
            "-gene12",
            "+gene3",
            "-gene4",
            "-gene6",
            "+gene13",
            "+gene14",
            "-gene15",
            "+gene16",
            "-gene17",
            "+gene18",
            "-gene21",
            "+gene22",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes2, "read3": genes1, "read4": genes2}, 3
        )
        # execution
        actual_genes = graph.get_genes_in_unitig([])
        # assertion
        self.assertEqual(actual_genes, [])

    def test___reverse_list_of_genes(self):
        # setup
        graph = GeneMerGraph({}, 3)
        list_of_genes = ["-gene6", "+gene13", "+gene14", "-gene15"]
        # execution
        actual_reversed_genes = graph.reverse_list_of_genes(list_of_genes)
        # assertion
        self.assertEqual(actual_reversed_genes, ["+gene15", "-gene14", "-gene13", "+gene6"])

    def test___reverse_list_of_genes_empty(self):
        # setup
        graph = GeneMerGraph({}, 3)
        list_of_genes = []
        # execution
        actual_reversed_genes = graph.reverse_list_of_genes(list_of_genes)
        # assertion
        self.assertEqual(actual_reversed_genes, [])

    def test___needleman_wunsch_both_empty(self):
        # setup
        graph = GeneMerGraph({}, 3)
        # execution
        alignment = graph.needleman_wunsch([], [])
        # assertion
        self.assertEqual(alignment, [])

    def test___needleman_wunsch_gap_in_middle(self):
        # setup
        graph = GeneMerGraph({}, 3)
        first = ["+gene1", "-gene2", "+gene3"]
        second = ["+gene1", "+gene3"]
        # execution
        actual_alignment = graph.needleman_wunsch(first, second)
        # assertion
        expected_alignment = [("+gene1", "+gene1"), ("-gene2", "*"), ("+gene3", "+gene3")]
        self.assertEqual(actual_alignment, expected_alignment)

    def test___needleman_wunsch_two_gaps_in_middle(self):
        # setup
        graph = GeneMerGraph({}, 3)
        first = ["+gene1", "-gene2", "+gene3", "-gene4"]
        second = ["+gene1", "-gene4"]
        # execution
        actual_alignment = graph.needleman_wunsch(first, second)
        # assertion
        expected_alignment = [
            ("+gene1", "+gene1"),
            ("-gene2", "*"),
            ("+gene3", "*"),
            ("-gene4", "-gene4"),
        ]
        self.assertEqual(actual_alignment, expected_alignment)

    def test___needleman_wunsch_four_gaps_in_middle(self):
        # setup
        graph = GeneMerGraph({}, 3)
        first = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
        second = ["+gene1", "-gene6"]
        # execution
        actual_alignment = graph.needleman_wunsch(first, second)
        # assertion
        expected_alignment = [
            ("+gene1", "+gene1"),
            ("-gene2", "*"),
            ("+gene3", "*"),
            ("-gene4", "*"),
            ("+gene5", "*"),
            ("-gene6", "-gene6"),
        ]
        self.assertEqual(actual_alignment, expected_alignment)

    def test___needleman_wunsch_four_SNPs_in_middle(self):
        # setup
        graph = GeneMerGraph({}, 3)
        first = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
        second = ["+gene1", "-gene7", "+gene8", "-gene9", "+gene10", "-gene6"]
        # execution
        actual_alignment = graph.needleman_wunsch(first, second)
        # assertion
        expected_alignment = [
            ("+gene1", "+gene1"),
            ("-gene2", "-gene7"),
            ("+gene3", "+gene8"),
            ("-gene4", "-gene9"),
            ("+gene5", "+gene10"),
            ("-gene6", "-gene6"),
        ]
        self.assertEqual(actual_alignment, expected_alignment)

    def test___reverse_gene_alignment_no_gaps(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment = [
            ("+gene1", "+gene1"),
            ("-gene2", "-gene7"),
            ("+gene3", "+gene8"),
            ("-gene4", "-gene9"),
            ("+gene5", "+gene10"),
            ("-gene6", "-gene6"),
        ]
        # execution
        actual_reversed_alignment = graph.reverse_gene_alignment(alignment)
        # assertion
        expected_reversed_alignment = [
            ("+gene6", "+gene6"),
            ("-gene5", "-gene10"),
            ("+gene4", "+gene9"),
            ("-gene3", "-gene8"),
            ("+gene2", "+gene7"),
            ("-gene1", "-gene1"),
        ]
        self.assertEqual(actual_reversed_alignment, expected_reversed_alignment)

    def test___reverse_gene_alignment_gaps(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment = [
            ("+gene1", "+gene1"),
            ("-gene2", "*"),
            ("+gene3", "+gene8"),
            ("-gene4", "*"),
            ("*", "+gene10"),
            ("-gene6", "-gene6"),
        ]
        # execution
        actual_reversed_alignment = graph.reverse_gene_alignment(alignment)
        # assertion
        expected_reversed_alignment = [
            ("+gene6", "+gene6"),
            ("*", "-gene10"),
            ("+gene4", "*"),
            ("-gene3", "-gene8"),
            ("+gene2", "*"),
            ("-gene1", "-gene1"),
        ]
        self.assertEqual(actual_reversed_alignment, expected_reversed_alignment)

    def test___count_snps_in_alignment_snps(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment = [
            ("+gene1", "+gene1"),
            ("-gene2", "-gene7"),
            ("+gene3", "+gene8"),
            ("-gene4", "-gene9"),
            ("+gene5", "+gene10"),
            ("-gene6", "-gene6"),
        ]
        # execution
        actual_snps = graph.count_snps_in_alignment(alignment)
        # assertion
        self.assertEqual(actual_snps, 4)

    def test___count_snps_in_alignment_no_snps(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment = [
            ("+gene1", "+gene1"),
            ("-gene2", "-gene2"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
            ("+gene5", "+gene5"),
            ("-gene6", "-gene6"),
        ]
        # execution
        actual_snps = graph.count_snps_in_alignment(alignment)
        # assertion
        self.assertEqual(actual_snps, 0)

    def test___count_indels_in_alignment_no_indels(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment = [
            ("+gene1", "+gene1"),
            ("-gene2", "-gene2"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
            ("+gene5", "+gene5"),
            ("-gene6", "-gene6"),
        ]
        # execution
        actual_snps = graph.count_indels_in_alignment(alignment)
        # assertion
        self.assertEqual(actual_snps, 0)

    def test___count_indels_in_alignment_indels(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment = [
            ("+gene1", "+gene1"),
            ("*", "-gene2"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
            ("+gene5", "*"),
            ("-gene6", "-gene6"),
        ]
        # execution
        actual_snps = graph.count_indels_in_alignment(alignment)
        # assertion
        self.assertEqual(actual_snps, 2)

    def test___collect_reads_in_path(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "-gene6",
            "+gene7",
            "+gene9",
            "-gene10",
            "+gene16",
            "-gene17",
            "+gene18",
            "-gene19",
            "+gene20",
        ]
        genes2 = [
            "+gene11",
            "-gene12",
            "+gene3",
            "-gene4",
            "-gene6",
            "+gene13",
            "+gene14",
            "-gene15",
            "+gene16",
            "-gene17",
            "+gene18",
            "-gene21",
            "+gene22",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes2, "read3": genes1, "read4": genes2}, 3
        )
        nodes = [n.__hash__() for n in graph.get_nodes_containing("gene15")]
        # execution
        actual_reads = graph.collect_reads_in_path(nodes)
        # assertion
        self.assertEqual(actual_reads, set(["read2", "read4"]))

    def test___reorient_alignment_fw(self):
        # setup
        graph = GeneMerGraph({}, 3)
        fw_alignment = [
            ("+gene1", "+gene1"),
            ("*", "-gene2"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
            ("+gene5", "*"),
            ("-gene6", "-gene6"),
        ]
        rv_alignment = [
            ("+gene6", "+gene6"),
            ("-gene5", "*"),
            ("+gene4", "+gene4"),
            ("-gene3", "-gene3"),
            ("*", "+gene2"),
            ("-gene1", "-gene1"),
        ]
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "-gene6"]
        gene_mers = []
        reverse_gene_mers = []
        for i in range(len(genes) - (3 - 1)):
            # take a slice of the list of Genes from index i to i + kmerSize
            gene_mer = genes[i : i + 3]
            gene_mers.append(tuple(gene_mer))
            reverse_gene_mers.append(tuple(graph.reverse_list_of_genes(gene_mer)))
        fw_gene_mers_in_path = Counter(gene_mers)
        bw_gene_mers_in_path = Counter(reverse_gene_mers)
        # execution
        actual_alignment = graph.reorient_alignment(
            [tuple(["+gene3", "-gene4", "-gene6"])],
            fw_gene_mers_in_path,
            bw_gene_mers_in_path,
            fw_alignment,
            rv_alignment,
        )
        # assertion
        self.assertEqual(actual_alignment, fw_alignment)

    def test___reorient_alignment_rv(self):
        # setup
        graph = GeneMerGraph({}, 3)
        fw_alignment = [
            ("+gene1", "+gene1"),
            ("*", "-gene2"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
            ("+gene5", "*"),
            ("-gene6", "-gene6"),
        ]
        rv_alignment = [
            ("+gene6", "+gene6"),
            ("-gene5", "*"),
            ("+gene4", "+gene4"),
            ("-gene3", "-gene3"),
            ("*", "+gene2"),
            ("-gene1", "-gene1"),
        ]
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "-gene6"]
        gene_mers = []
        reverse_gene_mers = []
        for i in range(len(genes) - (3 - 1)):
            # take a slice of the list of Genes from index i to i + kmerSize
            gene_mer = genes[i : i + 3]
            gene_mers.append(tuple(gene_mer))
            reverse_gene_mers.append(tuple(graph.reverse_list_of_genes(gene_mer)))
        fw_gene_mers_in_path = Counter(gene_mers)
        bw_gene_mers_in_path = Counter(reverse_gene_mers)
        # execution
        actual_alignment = graph.reorient_alignment(
            [tuple(["+gene6", "+gene4", "-gene3"])],
            fw_gene_mers_in_path,
            bw_gene_mers_in_path,
            fw_alignment,
            rv_alignment,
        )
        # assertion
        self.assertEqual(actual_alignment, rv_alignment)

    def test___slice_alignment_by_shared_elements_unique_genes(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment = [
            ("+gene1", "+gene1"),
            ("*", "-gene2"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
            ("+gene5", "*"),
            ("-gene6", "-gene6"),
        ]
        read_genes = ["+gene3", "-gene4", "-gene7"]
        # execution
        actual_alignment_subset, actual_start_index, actual_end_index = (
            graph.slice_alignment_by_shared_elements(alignment, read_genes)
        )
        self.assertEqual(actual_alignment_subset, [("+gene3", "+gene3"), ("-gene4", "-gene4")])
        self.assertEqual(actual_start_index, 0)
        self.assertEqual(actual_end_index, 1)

    def test___slice_alignment_by_shared_elements_all_shared_genes(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment = [
            ("+gene1", "+gene1"),
            ("*", "-gene2"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
            ("+gene5", "*"),
            ("-gene6", "-gene6"),
        ]
        read_genes = ["+gene1", "-gene2", "+gene3", "-gene4", "-gene6"]
        # execution
        actual_alignment_subset, actual_start_index, actual_end_index = (
            graph.slice_alignment_by_shared_elements(alignment, read_genes)
        )
        self.assertEqual(actual_alignment_subset, alignment)
        self.assertEqual(actual_start_index, 0)
        self.assertEqual(actual_end_index, 4)

    def test___slice_alignment_by_shared_elements_one_shared_gene(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment = [
            ("+gene1", "+gene1"),
            ("*", "-gene2"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
            ("+gene5", "*"),
            ("-gene6", "-gene6"),
        ]
        read_genes = ["+gene7", "-gene8", "+gene3", "-gene9", "-gene10"]
        # execution
        actual_alignment_subset, actual_start_index, actual_end_index = (
            graph.slice_alignment_by_shared_elements(alignment, read_genes)
        )
        self.assertEqual(actual_alignment_subset, [("+gene3", "+gene3")])
        self.assertEqual(actual_start_index, 2)
        self.assertEqual(actual_end_index, 2)

    def test___slice_alignment_by_shared_elements_no_shared_genes(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment = [
            ("+gene1", "+gene1"),
            ("*", "-gene2"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
            ("+gene5", "*"),
            ("-gene6", "-gene6"),
        ]
        read_genes = ["+gene7", "-gene8", "-gene9"]
        # execution
        actual_alignment_subset, actual_start_index, actual_end_index = (
            graph.slice_alignment_by_shared_elements(alignment, read_genes)
        )
        self.assertEqual(actual_alignment_subset, [])
        self.assertEqual(actual_start_index, None)
        self.assertEqual(actual_end_index, None)

    def test___slice_alignment_by_shared_elements_duplicated_genes_in_alignment(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment = [
            ("+gene1", "+gene1"),
            ("*", "+gene1"),
            ("+gene1", "+gene1"),
            ("-gene4", "-gene4"),
            ("+gene5", "*"),
            ("-gene4", "-gene4"),
        ]
        read_genes = ["-gene0", "-gene4", "+gene1", "-gene4", "+gene5"]
        # execution
        actual_alignment_subset, actual_start_index, actual_end_index = (
            graph.slice_alignment_by_shared_elements(alignment, read_genes)
        )
        self.assertEqual(actual_alignment_subset, [("+gene1", "+gene1"), ("-gene4", "-gene4")])
        self.assertEqual(actual_start_index, 2)
        self.assertEqual(actual_end_index, 3)

    def test___slice_alignment_by_shared_elements_different_path(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment = [
            ("+gene1", "+gene1"),
            ("*", "-gene2"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
            ("+gene5", "*"),
            ("-gene6", "-gene6"),
        ]
        read_genes = ["-gene2", "-gene7", "-gene4", "-gene8", "+gene9"]
        # execution
        actual_alignment_subset, actual_start_index, actual_end_index = (
            graph.slice_alignment_by_shared_elements(alignment, read_genes)
        )
        self.assertEqual(actual_alignment_subset, [])
        self.assertEqual(actual_start_index, None)
        self.assertEqual(actual_end_index, None)

    def test___slice_alignment_by_shared_elements_multiple_chunks(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment = [
            ("+gene1", "+gene1"),
            ("*", "-gene2"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
            ("+gene5", "*"),
            ("-gene6", "-gene6"),
            ("+gene7", "+gene7"),
        ]
        read_genes = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
        ]
        # execution
        actual_alignment_subset, actual_start_index, actual_end_index = (
            graph.slice_alignment_by_shared_elements(alignment, read_genes)
        )
        self.assertEqual(actual_alignment_subset, alignment)
        self.assertEqual(actual_start_index, 0)
        self.assertEqual(actual_end_index, 6)

    def test___correct_genes_on_read_all_shared(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment_subset = [
            ("+gene1", "+gene1"),
            ("*", "-gene2"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
            ("+gene5", "*"),
            ("-gene6", "-gene6"),
        ]
        genes_on_read = ["+gene1", "-gene2", "+gene3", "-gene4", "-gene6"]
        # execution
        actual_corrected_reads = graph.correct_genes_on_read(
            genes_on_read, 0, 4, alignment_subset, "read1"
        )
        # assertion
        self.assertEqual(actual_corrected_reads, ["+gene1", "+gene3", "-gene4", "+gene5", "-gene6"])

    def test___correct_genes_on_read_subset_shared(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment_subset = [
            ("*", "-gene2"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
        ]
        genes_on_read = ["+gene7", "-gene2", "+gene3", "-gene4", "-gene8"]
        # execution
        actual_corrected_reads = graph.correct_genes_on_read(
            genes_on_read, 1, 3, alignment_subset, "read1"
        )
        # assertion
        self.assertEqual(actual_corrected_reads, ["+gene7", "+gene3", "-gene4", "-gene8"])

    def test___correct_genes_on_read_duplicates(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment_subset = [
            ("+gene1", "+gene1"),
            ("-gene4", "-gene4"),
        ]
        genes_on_read = ["-gene0", "-gene4", "+gene1", "-gene4", "+gene5"]
        # execution
        actual_corrected_reads = graph.correct_genes_on_read(
            genes_on_read, 2, 3, alignment_subset, "read1"
        )
        # assertion
        self.assertEqual(actual_corrected_reads, ["-gene0", "-gene4", "+gene1", "-gene4", "+gene5"])

    def test___get_gene_position_prefix_middle(self):
        # setup
        graph = GeneMerGraph({}, 3)
        gene_positions = [
            (0, 1000),
            (1001, 2000),
            (2001, 3000),
            (3001, 4000),
            (4001, 5000),
            (5001, 6000),
        ]
        first_shared_read_index = 2
        # execution
        actual_prefix = graph.get_gene_position_prefix(gene_positions, first_shared_read_index)
        # assertion
        self.assertEqual(actual_prefix, [(0, 1000), (1001, 2000)])

    def test___get_gene_position_prefix_start(self):
        # setup
        graph = GeneMerGraph({}, 3)
        gene_positions = [
            (0, 1000),
            (1001, 2000),
            (2001, 3000),
            (3001, 4000),
            (4001, 5000),
            (5001, 6000),
        ]
        first_shared_read_index = 0
        # execution
        actual_prefix = graph.get_gene_position_prefix(gene_positions, first_shared_read_index)
        # assertion
        self.assertEqual(actual_prefix, [])

    def test___get_gene_position_prefix_end(self):
        # setup
        graph = GeneMerGraph({}, 3)
        gene_positions = [
            (0, 1000),
            (1001, 2000),
            (2001, 3000),
            (3001, 4000),
            (4001, 5000),
            (5001, 6000),
        ]
        first_shared_read_index = 5
        # execution
        actual_prefix = graph.get_gene_position_prefix(gene_positions, first_shared_read_index)
        # assertion
        self.assertEqual(
            actual_prefix, [(0, 1000), (1001, 2000), (2001, 3000), (3001, 4000), (4001, 5000)]
        )

    def test___get_gene_position_suffix_middle(self):
        # setup
        graph = GeneMerGraph({}, 3)
        gene_positions = [
            (0, 1000),
            (1001, 2000),
            (2001, 3000),
            (3001, 4000),
            (4001, 5000),
            (5001, 6000),
        ]
        last_shared_read_index = 3
        # execution
        actual_suffix = graph.get_gene_position_suffix(gene_positions, last_shared_read_index)
        # assertion
        self.assertEqual(actual_suffix, [(4001, 5000), (5001, 6000)])

    def test___get_gene_position_suffix_end(self):
        # setup
        graph = GeneMerGraph({}, 3)
        gene_positions = [
            (0, 1000),
            (1001, 2000),
            (2001, 3000),
            (3001, 4000),
            (4001, 5000),
            (5001, 6000),
        ]
        last_shared_read_index = 5
        # execution
        actual_suffix = graph.get_gene_position_suffix(gene_positions, last_shared_read_index)
        # assertion
        self.assertEqual(actual_suffix, [])

    def test___get_gene_position_suffix_start(self):
        # setup
        graph = GeneMerGraph({}, 3)
        gene_positions = [
            (0, 1000),
            (1001, 2000),
            (2001, 3000),
            (3001, 4000),
            (4001, 5000),
            (5001, 6000),
        ]
        last_shared_read_index = 0
        # execution
        actual_suffix = graph.get_gene_position_suffix(gene_positions, last_shared_read_index)
        # assertion
        self.assertEqual(
            actual_suffix, [(1001, 2000), (2001, 3000), (3001, 4000), (4001, 5000), (5001, 6000)]
        )

    def test___get_gene_position_core_middle(self):
        # setup
        graph = GeneMerGraph({}, 3)
        gene_positions = [
            (0, 1000),
            (1001, 2000),
            (2001, 3000),
            (3001, 4000),
            (4001, 5000),
            (5001, 6000),
        ]
        first_shared_read_index = 2
        last_shared_read_index = 4
        # execution
        actual_core = graph.get_gene_position_core(
            gene_positions, first_shared_read_index, last_shared_read_index
        )
        # assertion
        self.assertEqual(actual_core, [(2001, 3000), (3001, 4000), (4001, 5000)])

    def test___get_gene_position_core_ends(self):
        # setup
        graph = GeneMerGraph({}, 3)
        gene_positions = [
            (0, 1000),
            (1001, 2000),
            (2001, 3000),
            (3001, 4000),
            (4001, 5000),
            (5001, 6000),
        ]
        first_shared_read_index = 0
        last_shared_read_index = 5
        # execution
        actual_core = graph.get_gene_position_core(
            gene_positions, first_shared_read_index, last_shared_read_index
        )
        # assertion
        self.assertEqual(actual_core, gene_positions)

    def test___get_new_gene_position_core_single_gap(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment_subset = [
            ("+gene1", "+gene1"),
            ("*", "-gene2"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
            ("+gene5", "*"),
            ("-gene6", "-gene6"),
        ]
        core_gene_positions = [(0, 1000), (1001, 2000), (2001, 3000), (3001, 4000), (5001, 6000)]
        # execution
        actual_new_core_gene_positions = graph.get_new_gene_position_core(
            alignment_subset, core_gene_positions
        )
        # assertion
        self.assertEqual(
            actual_new_core_gene_positions,
            [(0, 1000), (2001, 3000), (3001, 4000), (None, None), (5001, 6000)],
        )

    def test___get_new_gene_position_core_two_single_gaps(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment_subset = [
            ("+gene1", "+gene1"),
            ("-gene2", "*"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
            ("+gene5", "*"),
            ("-gene6", "-gene6"),
        ]
        core_gene_positions = [(0, 1000), (2001, 3000), (3001, 4000), (5001, 6000)]
        # execution
        actual_new_core_gene_positions = graph.get_new_gene_position_core(
            alignment_subset, core_gene_positions
        )
        # assertion
        self.assertEqual(
            actual_new_core_gene_positions,
            [(0, 1000), (None, None), (2001, 3000), (3001, 4000), (None, None), (5001, 6000)],
        )

    def test___get_new_gene_position_core_double_gap(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment_subset = [
            ("+gene1", "+gene1"),
            ("-gene2", "*"),
            ("+gene3", "*"),
            ("-gene4", "-gene4"),
            ("+gene5", "*"),
            ("-gene6", "-gene6"),
        ]
        core_gene_positions = [(0, 1000), (3001, 4000), (5001, 6000)]
        # execution
        actual_new_core_gene_positions = graph.get_new_gene_position_core(
            alignment_subset, core_gene_positions
        )
        # assertion
        self.assertEqual(
            actual_new_core_gene_positions,
            [(0, 1000), (None, None), (None, None), (3001, 4000), (None, None), (5001, 6000)],
        )

    def test___get_new_gene_position_core_gap_at_start(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment_subset = [
            ("+gene1", "*"),
            ("-gene2", "-gene2"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
            ("+gene5", "+gene5"),
            ("-gene6", "-gene6"),
        ]
        core_gene_positions = [(1001, 2000), (2001, 3000), (3001, 4000), (4001, 5000), (5001, 6000)]
        # execution
        actual_new_core_gene_positions = graph.get_new_gene_position_core(
            alignment_subset, core_gene_positions
        )
        # assertion
        self.assertEqual(
            actual_new_core_gene_positions,
            [(None, None), (1001, 2000), (2001, 3000), (3001, 4000), (4001, 5000), (5001, 6000)],
        )

    def test___get_new_gene_position_core_gap_at_end(self):
        # setup
        graph = GeneMerGraph({}, 3)
        alignment_subset = [
            ("+gene1", "+gene1"),
            ("-gene2", "-gene2"),
            ("+gene3", "+gene3"),
            ("-gene4", "-gene4"),
            ("+gene5", "+gene5"),
            ("-gene6", "*"),
        ]
        core_gene_positions = [
            (0, 1000),
            (1001, 2000),
            (2001, 3000),
            (3001, 4000),
            (4001, 5000),
            (5001, 6000),
        ]
        # execution
        actual_new_core_gene_positions = graph.get_new_gene_position_core(
            alignment_subset, core_gene_positions
        )
        # assertion
        self.assertEqual(
            actual_new_core_gene_positions,
            [(0, 1000), (1001, 2000), (2001, 3000), (3001, 4000), (4001, 5000), (None, None)],
        )

    def test___join_gene_position_ends_with_core_all(self):
        # setup
        graph = GeneMerGraph({}, 3)
        position_prefix = [1, 2]
        position_suffix = [5, 6]
        core_gene_positions = [3, 4]
        # execution
        actual_new_positions = graph.join_gene_position_ends_with_core(
            position_prefix, position_suffix, core_gene_positions
        )
        # assertion
        self.assertEqual(actual_new_positions, [1, 2, 3, 4, 5, 6])

    def test___join_gene_position_ends_with_core_no_prefix(self):
        # setup
        graph = GeneMerGraph({}, 3)
        position_prefix = []
        position_suffix = [5, 6]
        core_gene_positions = [3, 4]
        # execution
        actual_new_positions = graph.join_gene_position_ends_with_core(
            position_prefix, position_suffix, core_gene_positions
        )
        # assertion
        self.assertEqual(actual_new_positions, [3, 4, 5, 6])

    def test___join_gene_position_ends_with_core_no_suffix(self):
        # setup
        graph = GeneMerGraph({}, 3)
        position_prefix = [1, 2]
        position_suffix = []
        core_gene_positions = [3, 4]
        # execution
        actual_new_positions = graph.join_gene_position_ends_with_core(
            position_prefix, position_suffix, core_gene_positions
        )
        # assertion
        self.assertEqual(actual_new_positions, [1, 2, 3, 4])

    def test___join_gene_position_ends_with_core_no_core(self):
        # setup
        graph = GeneMerGraph({}, 3)
        position_prefix = [1, 2]
        position_suffix = [5, 6]
        core_gene_positions = []
        # execution
        actual_new_positions = graph.join_gene_position_ends_with_core(
            position_prefix, position_suffix, core_gene_positions
        )
        # assertion
        self.assertEqual(actual_new_positions, [1, 2, 5, 6])

    def test___replace_invalid_gene_positions_middle(self):
        # setup
        graph = GeneMerGraph({}, 3)
        gene_positions = [(0, 5), (None, None), (None, None), (16, 20), (None, None), (26, 30)]
        fastq_data = {"read1": {"sequence": "ACTGTACTGTACTGTACTGTACTGTACTGTACTGT"}}
        # execution
        actual_new_positions = graph.replace_invalid_gene_positions(
            gene_positions, fastq_data, "read1"
        )
        # assertion
        expected_new_positions = [(0, 5), (5, 16), (5, 16), (16, 20), (20, 26), (26, 30)]
        self.assertEqual(actual_new_positions, expected_new_positions)

    def test___replace_invalid_gene_positions_start(self):
        # setup
        graph = GeneMerGraph({}, 3)
        gene_positions = [(None, None), (6, 15), (16, 20), (20, 26), (26, 30)]
        fastq_data = {"read1": {"sequence": "ACTGTACTGTACTGTACTGTACTGTACTGTACTGT"}}
        # execution
        actual_new_positions = graph.replace_invalid_gene_positions(
            gene_positions, fastq_data, "read1"
        )
        # assertion
        expected_new_positions = [(0, 6), (6, 15), (16, 20), (20, 26), (26, 30)]
        self.assertEqual(actual_new_positions, expected_new_positions)

    def test___replace_invalid_gene_positions_end(self):
        # setup
        graph = GeneMerGraph({}, 3)
        gene_positions = [(0, 5), (6, 20), (11, 15), (16, 20), (21, 25), (None, None)]
        fastq_data = {"read1": {"sequence": "ACTGTACTGTACTGTACTGTACTGTACTGTACTGT"}}
        # execution
        actual_new_positions = graph.replace_invalid_gene_positions(
            gene_positions, fastq_data, "read1"
        )
        # assertion
        expected_new_positions = [(0, 5), (6, 20), (11, 15), (16, 20), (21, 25), (25, 34)]
        self.assertEqual(actual_new_positions, expected_new_positions)

    def test___split_into_subpaths_linear(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene5",
            "+gene9",
            "-gene10",
            "+gene11",
        ]
        positions1 = [
            (0, 1),
            (2, 3),
            (4, 5),
            (6, 7),
            (8, 9),
            (10, 11),
            (12, 13),
            (14, 15),
            (16, 17),
            (18, 19),
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes1},
            3,
            {"read1": positions1, "read2": positions1, "read3": positions1},
        )
        nodesOfInterest = graph.get_nodes_containing("gene5")
        nodeHashesOfInterest = [n.__hash__() for n in nodesOfInterest]
        reads_with_gene = graph.collect_reads_in_path(nodeHashesOfInterest)
        gene_call_subset = {r: graph.get_reads()[r] for r in reads_with_gene}
        rc_reads = {}
        for r in gene_call_subset:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(gene_call_subset[r])
        gene_call_subset.update(rc_reads)
        suffix_tree = construct_suffix_tree(graph.get_readNodes())
        paths, path_coverages = graph.get_paths_for_gene(
            suffix_tree, gene_call_subset, nodeHashesOfInterest, 1, "gene5", 1
        )

        # execution
        actual_final_paths, _ = graph.split_into_subpaths("gene5", paths, path_coverages)
        # assertion
        self.assertEqual(len(actual_final_paths), 2)
        for k in actual_final_paths:
            self.assertEqual(len(actual_final_paths[k]), 3)
        expected_1 = set(["read1_8_9", "read2_8_9", "read3_8_9"])
        expected_2 = set(["read1_12_13", "read2_12_13", "read3_12_13"])
        self.assertTrue(
            any(len(set(actual_final_paths[k]).intersection(expected_1))) == 3
            for k in actual_final_paths
        )
        self.assertTrue(
            any(len(set(actual_final_paths[k]).intersection(expected_2))) == 3
            for k in actual_final_paths
        )

    def test___split_into_subpaths_triangle(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene7",
            "+gene5",
            "+gene7",
            "+gene5",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10",
            "+gene11",
        ]
        positions1 = [
            (0, 1),
            (2, 3),
            (4, 5),
            (6, 7),
            (8, 9),
            (10, 11),
            (12, 13),
            (14, 15),
            (16, 17),
            (18, 19),
            (20, 21),
            (22, 23),
            (24, 25),
            (26, 27),
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes1},
            3,
            {"read1": positions1, "read2": positions1, "read3": positions1},
        )
        nodesOfInterest = graph.get_nodes_containing("gene5")
        nodeHashesOfInterest = [n.__hash__() for n in nodesOfInterest]
        reads_with_gene = graph.collect_reads_in_path(nodeHashesOfInterest)
        gene_call_subset = {r: graph.get_reads()[r] for r in reads_with_gene}
        rc_reads = {}
        for r in gene_call_subset:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(gene_call_subset[r])
        gene_call_subset.update(rc_reads)
        suffix_tree = construct_suffix_tree(graph.get_readNodes())
        paths, path_coverages = graph.get_paths_for_gene(
            suffix_tree, gene_call_subset, nodeHashesOfInterest, 1, "gene5", 1
        )
        # execution
        actual_final_paths, _ = graph.split_into_subpaths("gene5", paths, path_coverages)
        # assertion
        self.assertEqual(len(actual_final_paths), 3)
        for k in actual_final_paths:
            self.assertEqual(len(actual_final_paths[k]), 3)
        expected_1 = set(["read1_16_17", "read2_16_17", "read3_16_17"])
        expected_2 = set(["read1_12_13", "read2_12_13", "read3_12_13"])
        expected_3 = set(["read1_8_9", "read2_8_9", "read3_8_9"])
        self.assertTrue(
            any(len(set(actual_final_paths[k]).intersection(expected_1))) == 3
            for k in actual_final_paths
        )
        self.assertTrue(
            any(len(set(actual_final_paths[k]).intersection(expected_2))) == 3
            for k in actual_final_paths
        )
        self.assertTrue(
            any(len(set(actual_final_paths[k]).intersection(expected_3))) == 3
            for k in actual_final_paths
        )

    def test___assess_connectivity(self):
        # setup
        import sourmash

        seq1 = "ATGGTCTCCGAGCTGCAGCGCCAGCTGGCGCTGCATCGGCAGACCCGCGGTGTAGGGTCTTCGTCGACTGCTT"
        seq2 = "ATGGTCTCCGAGCTGCAGCGCCAGCTTTCGCTGCATCGGCAGACCCGCGGTGTAGGGTCTTCGTCGACTGCTT"
        seq3 = "ATGAGTAGTAGGTCGTCGATCGTCAGCTGGATCTGAGATTCGGATTCGGCGGCTATCGGCTAGTCGACTGCTT"
        m1 = sourmash.MinHash(n=0, ksize=9, scaled=1)
        m1.add_sequence(seq1, force=True)
        m2 = sourmash.MinHash(n=0, ksize=9, scaled=1)
        m2.add_sequence(seq2, force=True)
        m3 = sourmash.MinHash(n=0, ksize=9, scaled=1)
        m3.add_sequence(seq3, force=True)
        minhash_dictionary = {(1, 2, 3): m1, (1, 4, 3): m2, (1, 5, 3): m3}
        path_dictionary = {
            (1, 2, 3): ["read1", "read2", "read3"],
            (1, 4, 3): ["read4", "read5", "read6"],
            (1, 5, 3): ["read7", "read8", "read9"],
        }
        graph = GeneMerGraph({}, 3)
        # execution
        actual_connectivity = graph.assess_connectivity(path_dictionary, minhash_dictionary, 0.9)
        # assertion
        self.assertEqual(actual_connectivity[(1, 2, 3)], {(1, 4, 3)})
        self.assertEqual(actual_connectivity[(1, 4, 3)], {(1, 2, 3)})
        self.assertEqual(actual_connectivity[(1, 5, 3)], set())

    def test___assess_connectivity_1(self):
        # setup
        import sourmash

        seq1 = "ATGGTCTCCGAGCTGCAGCGCCAGCTGGCGCTGCATCGGCAGACCCGCGGTGTAGGGTCTTCGTCGACTGCTT"
        seq2 = "ATGGTCTCCGAGCTGCAGCGCCAGCTTTCGCTGCATCGGCAGACCCGCGGTGTAGGGTCTTCGTCGACTGCTT"
        seq3 = "ATGAGTAGTAGGTCGTCGATCGTCAGCTGGATCTGAGATTCGGATTCGGCGGCTATCGGCTAGTCGACTGCTT"
        m1 = sourmash.MinHash(n=0, ksize=9, scaled=1)
        m1.add_sequence(seq1, force=True)
        m2 = sourmash.MinHash(n=0, ksize=9, scaled=1)
        m2.add_sequence(seq2, force=True)
        m3 = sourmash.MinHash(n=0, ksize=9, scaled=1)
        m3.add_sequence(seq3, force=True)
        minhash_dictionary = {(1, 2, 3): m1, (1, 4, 3): m2, (1, 5, 3): m3}
        path_dictionary = {
            (1, 2, 3): ["read1", "read2", "read3"],
            (1, 4, 3): ["read4", "read5", "read6"],
            (1, 5, 3): ["read7", "read8", "read9"],
        }
        graph = GeneMerGraph({}, 3)
        # execution
        actual_connectivity = graph.assess_connectivity(path_dictionary, minhash_dictionary, 1)
        # assertion
        self.assertEqual(actual_connectivity[(1, 2, 3)], set())
        self.assertEqual(actual_connectivity[(1, 4, 3)], set())
        self.assertEqual(actual_connectivity[(1, 5, 3)], set())

    def test___assess_connectivity_zero(self):
        # setup
        import sourmash

        seq1 = "ATGGTCTCCGAGCTGCAGCGCCAGCTGGCGCTGCATCGGCAGACCCGCGGTGTAGGGTCTTCGTCGACTGCTT"
        seq2 = "ATGGTCTCCGAGCTGCAGCGCCAGCTTTCGCTGCATCGGCAGACCCGCGGTGTAGGGTCTTCGTCGACTGCTT"
        seq3 = "ATGAGTAGTAGGTCGTCGATCGTCAGCTGGATCTGAGATTCGGATTCGGCGGCTATCGGCTAGTCGACTGCTT"
        m1 = sourmash.MinHash(n=0, ksize=9, scaled=1)
        m1.add_sequence(seq1, force=True)
        m2 = sourmash.MinHash(n=0, ksize=9, scaled=1)
        m2.add_sequence(seq2, force=True)
        m3 = sourmash.MinHash(n=0, ksize=9, scaled=1)
        m3.add_sequence(seq3, force=True)
        minhash_dictionary = {(1, 2, 3): m1, (1, 4, 3): m2, (1, 5, 3): m3}
        path_dictionary = {
            (1, 2, 3): ["read1", "read2", "read3"],
            (1, 4, 3): ["read4", "read5", "read6"],
            (1, 5, 3): ["read7", "read8", "read9"],
        }
        graph = GeneMerGraph({}, 3)
        # execution
        actual_connectivity = graph.assess_connectivity(path_dictionary, minhash_dictionary, 0)
        # assertion
        self.assertEqual(actual_connectivity[(1, 2, 3)], {(1, 4, 3), (1, 5, 3)})
        self.assertEqual(actual_connectivity[(1, 4, 3)], {(1, 2, 3), (1, 5, 3)})
        self.assertEqual(actual_connectivity[(1, 5, 3)], {(1, 4, 3), (1, 2, 3)})

    def test___cluster_paths_one(self):
        # setup
        graph = GeneMerGraph({}, 3)
        # execution
        actual_merged_paths = graph.cluster_paths(
            {(1, 2, 3): {(1, 4, 3)}, (1, 4, 3): {(1, 2, 3)}, (1, 5, 3): set()}
        )
        # assertion
        self.assertEqual(
            actual_merged_paths, {(1, 2, 3): {(1, 4, 3), (1, 2, 3)}, (1, 5, 3): {(1, 5, 3)}}
        )

    def test___cluster_paths_two(self):
        # setup
        graph = GeneMerGraph({}, 3)
        # execution
        actual_merged_paths = graph.cluster_paths(
            {(1, 2, 3): set(), (1, 4, 3): set(), (1, 5, 3): set()}
        )
        # assertion
        self.assertEqual(
            actual_merged_paths,
            {(1, 2, 3): {(1, 2, 3)}, (1, 4, 3): {(1, 4, 3)}, (1, 5, 3): {(1, 5, 3)}},
        )

    def test___cluster_paths_three(self):
        # setup
        graph = GeneMerGraph({}, 3)
        # execution
        actual_merged_paths = graph.cluster_paths(
            {
                (1, 2, 3): {(1, 4, 3), (1, 5, 3)},
                (1, 4, 3): {(1, 2, 3), (1, 5, 3)},
                (1, 5, 3): {(1, 4, 3), (1, 2, 3)},
            }
        )
        # assertion
        self.assertEqual(actual_merged_paths, {(1, 2, 3): {(1, 4, 3), (1, 5, 3), (1, 2, 3)}})

    def test___find_sublist_indices_one(self):
        # setup
        graph = GeneMerGraph({}, 3)
        # execution
        actual_indices = graph.find_sublist_indices(
            ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"], ["4", "5", "6"]
        )
        # assertion
        self.assertEqual(actual_indices, [(3, 5)])

    def test___find_sublist_indices_two(self):
        # setup
        graph = GeneMerGraph({}, 3)
        # execution
        actual_indices = graph.find_sublist_indices(
            ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"], ["11", "12", "13"]
        )
        # assertion
        self.assertEqual(actual_indices, [])

    def test___find_sublist_indices_three(self):
        # setup
        graph = GeneMerGraph({}, 3)
        # execution
        actual_indices = graph.find_sublist_indices(
            ["1", "2", "3", "4", "5", "6", "2", "3", "4", "10"], ["2", "3", "4"]
        )
        # assertion
        self.assertEqual(actual_indices, [(1, 3), (6, 8)])

    def test___find_sublist_indices_four(self):
        # setup
        graph = GeneMerGraph({}, 3)
        # execution
        actual_indices = graph.find_sublist_indices(
            ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"],
            ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"],
        )
        # assertion
        self.assertEqual(actual_indices, [(0, 9)])

    def test___find_sublist_indices_five(self):
        # setup
        graph = GeneMerGraph({}, 3)
        # execution
        actual_indices = graph.find_sublist_indices(["1", "1", "1", "1", "1"], ["1", "1", "1"])
        # assertion
        self.assertEqual(actual_indices, [(0, 2), (1, 3), (2, 4)])

    def test___make_intersection_matrix(self):
        # setup
        annotations = {
            "read1": [
                "+gene1",
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
                "+gene8",
                "+gene9",
                "+gene10",
            ],
            "read2": ["-gene4", "+gene5", "+gene6", "+gene7", "+gene8", "+gene9", "+gene10"],
            "read3": [
                "+gene1",
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
                "+gene8",
            ],
            "read4": ["+gene3", "-gene4", "+gene5", "+gene6", "+gene7", "+gene8"],
            "read5": [
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
                "+gene8",
                "+gene9",
            ],
            "read6": ["+gene7", "+gene8", "+gene9", "+gene10"],
            "read7": ["+gene3", "-gene4", "+gene5", "+gene6", "+gene7"],
        }
        graph = GeneMerGraph(annotations, 3)
        actual_matrix, actual_node_hashes = graph.make_intersection_matrix()
        # assertion
        expected_matrix = [
            [2, 2, 2, 2, 2, 2, 1, 1],
            [2, 3, 3, 3, 3, 3, 2, 1],
            [2, 3, 5, 5, 5, 4, 2, 1],
            [2, 3, 5, 6, 6, 5, 3, 2],
            [2, 3, 5, 6, 6, 5, 3, 2],
            [2, 3, 4, 5, 5, 5, 3, 2],
            [1, 2, 2, 3, 3, 3, 4, 3],
            [1, 1, 1, 2, 2, 2, 3, 3],
        ]
        self.assertEqual(actual_matrix, expected_matrix)
        self.assertEqual(len(actual_node_hashes), 8)

    def test___trim_fringe_nodes_linear(self):
        # setup
        annotations = {
            "read1": [
                "+gene1",
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
                "+gene8",
                "+gene9",
                "+gene10",
            ],
            "read2": ["-gene4", "+gene5", "+gene6", "+gene7", "+gene8", "+gene9", "+gene10"],
            "read3": [
                "+gene1",
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
                "+gene8",
            ],
            "read4": ["+gene3", "-gene4", "+gene5", "+gene6", "+gene7", "+gene8"],
            "read5": [
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
                "+gene8",
                "+gene9",
            ],
            "read6": ["+gene7", "+gene8", "+gene9", "+gene10"],
            "read7": ["+gene3", "-gene4", "+gene5", "+gene6", "+gene7"],
        }
        graph = GeneMerGraph(annotations, 3)
        matrix, node_hashes = graph.make_intersection_matrix()
        # execution
        trimmed_graph = graph.trim_fringe_nodes(5, matrix, node_hashes)
        # assertion
        self.assertEqual(len(trimmed_graph.get_nodes()), 4)
        self.assertTrue(
            all(len(n.get_list_of_reads()) in {5, 6} for n in trimmed_graph.all_nodes())
        )

    def test___trim_fringe_nodes_circle(self):
        # setup
        annotations = {
            "read1": [
                "-gene0",
                "+gene1",
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
                "+gene8",
                "+gene9",
                "+gene10",
                "-gene11",
            ],
            "read2": [
                "-gene0",
                "+gene1",
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
                "+gene8",
                "+gene9",
                "+gene10",
                "-gene11",
                "+gene12",
            ],
            "read3": [
                "+gene1",
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
                "+gene8",
                "+gene9",
                "+gene10",
                "-gene11",
                "+gene12",
            ],
            "read4": [
                "-gene0",
                "+gene1",
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
                "+gene8",
                "+gene9",
                "+gene10",
                "-gene11",
            ],
            "read5": [
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
                "+gene8",
                "+gene9",
                "+gene10",
                "-gene11",
                "+gene12",
            ],
            "read6": [
                "-gene0",
                "+gene1",
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene13",
                "+gene14",
                "+gene15",
                "+gene8",
                "+gene9",
                "+gene10",
                "-gene11",
                "+gene12",
            ],
            "read7": [
                "+gene1",
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene13",
                "+gene14",
                "+gene15",
                "+gene8",
                "+gene9",
                "+gene10",
                "-gene11",
            ],
            "read8": [
                "+gene1",
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene13",
                "+gene14",
                "+gene15",
                "+gene8",
                "+gene9",
                "+gene10",
                "-gene11",
            ],
            "read9": [
                "+gene1",
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene13",
                "+gene14",
                "+gene15",
                "+gene8",
                "+gene9",
                "+gene10",
                "-gene11",
            ],
            "read10": [
                "+gene1",
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene13",
                "+gene14",
                "+gene15",
                "+gene8",
                "+gene9",
                "+gene10",
                "-gene11",
            ],
        }
        graph = GeneMerGraph(annotations, 3)
        matrix, node_hashes = graph.make_intersection_matrix()
        # execution
        trimmed_graph = graph.trim_fringe_nodes(5, matrix, node_hashes)
        # assertion
        self.assertEqual(len(trimmed_graph.get_nodes()), 14)
        self.assertTrue(
            all(len(n.get_list_of_reads()) in {9, 10, 5} for n in trimmed_graph.all_nodes())
        )

    def test___trim_fringe_nodes_junction(self):
        # setup
        annotations = {
            "read1": [
                "+gene1",
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "-gene4",
                "+gene5",
                "+gene6",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
                "+gene8",
                "+gene9",
                "+gene10",
            ],
            "read2": [
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "-gene4",
                "+gene5",
                "+gene6",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
                "+gene8",
                "+gene9",
                "+gene10",
            ],
            "read3": [
                "+gene1",
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "-gene4",
                "+gene5",
                "+gene6",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
                "+gene8",
            ],
            "read4": [
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "-gene4",
                "+gene5",
                "+gene6",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
                "+gene8",
            ],
            "read5": [
                "-gene2",
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "-gene4",
                "+gene5",
                "+gene6",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
                "+gene8",
                "+gene9",
            ],
            "read6": ["+gene7", "+gene8", "+gene9", "+gene10"],
            "read7": [
                "+gene3",
                "-gene4",
                "+gene5",
                "+gene6",
                "-gene4",
                "+gene5",
                "+gene6",
                "-gene4",
                "+gene5",
                "+gene6",
                "+gene7",
            ],
        }
        graph = GeneMerGraph(annotations, 3)
        matrix, node_hashes = graph.make_intersection_matrix()
        # execution
        trimmed_graph = graph.trim_fringe_nodes(5, matrix, node_hashes)
        # assertion
        self.assertEqual(len(trimmed_graph.get_nodes()), 6)
        self.assertTrue(
            all(len(n.get_list_of_reads()) in {5, 6} for n in trimmed_graph.all_nodes())
        )

    def test___trim_fringe_nodes_complex(self):
        # setup
        import json

        call_file = "tests/complex_gene_calls_one.json"
        position_file = "tests/complex_gene_positions_one.json"
        with open(call_file) as i:
            calls = json.load(i)
        with open(position_file) as i:
            positions = json.load(i)
        filtered_calls = {}
        for r in calls:
            if any(g[1:] == "mphANG_0479861" for g in calls[r]):
                filtered_calls[r] = calls[r]
        graph = GeneMerGraph(filtered_calls, 3, positions)
        matrix, node_hashes = graph.make_intersection_matrix()
        # execution
        trimmed_graph = graph.trim_fringe_nodes(5, matrix, node_hashes)
        # assertion
        self.assertEqual(len(trimmed_graph.get_nodes()), 66)

    def test___get_closest_allele(self):
        # setup
        samfile = "tests/test_allele.sam"
        graph = GeneMerGraph({}, 3, {})
        # execution
        actual_validity, actual_references, actual_unique_reads = graph.get_closest_allele(
            samfile, "allele"
        )
        # assertion
        self.assertTrue(actual_validity)
        self.assertEqual(len(actual_references), 6)

    def test___path_finding_between_junctions(self):
        # setup
        import json

        with open("tests/test_path_calls.json") as i:
            calls = json.load(i)
        graph = GeneMerGraph(calls, 3)
        graph.filter_graph(3, 1)
        # get the potential start nodes
        potential_bubble_starts = graph.identify_potential_bubble_starts()
        max_distance = graph.get_kmerSize() * 3
        for component in graph.components():
            if component not in potential_bubble_starts:
                continue
            potential_bubble_starts_component = potential_bubble_starts[component]
            unique_paths = graph.get_all_paths_between_junctions_in_component(
                potential_bubble_starts_component, max_distance, 1
            )
            filtered_paths = graph.filter_paths_between_bubble_starts(unique_paths)
            print(filtered_paths)
            sorted_filtered_paths = sorted(filtered_paths, key=lambda x: len(x[0]), reverse=True)
            # assertion
            # Amira will never correct a path that starts and ends at the same node
            self.assertEqual(len(sorted_filtered_paths), 2)

    def test___get_minhashes_for_paths_same_path(self):
        # setup
        import json

        with open("tests/test_path_calls.json") as i:
            calls = json.load(i)
        with open("tests/test_path_positions.json") as i:
            positions = json.load(i)
        graph = GeneMerGraph(calls, 3, positions)
        fastq_data = parse_fastq("tests/test_1.fastq.gz")
        # get the potential start nodes
        potential_bubble_starts = graph.identify_potential_bubble_starts()
        max_distance = graph.get_kmerSize() * 3
        for component in graph.components():
            if component not in potential_bubble_starts:
                continue
            potential_bubble_starts_component = potential_bubble_starts[component]
            unique_paths = graph.get_all_paths_between_junctions_in_component(
                potential_bubble_starts_component, max_distance, 1
            )
            filtered_paths = graph.filter_paths_between_bubble_starts(unique_paths)
            sorted_filtered_paths = sorted(filtered_paths, key=lambda x: len(x[0]), reverse=True)
            # execution
            path_minimizers = graph.get_minhashes_for_paths(sorted_filtered_paths, fastq_data, 1)
            min1 = path_minimizers[tuple([n[0] for n in sorted_filtered_paths[0][0]])]
            min2 = path_minimizers[tuple([n[0] for n in sorted_filtered_paths[1][0]])]
            # assertion
            self.assertEqual(len(min1 & min2) / len(min1), 0.9155718701700154)
            self.assertEqual(len(min1 & min2) / len(min2), 0.9052531041069724)

    def test___get_subpaths_long_collapsed(self):
        # setup
        import json

        with open("tests/complex_gene_calls_three.json") as i:
            calls = json.load(i)
        with open("tests/complex_gene_positions_three.json") as i:
            positions = json.load(i)
        graph = GeneMerGraph(calls, 3, positions)
        nodesOfInterest = []
        for geneOfInterest in [
            "mphANG_0479861",
        ]:
            nodesOfInterest += graph.get_nodes_containing(geneOfInterest)
        nodeHashesOfInterest = set([n.__hash__() for n in nodesOfInterest])
        reads_with_gene = graph.collect_reads_in_path(nodeHashesOfInterest)
        gene_call_subset = {r: graph.get_reads()[r] for r in reads_with_gene}
        rc_reads = {}
        for r in gene_call_subset:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(gene_call_subset[r])
        gene_call_subset.update(rc_reads)
        # execution
        suffix_tree = construct_suffix_tree(graph.get_readNodes())
        paths, path_coverages = graph.get_paths_for_gene(
            suffix_tree, gene_call_subset, nodeHashesOfInterest, 1, geneOfInterest, 1
        )
        # assertion
        for p in paths:
            print(p, paths[p])
        self.assertEqual(len(paths), 4)
        self.assertTrue(all(paths[p] in {151, 101, 129, 131} for p in paths))

    def test___correct_low_coverage_paths(self):
        # setup
        import json

        with open("tests/complex_gene_calls_four.json") as i:
            calls = json.load(i)
        with open("tests/complex_gene_positions_four.json") as i:
            positions = json.load(i)
        graph = GeneMerGraph(calls, 5, positions)
        # execution
        actual_bubble_starts = graph.identify_potential_bubble_starts()
        # assertion
        for component in actual_bubble_starts:
            self.assertEqual(len(actual_bubble_starts[component]), 4)
            # execution
            potential_bubble_starts_component = actual_bubble_starts[component]
            # get all the paths of length <= max distance between all pairs of junctions
            unique_paths = graph.get_all_paths_between_junctions_in_component(
                potential_bubble_starts_component, 15, 1
            )
            self.assertEqual(len(unique_paths), 2)
