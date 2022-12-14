import unittest
import sys
sys.path.insert(0, "amira_prototype")

from construct_graph import GeneMerGraph
from construct_read import Read

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
            graph.add_node(g,
                        "read1")
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
        actual_returned_node = graph.add_node(geneMer,
                                            "read1")
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
        graph.add_node(geneMer1,
                    "read1")
        genes2 = ["+gene4", "-gene3", "+gene2"]
        read2 = Read("read2",
                    genes2)
        geneMer2 = [x for x in read2.get_geneMers(3)][0]
        # execution
        actual_returned_node2 = graph.add_node(geneMer2,
                                            "read2")
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
        graph.add_node(geneMer,
                    "read1")
        # execution
        actual_returned_node = graph.add_node(geneMer,
                                            "read2")
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
        graph.add_node(geneMer1,
                    "read1")
        genes2 = ["+gene4", "-gene3", "+gene2"]
        read2 = Read("read2",
                    genes2)
        geneMer2 = [x for x in read2.get_geneMers(3)][0]
        graph.add_node(geneMer2,
                    "read2")
        # execution
        actual_node1 = graph.get_node(geneMer1)
        actual_node1_read_list = [x for x in actual_node1.get_reads()]
        actual_node1_hash = actual_node1.__hash__()
        actual_node1_coverage = actual_node1.get_node_coverage()
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
        graph.add_node(geneMer1,
                    "read1")
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
            graph.add_node(g,
                        "read1")
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
            graph.add_node(g,
                        "read1")
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
            graph.add_node(g,
                        "read1")
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
            graph.add_node(g,
                        "read1")
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
            graph.add_node(g,
                        "read1")
        selectedGenes = ["gene2", "+gene6"]
        # assertion
        self.assertRaises(AssertionError, graph.get_nodes_containing, selectedGenes)
