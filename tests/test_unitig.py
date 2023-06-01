import unittest
import sys
sys.path.insert(0, "amira_prototype")

from construct_graph import GeneMerGraph
from construct_unitig import UnitigTools

class TestUnitigsConstructor(unittest.TestCase):

    def test___init_Unitigs(self):
        # setup
        graph = GeneMerGraph({"read1": ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]},
                            3)
        # execution
        actual_unitigs = UnitigTools(graph,
                                ["gene4"],
                                "tests/test.fastq.gz",
                                ".")
        actual_graph = actual_unitigs.get_graph()
        actual_selected_genes = actual_unitigs.get_selected_genes()
        # assertion
        expected_graph = graph
        expected_selected_genes = ["gene4"]
        self.assertEqual(actual_graph, expected_graph)
        self.assertEqual(actual_selected_genes, expected_selected_genes)

    def test___get_nodes_of_interest(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2},
                            3)
        unitigs = UnitigTools(graph,
                        [],
                        "tests/test.fastq.gz",
                        ".")
        # execution
        actual_nodes_of_interest = unitigs.get_nodes_of_interest("gene4")
        # assertion
        expected_nodes_of_interest = graph.get_nodes_containing("gene4")
        self.assertNotEqual(actual_nodes_of_interest, [])
        self.assertEqual(actual_nodes_of_interest, expected_nodes_of_interest)

    def test___get_all_nodes_containing_AMR_genes(self):
        # setup
        genes1 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3},
                            3)
        unitig = UnitigTools(graph,
                        ["gene7", "gene4", "gene1"],
                        "tests/test.fastq.gz",
                        ".")
        # execution
        actual_AMRNodes = unitig.get_all_nodes_containing_AMR_genes()
        # assertion
        expected_numberAMRNodes = 15
        self.assertEqual(len(actual_AMRNodes), expected_numberAMRNodes)

    # def test___get_reads_per_amr_node(self):
    #     # setup
    #     genes1 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
    #     genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
    #     genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
    #     graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3},
    #                         3)
    #     unitig = UnitigTools(graph,
    #                     ["gene7", "gene4", "gene1"],
    #                     "tests/test.fastq.gz",
    #                     ".")
    #     AMRNodes = unitig.get_all_nodes_containing_AMR_genes()
    #     # execution
    #     actual_reads_to_gene = unitig.get_reads_per_amr_node(AMRNodes)
    #     actual_reads = []
    #     for key, value in actual_reads_to_gene.items():
    #         actual_reads.append(value)
    #     # assertion
    #     expected_number_of_nodes = 15
    #     expected_reads = [['read1'],
    #                     ['read1'],
    #                     ['read1'],
    #                     ['read1'],
    #                     ['read1'],
    #                     ['read1'],
    #                     ['read1'],
    #                     ['read2'],
    #                     ['read2'],
    #                     ['read2'],
    #                     ['read1', 'read2'],
    #                     ['read1'], ['read1'],
    #                     ['read3'], ['read3']]
    #     self.assertEqual(len(actual_reads_to_gene), expected_number_of_nodes)
        # self.assertEqual(actual_reads, expected_reads)

    def test___cluster_reads_simple(self):
        # setup
        graph = GeneMerGraph({},
                            3)
        unitig = UnitigTools(graph,
                        [],
                        "tests/test.fastq.gz",
                        ".")
        reads_dict = {
            'read1': ['node1'],
            'read2': ['node1', 'node2'],
            'read3': ['node1'],
            'read4': ['node2'],
            'read5': ['node3']
        }
        # execution
        actual_clusters = unitig.cluster_anchor_reads(reads_dict)
        # assertion
        expected_clusters = {1: {'read1', 'read2', 'read3', 'read4'}, 2: {'read5'}}
        self.assertEqual(actual_clusters, expected_clusters)

    def test___cluster_reads_later_occurrence(self):
        # setup
        graph = GeneMerGraph({},
                            3)
        unitig = UnitigTools(graph,
                        [],
                        "tests/test.fastq.gz",
                        ".")
        reads_dict = {
            'read1': ['node1'],
            'read2': ['node1', 'node2'],
            'read3': ['node1'],
            'read4': ['node2', 'node5'],
            'read5': ['node3', 'node5'],
            'read7': ['node4'],
            'read8': ['node5']
        }
        # execution
        actual_clusters = unitig.cluster_anchor_reads(reads_dict)
        # assertion
        expected_clusters = {1: {'read1', 'read2', 'read3', 'read4', 'read5', 'read8'}, 2: {'read7'}}
        self.assertEqual(actual_clusters, expected_clusters)

    def test___get_AMR_anchors_and_junctions(self):
        # setup
        genes1 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3},
                            3)
        unitig = UnitigTools(graph,
                        ["gene7", "gene4", "gene1"],
                        "tests/test.fastq.gz",
                        ".")
        AMRNodes = unitig.get_all_nodes_containing_AMR_genes()
        # execution
        actual_AMRanchors, actual_AMRjunctions = unitig.get_AMR_anchors_and_junctions(AMRNodes)
        graph.generate_gml("test/test.gml",
                        3,
                        1,
                        1)
        # assertion
        expected_anchor_count = 4
        expected_junction_count = 2
        self.assertEqual(len(actual_AMRanchors), expected_anchor_count)
        self.assertEqual(len(actual_AMRjunctions), expected_junction_count)
