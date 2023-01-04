import unittest
import sys
sys.path.insert(0, "amira_prototype")

from construct_graph import GeneMerGraph
from construct_unitig import Unitigs

class TestUnitigsConstructor(unittest.TestCase):

    def test___init_Unitigs(self):
        # setup
        graph = GeneMerGraph({"read1": ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]},
                            3,
                            1,
                            1)
        # execution
        actual_unitigs = Unitigs(graph,
                                ["gene4"])
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
                            3,
                            1,
                            1)
        unitigs = Unitigs(graph,
                        ["gene4"])
        # execution
        actual_nodes_of_interest = unitigs.get_nodes_of_interest()
        # assertion
        expected_nodes_of_interest = graph.get_nodes_containing(["gene4"])
        self.assertNotEqual(actual_nodes_of_interest, [])
        self.assertEqual(actual_nodes_of_interest, expected_nodes_of_interest)