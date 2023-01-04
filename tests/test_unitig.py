import unittest
import sys
sys.path.insert(0, "amira_prototype")

from construct_graph import GeneMerGraph
from construct_unitig import Unitigs
from construct_edge import Edge
from construct_read import Read

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

    def test___get_existing_forward_node_from_node(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3,
                        1,
                        1)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g)
            nodes.append(node)
        mock_forward_edge = Edge(nodes[0],
                                nodes[1],
                                1,
                                1)
        mock_rc_forward_edge = Edge(nodes[0],
                                    nodes[1],
                                    -1,
                                    -1)
        graph.add_edges_to_graph(mock_forward_edge,
                                mock_rc_forward_edge)
        graph.add_edge_to_node(nodes[0],
                            mock_forward_edge)
        graph.add_edge_to_node(nodes[1],
                            mock_rc_forward_edge)
        unitig = Unitigs(graph,
                        [])
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = unitig.get_forward_node_from_node(nodes[0])
        # assertion
        expected_extend = True
        expected_targetNode = nodes[1]
        expected_targetDirection = 1
        self.assertEqual(actual_extend, expected_extend)
        self.assertEqual(actual_targetNode, expected_targetNode)
        self.assertEqual(actual_targetDirection, expected_targetDirection)

    def test___get_non_existing_forward_node_from_node(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3,
                        1,
                        1)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g)
            nodes.append(node)
        unitig = Unitigs(graph,
                        [])
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = unitig.get_forward_node_from_node(nodes[0])
        # assertion
        expected_extend = False
        expected_targetNode = None
        expected_targetDirection = None
        self.assertEqual(actual_extend, expected_extend)
        self.assertEqual(actual_targetNode, expected_targetNode)
        self.assertEqual(actual_targetDirection, expected_targetDirection)

    def test___get_existing_backward_node_from_node(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3,
                        1,
                        1)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g)
            nodes.append(node)
        mock_backward_edge = Edge(nodes[0],
                                nodes[1],
                                -1,
                                -1)
        mock_rc_backward_edge = Edge(nodes[0],
                                    nodes[1],
                                    1,
                                    1)
        graph.add_edges_to_graph(mock_backward_edge,
                                mock_rc_backward_edge)
        graph.add_edge_to_node(nodes[0],
                            mock_backward_edge)
        graph.add_edge_to_node(nodes[1],
                            mock_rc_backward_edge)
        unitig = Unitigs(graph,
                        [])
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = unitig.get_backward_node_from_node(nodes[0])
        # assertion
        expected_extend = True
        expected_targetNode = nodes[1]
        expected_targetDirection = -1
        self.assertEqual(actual_extend, expected_extend)
        self.assertEqual(actual_targetNode, expected_targetNode)
        self.assertEqual(actual_targetDirection, expected_targetDirection)

    def test___get_non_existing_backward_node_from_node(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        graph = GeneMerGraph({},
                        3,
                        1,
                        1)
        nodes = []
        for g in geneMers:
            node = graph.add_node(g)
            nodes.append(node)
        unitig = Unitigs(graph,
                        [])
        # execution
        actual_extend, actual_targetNode, actual_targetDirection = unitig.get_backward_node_from_node(nodes[0])
        # assertion
        expected_extend = False
        expected_targetNode = None
        expected_targetDirection = None
        self.assertEqual(actual_extend, expected_extend)
        self.assertEqual(actual_targetNode, expected_targetNode)
        self.assertEqual(actual_targetDirection, expected_targetDirection)