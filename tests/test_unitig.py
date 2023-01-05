import unittest
import sys
sys.path.insert(0, "amira_prototype")

from construct_graph import GeneMerGraph
from construct_unitig import Unitigs
from construct_edge import Edge
from construct_read import Read
from construct_node import Node
from construct_gene_mer import GeneMer
from construct_gene import Gene

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
                        [])
        # execution
        actual_nodes_of_interest = unitigs.get_nodes_of_interest("gene4")
        # assertion
        expected_nodes_of_interest = graph.get_nodes_containing("gene4")
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

    def test___get_forward_path_from_node(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4"]
        genes2 = ["-gene6", "+gene7", "-gene8", "+gene9"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2},
                            3,
                            1,
                            1)
        sourceGeneMer = GeneMer([Gene("+gene5"), Gene("-gene6"), Gene("+gene7")])
        sourceNode = graph.add_node(sourceGeneMer)
        sourceNode.increment_node_coverage()
        targetNode = graph.get_node(GeneMer([Gene("-gene6"), Gene("+gene7"), Gene("-gene8")]))
        mock_forward_edge = Edge(sourceNode,
                                targetNode,
                                1,
                                GeneMer([Gene("-gene6"), Gene("+gene7"), Gene("-gene8")]).get_geneMerDirection())
        mock_rc_forward_edge = Edge(targetNode,
                                    sourceNode,
                                    GeneMer([Gene("-gene6"), Gene("+gene7"), Gene("-gene8")]).get_geneMerDirection() * -1,
                                    -1)
        mock_forward_edge.increment_edge_coverage()
        mock_rc_forward_edge.increment_edge_coverage()
        graph.add_edges_to_graph(mock_forward_edge,
                                mock_rc_forward_edge)
        graph.add_edge_to_node(sourceNode,
                            mock_forward_edge)
        graph.add_edge_to_node(targetNode,
                            mock_rc_forward_edge)
        unitig = Unitigs(graph,
                        ["gene5"])
        # execution
        actual_forwardPath = unitig.get_forward_path_from_node(sourceNode)
        # assertion
        expected_forwardPath = [Node(GeneMer([Gene("-gene6"), Gene("+gene7"), Gene("-gene8")])).get_canonical_geneMer(),
                            Node(GeneMer([Gene("+gene7"), Gene("-gene8"), Gene("+gene9")])).get_canonical_geneMer()]
        self.assertEqual(actual_forwardPath, expected_forwardPath)

    def test___get_forward_path_from_first_node(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4"]
        genes2 = ["-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2},
                            3,
                            1,
                            1)
        sourceGeneMer = GeneMer([Gene("-gene0"), Gene("+gene1"), Gene("-gene2")])
        sourceNode = graph.add_node(sourceGeneMer)
        sourceNode.increment_node_coverage()
        targetNode = graph.get_node(GeneMer([Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]))
        mock_forward_edge = Edge(sourceNode,
                                targetNode,
                                1,
                                GeneMer([Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]).get_geneMerDirection())
        mock_rc_forward_edge = Edge(targetNode,
                                    sourceNode,
                                    GeneMer([Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]).get_geneMerDirection() * -1,
                                    -1)
        mock_forward_edge.increment_edge_coverage()
        mock_rc_forward_edge.increment_edge_coverage()
        graph.add_edges_to_graph(mock_forward_edge,
                                mock_rc_forward_edge)
        graph.add_edge_to_node(sourceNode,
                            mock_forward_edge)
        graph.add_edge_to_node(targetNode,
                            mock_rc_forward_edge)
        unitig = Unitigs(graph,
                        ["gene0"])
        # execution
        actual_forwardPath = unitig.get_forward_path_from_node(sourceNode)
        # assertion
        expected_forwardPath = [Node(GeneMer([Gene("+gene1"), Gene("-gene2"), Gene("+gene3")])).get_canonical_geneMer(),
                                Node(GeneMer([Gene("-gene2"), Gene("+gene3"), Gene("-gene4")])).get_canonical_geneMer(),
                                Node(GeneMer([Gene("+gene3"), Gene("-gene4"), Gene("+gene5")])).get_canonical_geneMer(),
                                Node(GeneMer([Gene("-gene4"), Gene("+gene5"), Gene("-gene6")])).get_canonical_geneMer()]
        self.assertEqual(actual_forwardPath, expected_forwardPath)

    def test___get_backward_path_from_node(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4"]
        genes2 = ["-gene6", "+gene7", "-gene8", "+gene9"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2},
                            3,
                            1,
                            1)
        targetGeneMer = GeneMer([Gene("+gene3"), Gene("-gene4"), Gene("+gene5")])
        targetNode = graph.add_node(targetGeneMer)
        targetNode.increment_node_coverage()
        sourceNode = graph.get_node(GeneMer([Gene("-gene2"), Gene("+gene3"), Gene("-gene4")]))
        mock_backward_edge = Edge(sourceNode,
                                targetNode,
                                GeneMer([Gene("-gene2"), Gene("+gene3"), Gene("-gene4")]).get_geneMerDirection(),
                                1)
        mock_rc_backward_edge = Edge(targetNode,
                                    sourceNode,
                                    -1,
                                    GeneMer([Gene("-gene2"), Gene("+gene3"), Gene("-gene4")]).get_geneMerDirection() * -1)
        mock_backward_edge.increment_edge_coverage()
        mock_rc_backward_edge.increment_edge_coverage()
        graph.add_edges_to_graph(mock_backward_edge,
                                mock_rc_backward_edge)
        graph.add_edge_to_node(sourceNode,
                            mock_backward_edge)
        graph.add_edge_to_node(targetNode,
                            mock_rc_backward_edge)
        graph.generate_gml("after")
        unitig = Unitigs(graph,
                        ["gene5"])
        # execution
        actual_backwardPath = unitig.get_backward_path_from_node(targetNode)
        # assertion
        expected_backwardPath = [Node(GeneMer([Gene("+gene1"), Gene("-gene2"), Gene("+gene3")])).get_canonical_geneMer(),
                                Node(GeneMer([Gene("-gene2"), Gene("+gene3"), Gene("-gene4")])).get_canonical_geneMer()]
        self.assertEqual(actual_backwardPath, expected_backwardPath)

    def test___get_backward_path_from_final_node(self):
        # setup
        genes1 = ["-gene4", "+gene5", "-gene6", "+gene7"]
        genes2 = ["+gene5","-gene6", "+gene7", "-gene8", "+gene9"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2},
                            3,
                            1,
                            1)
        unitig = Unitigs(graph,
                        ["gene9"])
        targetGeneMer = GeneMer([Gene("-gene8"), Gene("+gene9"), Gene("-gene10")])
        targetNode = graph.add_node(targetGeneMer)
        targetNode.increment_node_coverage()
        sourceNode = graph.get_node(GeneMer([Gene("+gene7"), Gene("-gene8"), Gene("+gene9")]))
        mock_backward_edge = Edge(sourceNode,
                                targetNode,
                                GeneMer([Gene("+gene7"), Gene("-gene8"), Gene("+gene9")]).get_geneMerDirection(),
                                1)
        mock_rc_backward_edge = Edge(targetNode,
                                    sourceNode,
                                    -1,
                                    GeneMer([Gene("+gene7"), Gene("-gene8"), Gene("+gene9")]).get_geneMerDirection() * -1)
        mock_backward_edge.increment_edge_coverage()
        mock_rc_backward_edge.increment_edge_coverage()
        graph.add_edges_to_graph(mock_backward_edge,
                                mock_rc_backward_edge)
        graph.add_edge_to_node(sourceNode,
                            mock_backward_edge)
        graph.add_edge_to_node(targetNode,
                            mock_rc_backward_edge)
        # execution
        actual_backwardPath = unitig.get_backward_path_from_node(targetNode)
        # assertion
        expected_backwardPath = [Node(GeneMer([Gene("-gene4"), Gene("+gene5"), Gene("-gene6")])).get_canonical_geneMer(),
                                Node(GeneMer([Gene("+gene5"), Gene("-gene6"), Gene("+gene7")])).get_canonical_geneMer(),
                                Node(GeneMer([Gene("-gene6"), Gene("+gene7"), Gene("-gene8")])).get_canonical_geneMer(),
                                Node(GeneMer([Gene("+gene7"), Gene("-gene8"), Gene("+gene9")])).get_canonical_geneMer()]
        self.assertEqual(actual_backwardPath, expected_backwardPath)

    def test___get_unitigs_of_interest(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2},
                            3,
                            1,
                            1)
        unitig = Unitigs(graph,
                        ["gene1","gene4", "gene7"])
        # execution
        actual_unitigs = unitig.get_unitigs_of_interest()
        actual_gene1_count = len(actual_unitigs["gene1"])
        actual_gene4_count = len(actual_unitigs["gene4"])
        actual_gene7_count = len(actual_unitigs["gene7"])
        # assertion
        expected_gene1_count = 1
        expected_gene4_count = 5
        expected_gene7_count = 4
        self.assertEqual(actual_gene1_count, expected_gene1_count)
        self.assertEqual(actual_gene4_count, expected_gene4_count)
        self.assertEqual(actual_gene7_count, expected_gene7_count)