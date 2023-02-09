import unittest
import sys
sys.path.insert(0, "amira_prototype")

from construct_graph import GeneMerGraph
from construct_unitig import Unitigs, UnitigTools
from construct_edge import Edge
from construct_read import Read
from construct_gene_mer import GeneMer
from construct_gene import Gene

class TestUnitigsConstructor(unittest.TestCase):

    def test___init_Unitigs(self):
        # setup
        graph = GeneMerGraph({"read1": ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]},
                            3)
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
                            3)
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
                        3)
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
                        3)
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
                        3)
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
                        3)
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
                            3)
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
        expected_forwardPath = [[[Gene("-gene6"), Gene("+gene7"), Gene("-gene8")],
                                [Gene("+gene7"), Gene("-gene8"), Gene("+gene9")]]]
        expected_forwardPath.append(list(reversed([list(reversed(e)) for e in expected_forwardPath[0]])))
        expected_forwardPath.append([[g.reverse_gene() for g in gmer] for gmer in expected_forwardPath[0]])
        expected_forwardPath.append([[g.reverse_gene() for g in gmer] for gmer in expected_forwardPath[1]])
        self.assertTrue(any(actual_forwardPath[0] == e for e in expected_forwardPath))

    def test___get_forward_path_from_first_node(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4"]
        genes2 = ["-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2},
                            3)
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
        expected_forwardPath = [[[Gene("+gene1"), Gene("-gene2"), Gene("+gene3")],
                                [Gene("-gene2"), Gene("+gene3"), Gene("-gene4")],
                                [Gene("+gene3"), Gene("-gene4"), Gene("+gene5")],
                                [Gene("-gene4"), Gene("+gene5"), Gene("-gene6")]]]
        expected_forwardPath.append(list(reversed([list(reversed(e)) for e in expected_forwardPath[0]])))
        expected_forwardPath.append([[g.reverse_gene() for g in gmer] for gmer in expected_forwardPath[0]])
        expected_forwardPath.append([[g.reverse_gene() for g in gmer] for gmer in expected_forwardPath[1]])
        self.assertTrue(any(actual_forwardPath[0] == e for e in expected_forwardPath))

    def test___get_backward_path_from_node(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4"]
        genes2 = ["-gene6", "+gene7", "-gene8", "+gene9"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2},
                            3)
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
        unitig = Unitigs(graph,
                        ["gene5"])
        # execution
        actual_backwardPath = unitig.get_backward_path_from_node(targetNode)
        # assertion
        expected_backwardPath = [[[Gene("+gene1"), Gene("-gene2"), Gene("+gene3")],
                                [Gene("-gene2"), Gene("+gene3"), Gene("-gene4")]]]
        expected_backwardPath.append(list(reversed([list(reversed(e)) for e in expected_backwardPath[0]])))
        expected_backwardPath.append([[g.reverse_gene() for g in gmer] for gmer in expected_backwardPath[0]])
        expected_backwardPath.append([[g.reverse_gene() for g in gmer] for gmer in expected_backwardPath[1]])
        self.assertTrue(any(actual_backwardPath[0] == e for e in expected_backwardPath))

    def test___get_backward_path_from_final_node(self):
        # setup
        genes1 = ["-gene4", "+gene5", "-gene6", "+gene7"]
        genes2 = ["+gene5","-gene6", "+gene7", "-gene8", "+gene9"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2},
                            3)
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
        expected_backwardPath = [[[Gene("-gene4"), Gene("+gene5"), Gene("-gene6")],
                                [Gene("+gene5"), Gene("-gene6"), Gene("+gene7")],
                                [Gene("-gene6"), Gene("+gene7"), Gene("-gene8")],
                                [Gene("+gene7"), Gene("-gene8"), Gene("+gene9")]]]
        expected_backwardPath.append(list(reversed([list(reversed(e)) for e in expected_backwardPath[0]])))
        expected_backwardPath.append([[g.reverse_gene() for g in gmer] for gmer in expected_backwardPath[0]])
        expected_backwardPath.append([[g.reverse_gene() for g in gmer] for gmer in expected_backwardPath[1]])
        self.assertTrue(any(actual_backwardPath[0] == e for e in expected_backwardPath))

    def test___get_unitig_for_node(self):
        # setup
        genes1 = ["-gene4", "+gene5", "-gene6", "+gene7"]
        genes2 = ["+gene5","-gene6", "+gene7", "-gene8", "+gene9"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2},
                            3)
        unitig = Unitigs(graph,
                        ["gene9"])
        node = unitig.get_nodes_of_interest("gene9")[0]
        # execution
        actual_unitig, actual_reads = unitig.get_unitig_for_node(node)
        # assertion
        expected_unitigs = [[[Gene("-gene4"), Gene("+gene5"), Gene("-gene6")],
                                [Gene("+gene5"), Gene("-gene6"), Gene("+gene7")],
                                [Gene("-gene6"), Gene("+gene7"), Gene("-gene8")],
                                [Gene("+gene7"), Gene("-gene8"), Gene("+gene9")]],
                            [[Gene("-gene9"), Gene("+gene8"), Gene("-gene7")],
                                [Gene("+gene8"), Gene("-gene7"), Gene("+gene6")],
                                [Gene("-gene7"), Gene("+gene6"), Gene("-gene5")],
                                [Gene("+gene6"), Gene("-gene5"), Gene("+gene4")]],
                            [[Gene("+gene4"), Gene("-gene5"), Gene("+gene6")],
                                [Gene("-gene5"), Gene("+gene6"), Gene("-gene7")],
                                [Gene("+gene6"), Gene("-gene7"), Gene("+gene8")],
                                [Gene("-gene7"), Gene("+gene8"), Gene("-gene9")]],
                            [[Gene("+gene9"), Gene("-gene8"), Gene("+gene7")],
                                [Gene("-gene8"), Gene("+gene7"), Gene("-gene6")],
                                [Gene("+gene7"), Gene("-gene6"), Gene("+gene5")],
                                [Gene("-gene6"), Gene("+gene5"), Gene("-gene4")]]]
        expected_reads =  [['read1'], ['read1', 'read2'], ['read2'], ['read2']]
        self.assertTrue(any(actual_unitig == expected_unitig for expected_unitig in expected_unitigs))
        self.assertTrue(actual_reads == expected_reads or actual_reads == list(reversed(expected_reads)))

    def test___hash_unitig(self):
        # setup
        genes1 = ["-gene4", "+gene5", "-gene6", "+gene7"]
        genes2 = ["+gene5","-gene6", "+gene7", "-gene8", "+gene9"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2},
                            3)
        unitig = Unitigs(graph,
                        ["gene9"])
        node = unitig.get_nodes_of_interest("gene9")[0]
        returned_unitig, returned_reads = unitig.get_unitig_for_node(node)
        # execution
        actual_unitig_hash = unitig.hash_unitig(returned_unitig,
                                                list(reversed(returned_unitig)))
        actual_reversed_unitig_hash =  unitig.hash_unitig(list(reversed(returned_unitig)),
                                                        returned_unitig)
        # assertion
        expected_geneMer_hash = [hash(tuple([gene.__hash__() for gene in geneMer])) for geneMer in returned_unitig]
        expected_reversed_geneMer_hash = [hash(tuple([gene.__hash__() for gene in geneMer])) for geneMer in list(reversed(returned_unitig))]
        expected_unitig_hash = sorted([hash(tuple(expected_geneMer_hash)), hash(tuple(expected_reversed_geneMer_hash))])[0]
        self.assertEqual(actual_unitig_hash, actual_reversed_unitig_hash)
        self.assertEqual(actual_unitig_hash, expected_unitig_hash)
        self.assertEqual(actual_reversed_unitig_hash, expected_unitig_hash)

    def test___get_unitigs_of_interest(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2},
                            3)
        unitig = Unitigs(graph,
                        ["gene1","gene4", "gene7"])
        # execution
        actual_unitigGenesOfInterest, actual_unitigsReadsOfInterest = unitig.get_unitigs_of_interest()
        actual_unitig_count = len(actual_unitigGenesOfInterest)
        actual_unitig_read_count = len(actual_unitigsReadsOfInterest)
        # assertion
        expected_unitig_count = 7
        expected_unitig_read_count = 7
        self.assertEqual(actual_unitig_count, expected_unitig_count)
        self.assertEqual(actual_unitig_read_count, expected_unitig_read_count)

    def test___visualise_unitigs(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2},
                            3)
        unitig = Unitigs(graph,
                        ["gene1","gene4", "gene7"])
        unitigTools = UnitigTools(graph,
                                ["gene1","gene4", "gene7"])
        unitig_mappings = {}
        unitigs = unitigTools.get_unitigGenesOfInterest()
        count = 1
        for u in unitigs:
            unitig_mappings[u] = count
            count += 1
        # execution
        unitigTools.visualise_unitigs({"gene1": [1500], "gene2": [3000], "gene3": [1000], "gene4": [90], "gene5": [600], "gene6": [2000], "gene7": [100], "gene8": [1000], "gene9": [4000], "gene10": [400]},
                            unitig_mappings,
                            "tests")
        # assertion
        import os
        self.assertTrue(os.path.exists("tests/context_plots.pdf"))
        os.remove("tests/context_plots.pdf")