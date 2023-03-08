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
                        ["gene2"])
        AMRNodes = unitig.get_all_nodes_containing_AMR_genes()
        # execution
        actual_targetNode, actual_targetDirection = unitig.get_forward_node_from_node(nodes[0],
                                                                                    set(AMRNodes.keys()))
        actual_targetNode = actual_targetNode[0]
        actual_targetDirection = actual_targetDirection[0]
        # assertion
        expected_targetNode = nodes[1]
        expected_targetDirection = 1
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
                        ["gene2"])
        AMRNodes = unitig.get_all_nodes_containing_AMR_genes()
        # execution
        actual_targetNode, actual_targetDirection = unitig.get_forward_node_from_node(nodes[0],
                                                                                    set(AMRNodes.keys()))
        # assertion
        expected_targetNode = None
        expected_targetDirection = None
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
                        ["gene2"])
        AMRNodes = unitig.get_all_nodes_containing_AMR_genes()
        # execution
        actual_targetNode, actual_targetDirection = unitig.get_backward_node_from_node(nodes[0],
                                                                                    set(AMRNodes.keys()))
        actual_targetNode = actual_targetNode[0]
        actual_targetDirection = actual_targetDirection[0]
        # assertion
        expected_targetNode = nodes[1]
        expected_targetDirection = -1
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
                        ["gene2"])
        AMRNodes = unitig.get_all_nodes_containing_AMR_genes()
        # execution
        actual_targetNode, actual_targetDirection = unitig.get_backward_node_from_node(nodes[0],
                                                                                    set(AMRNodes.keys()))
        # assertion
        expected_targetNode = None
        expected_targetDirection = None
        self.assertEqual(actual_targetNode, expected_targetNode)
        self.assertEqual(actual_targetDirection, expected_targetDirection)

   # def test___get_unitigs_of_interest(self):
    #    # setup
     #   genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
      #  genes2 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
       # graph = GeneMerGraph({"read1": genes1, "read2": genes2},
       #                     3)
       # unitig = Unitigs(graph,
      #                  ["gene1","gene4", "gene7"])
        # execution
       # actual_unitigGenesOfInterest, actual_unitigsReadsOfInterest = unitig.get_unitigs_of_interest()
       # actual_unitig_count = len(actual_unitigGenesOfInterest)
       # actual_unitig_read_count = len(actual_unitigsReadsOfInterest)
        # assertion
       # expected_unitig_count = 7
       # expected_unitig_read_count = 7
       # self.assertEqual(actual_unitig_count, expected_unitig_count)
       # self.assertEqual(actual_unitig_read_count, expected_unitig_read_count)

 #   def test___visualise_unitigs(self):
        # setup
#        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
#        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
#        graph = GeneMerGraph({"read1": genes1, "read2": genes2},
#                            3)
#        unitig = Unitigs(graph,
##                        ["gene1","gene4", "gene7"])
#        unitigTools = UnitigTools(graph,
#                                ["gene1","gene4", "gene7"])
#        unitig_mappings = {}
#        unitigs = unitigTools.get_unitigGenesOfInterest()
#        count = 1
#        for u in unitigs:
#            unitig_mappings[u] = count
#            count += 1
#        # execution
#        unitigTools.visualise_unitigs({"gene1": [1500], "gene2": [3000], "gene3": [1000], "gene4": [90], "gene5": [600], "gene6": [2000], "gene7": [100], "gene8": [1000], "gene9": [4000], "gene10": [400]},
#                            unitig_mappings,
#                            "tests")
 #       # assertion
#        import os
 #       self.assertTrue(os.path.exists("tests/context_plots.pdf"))
 #       os.remove("tests/context_plots.pdf")

    def test___get_all_nodes_containing_AMR_genes(self):
        # setup
        genes1 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3},
                            3)
        unitig = Unitigs(graph,
                        ["gene7", "gene4", "gene1"])
        # execution
        actual_AMRNodes = unitig.get_all_nodes_containing_AMR_genes()
        # assertion
        expected_numberAMRNodes = 15
        self.assertEqual(len(actual_AMRNodes), expected_numberAMRNodes)

    def test___get_AMR_anchors(self):
        # setup
        genes1 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3},
                            3)
        unitig = Unitigs(graph,
                        ["gene7", "gene4", "gene1"])
        AMRNodes = unitig.get_all_nodes_containing_AMR_genes()
        # execution
        actual_AMRanchors = unitig.get_AMR_anchors(AMRNodes)
        # assertion
        expected_anchor_count = 6
        self.assertEqual(len(actual_AMRanchors), expected_anchor_count)

    def test___get_unitigs_of_interest(self):
        # setup
        genes1 = ["-gene6", "+gene10", "+gene9", "-gene6", "+gene3", "-gene7", "+gene5", "-gene6", "+gene3", "-gene7", "-gene6", "+gene3", "-gene7", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5", "+gene3", "-gene4", "+gene5"]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        genes4 = ["+gene7", "-gene2", "+gene0"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3, "read4": genes4},
                            3)
        unitig = Unitigs(graph,
                        ["gene7", "gene4", "gene1"])
        # execution
        graph.generate_gml("test",
                    3,
                    1,
                    1)
        actual_unitigGenesOfInterest, actual_unitigsReadsOfInterest, actual_uniqueReads = unitig.get_unitigs_of_interest()
        actual_unitig_count = len(actual_unitigGenesOfInterest)
        actual_unitig_read_count = len(actual_unitigsReadsOfInterest)
        # assertion
        expected_unitig_count = 2
        expected_unitig_read_count = 7
        self.assertEqual(actual_unitig_count, expected_unitig_count)
        self.assertEqual(actual_unitig_read_count, expected_unitig_read_count)