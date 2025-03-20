import unittest

from amira.construct_graph import GeneMerGraph
from amira.path_finding_utils import cluster_upstream_adjacent_paths, construct_suffix_tree


class TestPathFindingConstructor(unittest.TestCase):

    def test___find_full_paths_linear_simple(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "-gene6",
            "+gene7",
        ]
        genes2 = [
            "-gene2",
            "+gene3",
            "-gene4",
        ]
        genes3 = [
            "-gene4",
            "-gene6",
            "+gene7",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes3}, 3
        )
        nodesOfInterest = graph.get_nodes_containing("gene4")
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "gene4",
            1,
        )
        # assertion
        self.assertEqual(len(full_blocks), 1)
        for k in full_blocks:
            self.assertEqual(full_blocks[k], 4)

    def test___find_full_paths_no_adjacent_paths(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene4",
            "-gene4",
            "-gene4",
            "+gene7",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes1, "read4": genes1}, 3
        )
        nodesOfInterest = graph.get_nodes_containing("gene4")
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "gene4",
            1,
        )
        # assertion
        self.assertEqual(len(full_blocks), 1)
        for k in full_blocks:
            self.assertEqual(full_blocks[k], 4)

    def test___find_full_paths_one_adjacent_path(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene4", "-gene4", "-gene4", "+gene7", "-gene8"]
        genes2 = [
            "+gene-1",
            "-gene0",
            "+gene1",
            "-gene2",
            "+gene4",
            "-gene4",
            "-gene4",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10" "+gene11",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes2}, 3
        )
        nodesOfInterest = graph.get_nodes_containing("gene4")
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "gene4",
            1,
        )
        # assertion
        self.assertEqual(len(full_blocks), 1)
        for k in full_blocks:
            self.assertEqual(full_blocks[k], 4)

    def test___find_full_paths_linear_path_duplicates_simple(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene4",
            "-gene6",
            "+gene7",
        ]
        genes2 = [
            "-gene2",
            "+gene3",
            "-gene4",
        ]
        genes3 = [
            "+gene4",
            "-gene6",
            "+gene7",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes3}, 3
        )
        nodesOfInterest = graph.get_nodes_containing("gene4")
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "gene4",
            1,
        )
        # assertion
        self.assertEqual(len(full_blocks), 1)
        for k in full_blocks:
            self.assertEqual(len(k), 6)
            self.assertEqual(full_blocks[k], 2)

    def test___find_full_paths_linear_path_contained(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene7", "-gene8"]
        genes2 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes2}, 3
        )
        nodesOfInterest = graph.get_nodes_containing("gene4")
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "gene4",
            1,
        )
        # assertion
        self.assertEqual(len(full_blocks), 2)
        for k in full_blocks:
            self.assertTrue(len(k) == 5 or len(k) == 8)
            self.assertEqual(full_blocks[k], 2)

    def test___find_full_paths_linear_path_contained_two(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
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
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene3",
            "-gene4",
            "+gene5",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes2}, 3
        )
        nodesOfInterest = graph.get_nodes_containing("gene4")
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "gene4",
            1,
        )
        # assertion
        self.assertEqual(len(full_blocks), 2)
        for k in full_blocks:
            self.assertTrue(len(k) == 13 or len(k) == 10)
            self.assertEqual(full_blocks[k], 2)

    def test___find_full_paths_terminate_at_junction(self):
        # setup
        genes1 = [
            "-trbC",
            "-trbB",
            "-group_1081",
            "-group_6156",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "-group_5175",
            "+group_5625",
        ]
        genes2 = [
            "-alkB",
            "-ada",
            "-apbE",
            "-ompC",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+rcsD",
            "+rcsB",
            "-rcsC",
            "+atoS",
            "+atoC",
            "+atoD",
            "+atoA",
            "+atoE",
            "+atoB",
            "-yfaP",
            "-yfaQ",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes2}, 3
        )
        nodesOfInterest = graph.get_nodes_containing("blaCMY54NG_0488491")
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "blaCMY54NG_0488491",
            1,
        )
        # assertion
        self.assertEqual(len(full_blocks), 2)
        for k in full_blocks:
            self.assertTrue(len(k) == 5 or len(k) == 20)
            self.assertEqual(full_blocks[k], 2)

    def test___find_full_paths_terminate_and_start_at_junction(self):
        # setup
        genes1 = [
            "-trbC",
            "-trbB",
            "-group_1081",
            "-group_6156",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-group_5175",
            "+group_5625",
        ]
        genes2 = [
            "-alkB",
            "-ada",
            "-apbE",
            "-ompC",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+rcsD",
            "+rcsB",
            "-rcsC",
            "+atoS",
            "+atoC",
            "+atoD",
            "+atoA",
            "+atoE",
            "+atoB",
            "-yfaP",
            "-yfaQ",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes2}, 3
        )
        nodesOfInterest = graph.get_nodes_containing("blaCMY54NG_0488491")
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "blaCMY54NG_0488491",
            1,
        )
        # assertion
        self.assertEqual(len(full_blocks), 2)
        for k in full_blocks:
            self.assertTrue(len(k) == 20 or len(k) == 5)
            self.assertEqual(full_blocks[k], 2)

    def test___find_full_paths_singleton(self):
        # setup
        genes1 = [
            "+gene7",
            "-gene4",
            "-gene13",
        ]
        graph = GeneMerGraph({"read1": genes1, "read2": genes1}, 3)
        nodesOfInterest = graph.get_nodes_containing("gene7")
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "gene7",
            1,
        )
        # assertion
        self.assertEqual(len(nodeAnchors), 1)
        self.assertEqual(len(full_blocks), 0)

    def test___get_singleton_paths(self):
        # setup
        genes1 = [
            "+gene7",
            "-gene4",
            "-gene13",
        ]
        graph = GeneMerGraph({"read1": genes1, "read2": genes1}, 3)
        nodesOfInterest = graph.get_nodes_containing("gene7")
        tree = construct_suffix_tree(graph.get_readNodes())
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        # execution
        full_blocks, seen_nodes, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "gene7",
            1,
        )
        graph.get_singleton_paths(seen_nodes, nodeAnchors, full_blocks, full_block_coverages)
        # assertion
        self.assertEqual(len(nodeAnchors), 1)
        self.assertEqual(len(full_blocks), 1)
        for f in full_blocks:
            self.assertEqual(full_blocks[f], 2)

    def test___find_full_paths_branching_path(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene5",
            "-gene6",
            "+gene7",
            "-gene4",
            "-gene6",
            "+gene7",
            "-gene10",
            "-gene11",
        ]
        genes2 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene5",
            "-gene6",
            "+gene7",
            "-gene4",
            "-gene13",
            "+gene14",
            "-gene15",
            "-gene16",
        ]
        genes3 = [
            "+gene7",
            "-gene4",
            "-gene13",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes2, "read5": genes3}, 3
        )
        nodesOfInterest = graph.get_nodes_containing("gene7")
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            2,
            calls,
            "gene7",
            1,
        )
        # assertion
        self.assertEqual(len(full_blocks), 2)
        for k in full_blocks:
            self.assertTrue(len(k) == 8 or len(k) == 3)
            self.assertTrue(full_blocks[k] == 2 or full_blocks[k] == 3)

    def test___find_full_paths_triangle(self):
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
        genes2 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10",
            "+gene11",
        ]
        genes3 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene7",
            "+gene5",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10",
            "+gene11",
        ]
        graph = GeneMerGraph(
            {
                "read1": genes1,
                "read2": genes1,
                "read3": genes2,
                "read4": genes2,
                "read5": genes3,
                "read6": genes3,
            },
            3,
        )
        nodesOfInterest = graph.get_nodes_containing("gene5")
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "gene5",
            1,
        )
        # assertion
        self.assertEqual(len(full_blocks), 3)
        for k in full_blocks:
            self.assertTrue(len(k) == 5 or len(k) == 7 or len(k) == 9)
            self.assertEqual(full_blocks[k], 2)

    def test___find_full_paths_linear_path_duplicates_long_reads(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene9",
            "+gene4",
            "-gene10",
            "+gene11",
            "-gene12",
        ]
        genes2 = [
            "-gene2",
            "+gene3",
            "-gene4",
        ]
        genes3 = [
            "+gene4",
            "-gene10",
            "+gene11",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes3}, 3
        )
        nodesOfInterest = graph.get_nodes_containing("gene4")
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "gene4",
            1,
        )
        # assertion
        self.assertEqual(len(full_blocks), 1)
        for k in full_blocks:
            self.assertEqual(len(k), 11)
            self.assertEqual(full_blocks[k], 2)

    def test___find_full_paths_complex_one(self):
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
        nodesOfInterest = graph.get_nodes_containing("mphANG_0479861")
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "mphANG_0479861",
            1,
        )
        # assertion
        self.assertEqual(len(full_blocks), 2)
        for k in full_blocks:
            self.assertTrue(len(k) == 16 or len(k) == 6)
            self.assertTrue(full_blocks[k] == 7 or full_blocks[k] == 52)

    def test___find_full_paths_diverging_paths_at_terminals(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10",
            "+gene11",
            "-gene12",
            "+gene13",
            "-gene14",
            "+gene15",
        ]
        genes2 = [
            "+gene16",
            "-gene17",
            "+gene18",
            "-gene19",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10",
            "+gene11",
            "-gene20",
            "+gene21",
            "-gene22",
            "+gene23",
        ]
        graph = GeneMerGraph(
            {
                "read1": genes1,
                "read2": genes1,
                "read3": genes1,
                "read4": genes2,
                "read5": genes2,
                "read6": genes2,
            },
            3,
        )
        nodesOfInterest = graph.get_nodes_containing("gene8")
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "gene8",
            1,
        )
        # assertion
        self.assertEqual(len(full_blocks), 2)
        for k in full_blocks:
            self.assertEqual(len(k), 9)
            self.assertEqual(full_blocks[k], 3)

    def test___find_full_paths_edge_case_one(self):
        # setup
        import json

        with open("tests/complex_gene_calls_five.json") as i:
            calls = json.load(i)
        with open("tests/complex_gene_positions_five.json") as i:
            positions = json.load(i)
        graph = GeneMerGraph(calls, 3, positions)
        nodesOfInterest = []
        for geneOfInterest in [
            "blaCTXM110NG_0489052",
        ]:
            nodesOfInterest += graph.get_nodes_containing(geneOfInterest)
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "blaCTXM110NG_0489052",
            1,
        )
        # assertion
        self.assertEqual(len(full_blocks), 1)
        for f in full_blocks:
            self.assertEqual(len(f), 1)
            self.assertEqual(full_blocks[f], 44)

    def test___find_full_paths_variant(self):
        # setup
        import json

        with open("tests/complex_gene_calls_six.json") as i:
            calls = json.load(i)
        with open("tests/complex_gene_positions_six.json") as i:
            positions = json.load(i)
        graph = GeneMerGraph(calls, 3, positions)
        nodesOfInterest = []
        for geneOfInterest in ["blaTEM239NG_0766451"]:
            nodesOfInterest += graph.get_nodes_containing(geneOfInterest)
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "blaTEM239NG_0766451",
            1,
        )
        # assertion
        self.assertEqual(len(full_blocks), 2)
        for f in full_blocks:
            self.assertEqual(len(f), 2)
            self.assertTrue(full_blocks[f] == 30 or full_blocks[f] == 9)

    def test___find_full_paths_multi_tandem(self):
        # setup
        read1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene9",
        ]
        read2 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene9",
        ]
        read3 = [
            "-gene2",
            "+gene3",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "+gene5",
            "-gene6",
        ]
        read4 = [
            "+gene3",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "+gene5",
        ]
        read5 = ["+gene1", "-gene2", "+gene3", "-gene4", "-gene4", "-gene4"]
        read6 = ["-gene4", "-gene4", "-gene4", "-gene4", "-gene4", "+gene5", "-gene6"]
        read7 = ["-gene4", "-gene4", "-gene4"]
        read8 = ["+gene3", "-gene4", "-gene4", "-gene4", "-gene4", "-gene4", "-gene4", "+gene5"]
        read9 = ["-gene10", "+gene4", "-gene11"]
        read10 = [
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
        ]
        calls = {
            "read1": read1,
            "read2": read2,
            "read3": GeneMerGraph({}, 3).reverse_list_of_genes(read3),
            "read4": GeneMerGraph({}, 3).reverse_list_of_genes(read4),
            "read5": read5,
            "read6": read6,
            "read7": read7,
            "read8": read8,
            "read9": read9,
            "read10": read10,
            "read11": read10,
            "read12": read10,
            "read13": read10,
        }
        graph = GeneMerGraph(calls, 3)
        nodesOfInterest = []
        for geneOfInterest in [
            "gene4",
        ]:
            nodesOfInterest += graph.get_nodes_containing(geneOfInterest)
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "gene4",
            1,
        )
        self.assertEqual(len(full_blocks), 2)
        for f in full_blocks:
            self.assertTrue(len(f) == 9 or len(f) == 8)
            self.assertTrue(full_blocks[f] == 6 or full_blocks[f] == 3)

    def test___find_full_paths_end_with_self_loop(self):
        # setup
        import json

        with open("tests/complex_gene_calls_seven.json") as i:
            calls = json.load(i)
        with open("tests/complex_gene_positions_seven.json") as i:
            positions = json.load(i)
        graph = GeneMerGraph(calls, 3, positions)
        nodesOfInterest = []
        for geneOfInterest in [
            "blaIMI9NG_0491711",
        ]:
            nodesOfInterest += graph.get_nodes_containing(geneOfInterest)
        all_reads = graph.collect_reads_in_path([n.__hash__() for n in nodesOfInterest])
        tree = construct_suffix_tree({r: graph.get_readNodes()[r] for r in all_reads})
        calls = {r: graph.get_reads()[r] for r in all_reads}
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            3,
            calls,
            "blaIMI9NG_0491711",
            1,
        )
        # assertion
        self.assertEqual(len(full_blocks), 1)
        for f in full_blocks:
            self.assertEqual(len(f), 5)
            self.assertEqual(full_blocks[f], 4)

    def test___cluster_adjacent_paths(self):
        # setup
        adjacent_paths = {
            (0, 1, 2, 3, 4): {"read1"},
            (1, 2, 3, 4): {"read2"},
            (2, 3, 4): {"read3"},
            (5, 6, 3, 4): {"read4", "read5"},
            (6, 3, 4): {"read6"},
            (5, 3, 2, 4): {"read7"},
            (3, 4): {"read8"},
        }
        # execution
        actual_clusters = cluster_upstream_adjacent_paths(adjacent_paths)
        # assertion
        self.assertEqual(len(actual_clusters), 3)
        self.assertTrue((2, 3, 4) in actual_clusters)
        self.assertTrue((6, 3, 4) in actual_clusters)
        self.assertTrue((5, 3, 2, 4) in actual_clusters)

    def test___cluster_adjacent_paths_overlapping(self):
        # setup
        adjacent_paths = {
            (0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 12): {"read1", "read2"},
            (5, 1, 2, 3, 4, 7, 8, 9, 10, 11, 12): {"read3", "read4", "read5"},
            (5, 6, 2, 3, 4, 7, 8, 9, 10, 11, 12): {"read6", "read7"},
            (2, 3, 4, 7, 8, 9, 10): {"read8"},
        }
        # execution
        actual_clusters = cluster_upstream_adjacent_paths(adjacent_paths)
        # assertion
        self.assertEqual(len(actual_clusters), 4)
        self.assertTrue((0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 12) in actual_clusters)
        self.assertTrue((5, 1, 2, 3, 4, 7, 8, 9, 10, 11, 12) in actual_clusters)
        self.assertTrue((5, 6, 2, 3, 4, 7, 8, 9, 10, 11, 12) in actual_clusters)
        self.assertTrue((2, 3, 4, 7, 8, 9, 10) in actual_clusters)

    def test___partial_internal_block(self):
        # setup
        import json

        with open("tests/complex_gene_calls_eight.json") as i:
            calls = json.load(i)
        with open("tests/complex_gene_positions_eight.json") as i:
            positions = json.load(i)
        graph = GeneMerGraph(calls, 3, positions)
        nodesOfInterest = []
        for geneOfInterest in ["dfrA17"]:
            nodesOfInterest += graph.get_nodes_containing(geneOfInterest)
        tree = construct_suffix_tree(graph.get_readNodes())
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        node_mapping = {n.__hash__(): n for n in nodesOfInterest}
        nodeAnchors = graph.get_AMR_anchors(node_mapping)
        # execution
        full_blocks, _, full_block_coverages = graph.get_full_paths(
            tree,
            graph.get_readNodes(),
            nodeAnchors,
            1,
            calls,
            "dfrA17",
            1,
        )
        finalAllelesOfInterest, copy_numbers = graph.split_into_subpaths(
            "dfrA17", full_blocks, full_block_coverages, 188, {}
        )
        # assertion
        self.assertEqual(len(full_blocks), 3)
        self.assertEqual(len(finalAllelesOfInterest), 2)
        for f in full_blocks:
            self.assertTrue(len(f) == 6 or len(f) == 3)
        for f in finalAllelesOfInterest:
            self.assertTrue(
                len(finalAllelesOfInterest[f]) == 69 or len(finalAllelesOfInterest[f]) == 52
            )
