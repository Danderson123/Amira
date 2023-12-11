import unittest

from amira_prototype.construct_graph import GeneMerGraph
from amira_prototype.construct_unitig import UnitigTools


class TestUnitigsConstructor(unittest.TestCase):
    def test___init_Unitigs(self):
        # setup
        graph = GeneMerGraph(
            {"read1": ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]}, 3
        )
        # execution
        actual_unitigs = UnitigTools(graph, ["gene4"], "tests/test.fastq.gz", ".")
        actual_graph = actual_unitigs.get_graph()
        actual_selected_genes = actual_unitigs.get_selected_genes()
        # assertion
        expected_graph = graph
        expected_selected_genes = ["gene4"]
        self.assertEqual(actual_graph, expected_graph)
        self.assertEqual(actual_selected_genes, expected_selected_genes)

    def test___get_nodes_of_interest(self):
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
        unitigs = UnitigTools(graph, [], "tests/test.fastq.gz", ".")
        # execution
        actual_nodes_of_interest = unitigs.get_nodes_of_interest("gene4")
        # assertion
        expected_nodes_of_interest = graph.get_nodes_containing("gene4")
        self.assertNotEqual(actual_nodes_of_interest, [])
        self.assertEqual(actual_nodes_of_interest, expected_nodes_of_interest)

    def test___get_all_nodes_containing_AMR_genes(self):
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
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3}, 3)
        unitig = UnitigTools(graph, ["gene7", "gene4", "gene1"], "tests/test.fastq.gz", ".")
        # execution
        actual_AMRNodes = unitig.get_all_nodes_containing_AMR_genes()
        # assertion
        expected_numberAMRNodes = 15
        self.assertEqual(len(actual_AMRNodes), expected_numberAMRNodes)

    def test___cluster_anchor_reads_simple(self):
        # setup
        graph = GeneMerGraph({}, 3)
        unitig = UnitigTools(graph, [], "tests/test.fastq.gz", ".")
        reads_dict = {
            "read1": ["node1"],
            "read2": ["node1", "node2"],
            "read3": ["node1"],
            "read4": ["node2"],
            "read5": ["node3"],
        }
        # execution
        actual_clusters = unitig.cluster_anchor_reads(reads_dict)
        # assertion
        expected_clusters = {1: {"read1", "read2", "read3", "read4"}, 2: {"read5"}}
        self.assertEqual(actual_clusters, expected_clusters)

    def test___cluster_anchor_reads_later_occurrence(self):
        # setup
        graph = GeneMerGraph({}, 3)
        unitig = UnitigTools(graph, [], "tests/test.fastq.gz", ".")
        reads_dict = {
            "read1": ["node1"],
            "read2": ["node1", "node2"],
            "read3": ["node1"],
            "read4": ["node2", "node5"],
            "read5": ["node3", "node5"],
            "read7": ["node4"],
            "read8": ["node5"],
        }
        # execution
        actual_clusters = unitig.cluster_anchor_reads(reads_dict)
        # assertion
        expected_clusters = {
            1: {"read1", "read2", "read3", "read4", "read5", "read8"},
            2: {"read7"},
        }
        self.assertEqual(actual_clusters, expected_clusters)

    def test___get_AMR_anchors_and_junctions(self):
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
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3}, 3)
        unitig = UnitigTools(graph, ["gene7", "gene4", "gene1"], "tests/test.fastq.gz", ".")
        AMRNodes = unitig.get_all_nodes_containing_AMR_genes()
        # execution
        actual_AMRanchors, actual_AMRjunctions = unitig.get_AMR_anchors_and_junctions(AMRNodes)
        # assertion
        expected_anchor_count = 4
        expected_junction_count = 2
        self.assertEqual(len(actual_AMRanchors), expected_anchor_count)
        self.assertEqual(len(actual_AMRjunctions), expected_junction_count)

    def test___get_anchoring_reads_one(self):
        # setup
        genes1 = [
            "+gene9",
            "-gene6",
            "+gene7",
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
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3}, 3)
        unitig = UnitigTools(graph, ["gene7"], "tests/test.fastq.gz", ".")
        AMRNodes = unitig.get_all_nodes_containing_AMR_genes()
        nodeAnchors, nodeJunctions = unitig.get_AMR_anchors_and_junctions(AMRNodes)
        # execution
        actual_anchorReads = unitig.get_anchoring_reads(nodeAnchors)
        # assertion
        self.assertEqual(len(AMRNodes), 10)
        self.assertEqual(len(nodeAnchors), 4)
        self.assertEqual(len(nodeJunctions), 1)
        self.assertEqual(len(actual_anchorReads), 2)
        self.assertEqual(len(set(actual_anchorReads["read1"])), 4)
        self.assertEqual(len(set(actual_anchorReads["read2"])), 2)

    def test___get_anchoring_reads_two(self):
        # setup
        genes1 = [
            "+gene9",
            "-gene6",
            "+gene7",
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
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3}, 3)
        unitig = UnitigTools(graph, ["gene1"], "tests/test.fastq.gz", ".")
        AMRNodes = unitig.get_all_nodes_containing_AMR_genes()
        nodeAnchors, nodeJunctions = unitig.get_AMR_anchors_and_junctions(AMRNodes)
        # execution
        actual_anchorReads = unitig.get_anchoring_reads(nodeAnchors)
        # assertion
        self.assertEqual(len(AMRNodes), 2)
        self.assertEqual(len(nodeAnchors), 2)
        self.assertEqual(len(nodeJunctions), 0)
        self.assertEqual(len(actual_anchorReads), 1)
        self.assertEqual(len(set(actual_anchorReads["read3"])), 2)

    def test___get_anchoring_reads_three(self):
        # setup
        genes1 = [
            "+gene9",
            "-gene6",
            "+gene7",
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
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3}, 3)
        unitig = UnitigTools(graph, ["gene4"], "tests/test.fastq.gz", ".")
        AMRNodes = unitig.get_all_nodes_containing_AMR_genes()
        nodeAnchors, nodeJunctions = unitig.get_AMR_anchors_and_junctions(AMRNodes)
        # execution
        actual_anchorReads = unitig.get_anchoring_reads(nodeAnchors)
        # assertion
        self.assertEqual(len(AMRNodes), 6)
        self.assertEqual(len(nodeAnchors), 3)
        self.assertEqual(len(nodeJunctions), 1)
        self.assertEqual(len(actual_anchorReads), 1)
        self.assertEqual(len(set(actual_anchorReads["read1"])), 3)

    def test___get_anchoring_reads_multiple_genes(self):
        # setup
        genes1 = [
            "+gene9",
            "-gene6",
            "+gene7",
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
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3}, 3)
        unitig = UnitigTools(graph, ["gene1", "gene4", "gene7"], "tests/test.fastq.gz", ".")
        AMRNodes = unitig.get_all_nodes_containing_AMR_genes()
        nodeAnchors, nodeJunctions = unitig.get_AMR_anchors_and_junctions(AMRNodes)
        # execution
        actual_anchorReads = unitig.get_anchoring_reads(nodeAnchors)
        # assertion
        self.assertEqual(len(AMRNodes), 16)
        self.assertEqual(len(nodeAnchors), 5)
        self.assertEqual(len(nodeJunctions), 2)
        self.assertEqual(len(actual_anchorReads), 2)
        self.assertEqual(len(set(actual_anchorReads["read1"])), 3)
        self.assertEqual(len(set(actual_anchorReads["read3"])), 2)

    def test___assess_resolvability(self):
        # setup
        # setup
        genes1 = [
            "+gene9",
            "-gene6",
            "+gene7",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene10",
            "+gene9",
            "+gene5",
            "-gene6",
            "+gene3",
            "-gene7",
            "-gene6",
            "+gene3",
            "-gene7",
        ]
        genes2 = ["+gene9", "-gene6", "+gene7", "+gene3", "-gene4", "+gene5"]
        genes3 = ["-gene0", "+gene1", "-gene2", "+gene3"]
        graph = GeneMerGraph({"read1": genes1, "read2": genes2, "read3": genes3}, 3)
        unitig = UnitigTools(graph, ["gene1", "gene4", "gene7"], "tests/test.fastq.gz", ".")
        AMRNodes = unitig.get_all_nodes_containing_AMR_genes()
        nodeAnchors, nodeJunctions = unitig.get_AMR_anchors_and_junctions(AMRNodes)
        anchorReads = unitig.get_anchoring_reads(nodeAnchors)
        clusters = unitig.cluster_anchor_reads(anchorReads)
        # execution
        actual_easy, actual_intermediate, actual_difficult, clusterId = unitig.assess_resolvability(
            clusters, nodeJunctions
        )
        # assertion
        self.assertEqual(len(AMRNodes), 10)
        self.assertEqual(len(nodeAnchors), 4)
        self.assertEqual(len(nodeJunctions), 1)
        self.assertEqual(len(anchorReads), 2)
        self.assertEqual(len(clusters), 2)
        self.assertEqual(len(actual_easy), 1)
        self.assertEqual(len(actual_intermediate), 1)
        self.assertEqual(len(actual_difficult), 0)

    def test___contains_sublist_contains(self):
        # setup
        graph = GeneMerGraph({}, 3)
        unitig = UnitigTools(graph, [], "tests/test.fastq.gz", ".")
        allNodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        # assertion
        self.assertTrue(unitig.contains_sublist(allNodes, [1, 2, 3, 4]))
        self.assertTrue(unitig.contains_sublist(allNodes, [4, 3, 2, 1]))
        self.assertTrue(unitig.contains_sublist(allNodes, [5, 6, 7, 8]))
        self.assertTrue(unitig.contains_sublist(allNodes, [8, 7, 6, 5]))
        self.assertTrue(unitig.contains_sublist(allNodes, [7, 8, 9, 10]))
        self.assertTrue(unitig.contains_sublist(allNodes, [10, 9, 8, 7]))
        self.assertTrue(unitig.contains_sublist(allNodes, [1]))
        self.assertTrue(unitig.contains_sublist(allNodes, [10]))

    def test___contains_sublist_does_not_contain(self):
        # setup
        graph = GeneMerGraph({}, 3)
        unitig = UnitigTools(graph, [], "tests/test.fastq.gz", ".")
        allNodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        # assertion
        self.assertFalse(unitig.contains_sublist(allNodes, [11, 12, 13, 14]))
        self.assertFalse(unitig.contains_sublist(allNodes, [6, 7, 11, 12]))
        self.assertFalse(unitig.contains_sublist(allNodes, [5, 11, 7, 6]))

    def test___contains_sublist_overlaps(self):
        # setup
        graph = GeneMerGraph({}, 3)
        unitig = UnitigTools(graph, [], "tests/test.fastq.gz", ".")
        allNodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        # assertion
        self.assertFalse(unitig.contains_sublist(allNodes, [10, 11, 12, 13]))
        self.assertFalse(unitig.contains_sublist(allNodes, [9, 10, 11, 12]))
        self.assertFalse(unitig.contains_sublist(allNodes, [8, 9, 10, 11]))
        self.assertFalse(unitig.contains_sublist(allNodes, [3, 2, 1, 0]))
        self.assertFalse(unitig.contains_sublist(allNodes, [2, 1, 0]))
        self.assertFalse(unitig.contains_sublist(allNodes, [1, 0]))

    def test_basic(self):
        # setup
        graph = GeneMerGraph({}, 3)
        unitig = UnitigTools(graph, [], "tests/test.fastq.gz", ".")
        # assertion
        self.assertEqual(
            unitig.get_unitigs_on_read([1, 2, 3, 4, 5, 6, 7, 8, 9], {4, 7}), [[4, 5, 6, 7]]
        )
        self.assertEqual(
            unitig.get_unitigs_on_read([9, 8, 7, 6, 5, 4, 3, 2, 1], {4, 7}), [[7, 6, 5, 4]]
        )

    def test_advanced(self):
        # setup
        graph = GeneMerGraph({}, 3)
        unitig = UnitigTools(graph, [], "tests/test.fastq.gz", ".")
        # assertion
        self.assertEqual(
            unitig.get_unitigs_on_read([0, 1, 2, 3, 4, 5], {1, 2, 5}), [[1, 2], [2, 3, 4, 5]]
        )
        self.assertEqual(
            unitig.get_unitigs_on_read([5, 4, 3, 2, 1, 0], {1, 2, 5}), [[5, 4, 3, 2], [2, 1]]
        )

    def test_no_delimiters(self):
        # setup
        graph = GeneMerGraph({}, 3)
        unitig = UnitigTools(graph, [], "tests/test.fastq.gz", ".")
        # assertion
        self.assertEqual(unitig.get_unitigs_on_read([1, 2, 3], {4, 7}), [])
        self.assertEqual(unitig.get_unitigs_on_read([3, 2, 1], {4, 7}), [])

    def test_all_delimiters(self):
        # setup
        graph = GeneMerGraph({}, 3)
        unitig = UnitigTools(graph, [], "tests/test.fastq.gz", ".")
        # assertion
        self.assertEqual(unitig.get_unitigs_on_read([4, 7, 4, 7], {4, 7}), [[4, 7], [7, 4], [4, 7]])
        self.assertEqual(unitig.get_unitigs_on_read([7, 4, 7, 4], {4, 7}), [[7, 4], [4, 7], [7, 4]])

    def test_single_element(self):
        # setup
        graph = GeneMerGraph({}, 3)
        unitig = UnitigTools(graph, [], "tests/test.fastq.gz", ".")
        # assertion
        self.assertEqual(unitig.get_unitigs_on_read([4], {4}), [[4]])

    def test_empty_list(self):
        # setup
        graph = GeneMerGraph({}, 3)
        unitig = UnitigTools(graph, [], "tests/test.fastq.gz", ".")
        # assertion
        self.assertEqual(unitig.get_unitigs_on_read([], {4, 7}), [])
