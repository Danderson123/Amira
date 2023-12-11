import unittest

from amira_prototype.construct_gene import Gene
from amira_prototype.construct_read import Read, convert_genes


class TestReadConstructor(unittest.TestCase):
    def test___init_Read(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
        # execution
        actual_read = Read("read1", genes)
        actual_readId = actual_read.get_readId()
        actual_numberOfGenes = actual_read.get_number_of_genes()
        actual_listOfGenes = actual_read.get_genes()
        # assertion
        expected_readId = "read1"
        expected_numberOfGenes = 6
        expected_listOfGenes = [Gene(g) for g in genes]
        self.assertEqual(actual_readId, expected_readId)
        self.assertEqual(actual_numberOfGenes, expected_numberOfGenes)
        self.assertEqual(actual_listOfGenes, expected_listOfGenes)

    def test___convert_non_empty_genes(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
        # execution
        actual_converted_genes = convert_genes(genes)
        # assertion
        expected_converted_genes = [Gene(g) for g in genes]
        self.assertEqual(actual_converted_genes, expected_converted_genes)

    def test___convert_empty_genes(self):
        # setup
        genes = []
        # execution
        actual_converted_genes = convert_genes(genes)
        # assertion
        expected_converted_genes = []
        self.assertEqual(actual_converted_genes, expected_converted_genes)

    def test_gene_geneMers_size_1(self):
        # setup
        annotatedGenes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5"]
        read = Read("read1", annotatedGenes)
        # execution
        actual_geneMers = [x for x in read.get_geneMers(1)]
        actual_geneObjects = [gmer.get_canonical_geneMer() for gmer in actual_geneMers]
        actual_geneStrings = [[g.get_name() for g in gmer] for gmer in actual_geneObjects]
        actual_number_of_geneMers = len(actual_geneMers)
        # assertion
        expected_number_of_geneMers = len(annotatedGenes)
        self.assertEqual(actual_number_of_geneMers, expected_number_of_geneMers)
        self.assertTrue(
            actual_geneStrings == [["gene1"], ["gene2"], ["gene3"], ["gene4"], ["gene5"]]
            or actual_geneStrings
            == list(reversed([["gene1"], ["gene2"], ["gene3"], ["gene4"], ["gene5"]]))
        )

    def test_gene_geneMers_size_2(self):
        # setup
        annotatedGenes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5"]
        read = Read("read1", annotatedGenes)
        # execution
        actual_geneMers = [x for x in read.get_geneMers(2)]
        actual_geneObjects = [gmer.get_canonical_geneMer() for gmer in actual_geneMers]
        actual_geneStrings = [[g.get_name() for g in gmer] for gmer in actual_geneObjects]
        actual_number_of_geneMers = len(actual_geneMers)
        # assertion
        expected_number_of_geneMers = len(annotatedGenes) - (2 - 1)
        self.assertEqual(actual_number_of_geneMers, expected_number_of_geneMers)
        self.assertTrue(
            all(any(test_gene[1:] in g for g in actual_geneStrings) for test_gene in annotatedGenes)
        )

    def test_gene_geneMers_size_3(self):
        # setup
        annotatedGenes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5"]
        read = Read("read1", annotatedGenes)
        # execution
        actual_geneMers = [x for x in read.get_geneMers(3)]
        actual_geneObjects = [gmer.get_canonical_geneMer() for gmer in actual_geneMers]
        actual_geneStrings = [[g.get_name() for g in gmer] for gmer in actual_geneObjects]
        actual_number_of_geneMers = len(actual_geneMers)
        # assertion
        expected_number_of_geneMers = len(annotatedGenes) - (3 - 1)
        self.assertEqual(actual_number_of_geneMers, expected_number_of_geneMers)
        self.assertTrue(
            all(any(test_gene[1:] in g for g in actual_geneStrings) for test_gene in annotatedGenes)
        )

    def test_gene_geneMers_size_4(self):
        # setup
        annotatedGenes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5"]
        read = Read("read1", annotatedGenes)
        # execution
        actual_geneMers = [x for x in read.get_geneMers(4)]
        actual_geneObjects = [gmer.get_canonical_geneMer() for gmer in actual_geneMers]
        actual_geneStrings = [[g.get_name() for g in gmer] for gmer in actual_geneObjects]
        actual_number_of_geneMers = len(actual_geneMers)
        # assertion
        expected_number_of_geneMers = len(annotatedGenes) - (4 - 1)
        self.assertEqual(actual_number_of_geneMers, expected_number_of_geneMers)
        self.assertTrue(
            all(any(test_gene[1:] in g for g in actual_geneStrings) for test_gene in annotatedGenes)
        )

    def test_gene_geneMers_size_number_of_genes(self):
        # setup
        annotatedGenes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5"]
        read = Read("read1", annotatedGenes)
        # execution
        actual_geneMers = [x for x in read.get_geneMers(5)]
        actual_geneObjects = [gmer.get_canonical_geneMer() for gmer in actual_geneMers]
        actual_geneStrings = [[g.get_name() for g in gmer] for gmer in actual_geneObjects]
        actual_number_of_geneMers = len(actual_geneMers)
        # assertion
        expected_number_of_geneMers = 1
        self.assertEqual(actual_number_of_geneMers, expected_number_of_geneMers)
        self.assertTrue(
            all(any(test_gene[1:] in g for g in actual_geneStrings) for test_gene in annotatedGenes)
        )

    def test_gene_geneMers_size_more_than_number_of_genes(self):
        # setup
        annotatedGenes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5"]
        read = Read("read1", annotatedGenes)
        # execution
        actual_geneMers = [x for x in read.get_geneMers(6)]
        actual_number_of_geneMers = len(actual_geneMers)
        # assertion
        expected_geneMers = []
        expected_number_of_geneMers = 0
        self.assertEqual(actual_number_of_geneMers, expected_number_of_geneMers)
        self.assertEqual(actual_geneMers, expected_geneMers)

    def test_gene_geneMers_reverse_forward(self):
        # setup
        read1 = Read("read1", ["+gene1", "-gene2", "+gene3"])
        read2 = Read("read2", ["+gene3", "-gene2", "+gene1"])
        # execution
        actual_geneMers1 = [x for x in read1.get_geneMers(3)]
        actual_geneObjects1 = [gmer.get_canonical_geneMer() for gmer in actual_geneMers1]
        actual_geneStrings1 = [[g.get_name() for g in gmer] for gmer in actual_geneObjects1]
        actual_geneMers2 = [x for x in read2.get_geneMers(3)]
        actual_geneObjects2 = [gmer.get_canonical_geneMer() for gmer in actual_geneMers2]
        actual_geneStrings2 = [[g.get_name() for g in gmer] for gmer in actual_geneObjects2]
        # assertion
        self.assertNotEqual(actual_geneMers1, actual_geneMers2)
        self.assertNotEqual(actual_geneObjects1, actual_geneObjects2)
        self.assertNotEqual(actual_geneStrings1, actual_geneStrings2)
