import unittest

from amira_prototype.construct_gene import convert_string_strand_to_int, convert_int_strand_to_string, reverse_strand, Gene

class TestGeneConstructor(unittest.TestCase):

    def test___init___construction_positive_strand(self):
        # setup
        gene = Gene("+gene1")
        # execution
        actual_name = gene.get_name()
        actual_strand = gene.get_strand()
        # assertion
        expected_name = "gene1"
        expected_strand = 1
        self.assertEqual(actual_name, expected_name)
        self.assertEqual(actual_strand, expected_strand)

    def test___init___construction_negative_strand(self):
        # setup
        gene = Gene("-gene2")
        # execution
        actual_name = gene.get_name()
        actual_strand = gene.get_strand()
        # assertion
        expected_name = "gene2"
        expected_strand = -1
        self.assertEqual(actual_name, expected_name)
        self.assertEqual(actual_strand, expected_strand)

    def test___split_gene_and_positive_strand(self):
        # setup
        gene = Gene("+gene1")
        # execution
        actual_name = gene.get_name()
        actual_strand = gene.get_strand()
        # assertion
        expected_name = "gene1"
        expected_strand = 1
        self.assertEqual(actual_name, expected_name)
        self.assertEqual(actual_strand, expected_strand)

    def test___split_gene_and_negative_strand(self):
        #setup
        gene = Gene("-gene2")
        # execution
        actual_name = gene.get_name()
        actual_strand = gene.get_strand()
        # assertion
        expected_name = "gene2"
        expected_strand = -1
        self.assertEqual(actual_name, expected_name)
        self.assertEqual(actual_strand, expected_strand)

    def test__split_gene_and_positive_strand_positive_strand_char_in_name(self):
        # setup
        gene = Gene("+gene+1")
        # execution
        actual_name = gene.get_name()
        actual_strand = gene.get_strand()
        # assertion
        expected_name = "gene+1"
        expected_strand = 1
        self.assertEqual(actual_name, expected_name)
        self.assertEqual(actual_strand, expected_strand)

    def test__split_gene_and_negative_strand_positive_strand_char_in_name(self):
        # setup
        gene = Gene("-gene+1")
        # execution
        actual_name = gene.get_name()
        actual_strand = gene.get_strand()
        # assertion
        expected_name = "gene+1"
        expected_strand = -1
        self.assertEqual(actual_name, expected_name)
        self.assertEqual(actual_strand, expected_strand)

    def test__split_gene_and_positive_strand_negative_strand_char_in_name(self):
        # setup
        gene = Gene("+gene-2")
        # execution
        actual_name = gene.get_name()
        actual_strand = gene.get_strand()
        # assertion
        expected_name = "gene-2"
        expected_strand = 1
        self.assertEqual(actual_name, expected_name)
        self.assertEqual(actual_strand, expected_strand)

    def test__split_gene_and_negative_strand_negative_strand_char_in_name(self):
        # setup
        gene = Gene("-gene-2")
        # execution
        actual_name = gene.get_name()
        actual_strand = gene.get_strand()
        # assertion
        expected_name = "gene-2"
        expected_strand = -1
        self.assertEqual(actual_name, expected_name)
        self.assertEqual(actual_strand, expected_strand)

    def test__split_gene_and_strand_missing(self):
        self.assertRaises(AssertionError, Gene, "gene1")

    def test__split_gene_missing_and_strand(self):
        self.assertRaises(AssertionError, Gene, "+")
        self.assertRaises(AssertionError, Gene, "-")

    def test__split_gene_missing_and_strand_missing(self):
        self.assertRaises(AssertionError, Gene, "")

    def test___reverse_positive_strand(self):
        # execution
        actual_reverse_strand = reverse_strand(+1)
        # assertion
        expected_reverse_strand = -1
        self.assertEqual(actual_reverse_strand, expected_reverse_strand)

    def test___reverse_negative_strand(self):
        # execution
        actual_reverse_strand = reverse_strand(-1)
        # assertion
        expected_reverse_strand = 1
        self.assertEqual(actual_reverse_strand, expected_reverse_strand)

    def test___reverse_zero_strand(self):
        self.assertRaises(AssertionError, reverse_strand, 0)

    def test___reverse_positive_string_strand(self):
        self.assertRaises(AssertionError, reverse_strand, "+")

    def test___reverse_negative_string_strand(self):
        self.assertRaises(AssertionError, reverse_strand, "-")

    def test___convert_positive_int_strand_to_string(self):
        # execution
        actual_string_strand = convert_int_strand_to_string(1)
        # assertion
        expected_string_strand = "+"
        self.assertEqual(actual_string_strand, expected_string_strand)

    def test___convert_negative_int_strand_to_string(self):
        # execution
        actual_string_strand = convert_int_strand_to_string(-1)
        # assertion
        expected_string_strand = "-"
        self.assertEqual(actual_string_strand, expected_string_strand)

    def test___convert_zero_int_strand_to_string(self):
        self.assertRaises(AssertionError, convert_int_strand_to_string, 0)

    def test___convert_positive_string_strand_to_string(self):
        self.assertRaises(AssertionError, convert_int_strand_to_string, "+")

    def test___convert_positive_string_strand_to_string(self):
        self.assertRaises(AssertionError, convert_int_strand_to_string, "-")

    def test___convert_positive_string_strand_to_int(self):
        # execution
        actual_int_strand = convert_string_strand_to_int("+")
        # assertion
        expected_string_strand = 1
        self.assertEqual(actual_int_strand, expected_string_strand)

    def test___convert_negative_string_strand_to_int(self):
        # execution
        actual_int_strand = convert_string_strand_to_int("-")
        # assertion
        expected_string_strand = -1
        self.assertEqual(actual_int_strand, expected_string_strand)

    def test___convert_positive_int_strand_to_int(self):
        self.assertRaises(AssertionError, convert_string_strand_to_int, 1)

    def test___convert_negative_int_strand_to_int(self):
        self.assertRaises(AssertionError, convert_string_strand_to_int, -1)

    def test___convert_empty_strand_to_int(self):
        self.assertRaises(AssertionError, convert_string_strand_to_int, "")

    def test___reverse_positive_gene(self):
        # setup
        gene = Gene("+gene1")
        # execution
        reverse_gene = gene.reverse_gene()
        actual_name = reverse_gene.get_name()
        actual_strand = reverse_gene.get_strand()
        # assertion
        expected_name = "gene1"
        expected_strand = -1
        self.assertEqual(actual_name, expected_name)
        self.assertEqual(actual_strand, expected_strand)

    def test___reverse_negative_gene(self):
        # setup
        gene = Gene("-gene2")
        # execution
        reverse_gene = gene.reverse_gene()
        actual_name = reverse_gene.get_name()
        actual_strand = reverse_gene.get_strand()
        # assertion
        expected_name = "gene2"
        expected_strand = 1
        self.assertEqual(actual_name, expected_name)
        self.assertEqual(actual_strand, expected_strand)

    def test___eq_positive_gene(self):
        # setup
        gene1 = Gene("+gene1")
        gene2 = Gene("+gene1")
        # execution
        actual_gene1_gene2_comparison = gene1.__eq__(gene2)
        actual_gene2_gene1_comparison = gene2.__eq__(gene1)
        # assertion
        expected_gene1_gene2_comparison = True
        expected_gene2_gene1_comparison = True
        self.assertEqual(actual_gene1_gene2_comparison, expected_gene1_gene2_comparison)
        self.assertEqual(actual_gene2_gene1_comparison, expected_gene2_gene1_comparison)

    def test___eq_negative_gene(self):
        # setup
        gene1 = Gene("-gene2")
        gene2 = Gene("-gene2")
        # execution
        actual_gene1_gene2_comparison = gene1.__eq__(gene2)
        actual_gene2_gene1_comparison = gene2.__eq__(gene1)
        # assertion
        expected_gene1_gene2_comparison = True
        expected_gene2_gene1_comparison = True
        self.assertEqual(actual_gene1_gene2_comparison, expected_gene1_gene2_comparison)
        self.assertEqual(actual_gene2_gene1_comparison, expected_gene2_gene1_comparison)

    def test___eq_positive_negative_same_gene(self):
        # setup
        gene1 = Gene("+gene1")
        gene2 = Gene("-gene1")
        # execution
        actual_gene1_gene2_comparison = gene1.__eq__(gene2)
        actual_gene2_gene1_comparison = gene2.__eq__(gene1)
        # assertion
        expected_gene1_gene2_comparison = True
        expected_gene2_gene1_comparison = True
        self.assertEqual(actual_gene1_gene2_comparison, expected_gene1_gene2_comparison)
        self.assertEqual(actual_gene2_gene1_comparison, expected_gene2_gene1_comparison)

    def test___negative_positive_same_gene(self):
        # setup
        gene1 = Gene("-gene1")
        gene2 = Gene("+gene1")
        # execution
        actual_gene1_gene2_comparison = gene1.__eq__(gene2)
        actual_gene2_gene1_comparison = gene2.__eq__(gene1)
        # assertion
        expected_gene1_gene2_comparison = True
        expected_gene2_gene1_comparison = True
        self.assertEqual(actual_gene1_gene2_comparison, expected_gene1_gene2_comparison)
        self.assertEqual(actual_gene2_gene1_comparison, expected_gene2_gene1_comparison)

    def test___positive_strand_different_gene(self):
        # setup
        gene1 = Gene("+gene1")
        gene2 = Gene("+gene2")
        # execution
        actual_gene1_gene2_comparison = gene1.__eq__(gene2)
        actual_gene2_gene1_comparison = gene2.__eq__(gene1)
        # assertion
        expected_gene1_gene2_comparison = False
        expected_gene2_gene1_comparison = False
        self.assertEqual(actual_gene1_gene2_comparison, expected_gene1_gene2_comparison)
        self.assertEqual(actual_gene2_gene1_comparison, expected_gene2_gene1_comparison)

    def test___negative_strand_different_gene(self):
        # setup
        gene1 = Gene("-gene1")
        gene2 = Gene("-gene2")
        # execution
        actual_gene1_gene2_comparison = gene1.__eq__(gene2)
        actual_gene2_gene1_comparison = gene2.__eq__(gene1)
        # assertion
        expected_gene1_gene2_comparison = False
        expected_gene2_gene1_comparison = False
        self.assertEqual(actual_gene1_gene2_comparison, expected_gene1_gene2_comparison)
        self.assertEqual(actual_gene2_gene1_comparison, expected_gene2_gene1_comparison)

    def test___hash_positive_gene(self):
        # setup
        gene = Gene("+gene1")
        otherGene = Gene("+gene1")
        # execution
        actual_geneHash = gene.__hash__()
        actual_otherGeneHash = otherGene.__hash__()
        # assertion
        self.assertEqual(actual_geneHash, actual_otherGeneHash)

    def test___hash_negative_gene(self):
        # setup
        gene = Gene("-gene2")
        otherGene = Gene("-gene2")
        # execution
        actual_geneHash = gene.__hash__()
        actual_otherGeneHash = otherGene.__hash__()
        # assertion
        self.assertEqual(actual_geneHash, actual_otherGeneHash)

    def test___hash_positive_negative_gene(self):
        # setup
        gene = Gene("+gene1")
        otherGene = Gene("-gene1")
        # execution
        actual_geneHash = gene.__hash__()
        actual_otherGeneHash = otherGene.__hash__()
        # assertion
        self.assertNotEqual(actual_geneHash, actual_otherGeneHash)
