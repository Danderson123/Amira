import unittest
import sys
sys.path.insert(0, "amira_prototype")

from construct_gene import convert_string_strand_to_int, convert_int_strand_to_string, reverse_strand, Gene

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

    def test___eq_and_hash_properties(self):
        """
        This test is mostly to support a discussion in a PR, as I think there might be a bug with Gene.__eq__() and Gene.__hash__().
        It might not make it into the final set of amira tests.
        The properties to satisfy when implementing __eq__() and __hash__() in python classes are:
            If a == b then hash(a) == hash(b)
            If hash(a) == hash(b), then a might equal b
            If hash(a) != hash(b), then a != b
        For more details, see https://eng.lyft.com/hashing-and-equality-in-python-2ea8c738fb9d
        This test checks properties 1 and 3. Property 2 is cumbersome to test, but can be done easily with mocking,
        but we are skipping this here.
        The consequences of this is that the gene might not be correctly used in a python set or dictionary (see
        next tests)
        """
        gene_plus_strand = Gene("+gene")
        gene_minus_strand = Gene("-gene")
        another_gene = Gene("+another_gene")

        # property 3
        # I am testing property 3 first because it succeeds
        self.assertNotEqual(hash(gene_plus_strand), hash(another_gene))
        self.assertNotEqual(gene_plus_strand, another_gene)

        # property one
        self.assertEqual(gene_plus_strand, gene_minus_strand)
        # Fixme: this fails: if two objects are equal, they must hash to a same hash value
        self.assertEqual(hash(gene_plus_strand), hash(gene_minus_strand))


    def test___hash_and_eq_in_a_set(self):
        """
        This test inserts one gene in a set and queries if we find it.
        """
        gene_plus_strand = Gene("+gene")
        gene_minus_strand = Gene("-gene")
        genes = {gene_plus_strand}

        # assert that the two genes are equal
        self.assertEqual(gene_plus_strand, gene_minus_strand)

        # if they are equal, i.e. the same object, and as we added one of the objects to the set, we should be
        # able to find both:
        self.assertIn(gene_plus_strand, genes)  # this works
        self.assertIn(gene_minus_strand, genes)  # Fixme: this fails but should work: is the same gene, just on the minus strand

    def test___hash_and_eq_in_a_dictionary___counting_use_case(self):
        """
        This test counts how many times we've seen a gene
        """
        gene_plus_strand = Gene("+gene")
        gene_minus_strand = Gene("-gene")
        another_gene = Gene("+another_gene")

        list_of_genes = [gene_minus_strand, another_gene, gene_plus_strand, gene_minus_strand, gene_minus_strand,
                         another_gene, another_gene, gene_plus_strand]

        # count the number of genes
        from collections import defaultdict
        genes_to_count = defaultdict(int)
        for gene in list_of_genes:
            genes_to_count[gene] += 1

        # Fixme: we should have 5 counts of "gene" and 3 counts of "another_gene"
        self.assertEqual(3, genes_to_count[another_gene])  # this works
        self.assertEqual(5, genes_to_count[gene_plus_strand])  # Fixme: this fails
