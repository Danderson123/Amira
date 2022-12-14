import unittest
import sys
sys.path.insert(0, "amira_prototype")

from construct_gene_mer import define_rc_geneMer, sort_geneMers, choose_canonical_geneMer, define_geneMer, GeneMer
from construct_gene import Gene


class TestGeneMerConstructor(unittest.TestCase):

    def test___init_genemer(self):
        # setup
        genes = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
        # execution
        gene_mer = GeneMer(genes)
        actual_canonical_gene_mer = gene_mer.get_canonical_geneMer()
        actual_reverse_gene_mer = gene_mer.get_rc_geneMer()
        actual_gene_mer_size = gene_mer.get_geneMer_size()
        actual_gene_mer_direction = gene_mer.get_geneMerDirection()
        # assertion
        if actual_gene_mer_direction == 1:
            expected_gene_mer_direction = 1
            expected_canonical_gene_mer = genes
            expected_reverse_gene_mer = define_rc_geneMer(genes)
            expected_gene_mer_size = 3
        else:
            expected_gene_mer_direction = -1
            expected_reverse_gene_mer = genes
            expected_canonical_gene_mer = define_rc_geneMer(genes)
            expected_gene_mer_size = 3
        self.assertEqual(actual_canonical_gene_mer, expected_canonical_gene_mer)
        self.assertEqual(actual_reverse_gene_mer, expected_reverse_gene_mer)
        self.assertEqual(actual_gene_mer_direction, expected_gene_mer_direction)
        self.assertEqual(actual_gene_mer_size, expected_gene_mer_size)

    def test___define_rc_geneMer(self):
        # setup
        geneMer = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
        # execution
        actual_rcGeneMer = define_rc_geneMer(geneMer)
        # assertion
        expected_rcGeneMer = [Gene("-gene3"), Gene("+gene2"), Gene("-gene1")]
        self.assertEqual(actual_rcGeneMer, expected_rcGeneMer)

    def test__define_rc_geneMer_non_gene_input(self):
        # setup
        geneMer = ["+gene1", "-gene2", "+gene3"]
        # assertion
        self.assertRaises(AssertionError, define_rc_geneMer, geneMer)

    def test__define_rc_geneMer_empty_input(self):
        # setup
        geneMer = []
        # execution
        actual_rcGeneMer = define_rc_geneMer(geneMer)
        # assertion
        expected_geneMer = []
        self.assertEqual(actual_rcGeneMer, expected_geneMer)

    def test___sort_geneMers(self):
        # setup
        geneMer = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
        rcGeneMer = [Gene("-gene3"), Gene("+gene2"), Gene("-gene1")]
        # execution
        actual_geneMerHashes, actual_rcGeneMerHashes, actual_sortedGeneMerhashes = sort_geneMers(geneMer,
                                                                                            rcGeneMer)
        # assertion
        self.assertTrue((actual_geneMerHashes == actual_sortedGeneMerhashes[0] and actual_rcGeneMerHashes == actual_sortedGeneMerhashes[1]) or (actual_geneMerHashes == actual_sortedGeneMerhashes[1] and actual_rcGeneMerHashes == actual_sortedGeneMerhashes[0]))

    def test___choose_canonical_geneMer_size_3(self):
        # setup
        geneMer = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
        rcGeneMer = [Gene("-gene3"), Gene("+gene2"), Gene("-gene1")]
        geneMerHashes, rcGeneMerHashes, sortedGeneMerhashes = sort_geneMers(geneMer,
                                                                            rcGeneMer)
        # execution
        actual_canonicalGeneMer, actual_reverseCanonicalGeneMer = choose_canonical_geneMer(geneMer,
                                                                                        geneMerHashes,
                                                                                        rcGeneMer,
                                                                                        rcGeneMerHashes,
                                                                                        sortedGeneMerhashes)
        # assertion
        self.assertNotEqual(actual_canonicalGeneMer, [])
        self.assertNotEqual(actual_reverseCanonicalGeneMer, [])
        self.assertNotEqual(actual_canonicalGeneMer, actual_reverseCanonicalGeneMer)
        self.assertTrue((actual_canonicalGeneMer == geneMer and actual_reverseCanonicalGeneMer == rcGeneMer) or (actual_canonicalGeneMer == rcGeneMer and actual_reverseCanonicalGeneMer == geneMer))

    def test___choose_canonical_geneMer_size_1(self):
        # setup
        geneMer = [Gene("+gene1")]
        rcGeneMer = [Gene("-gene1")]
        geneMerHashes, rcGeneMerHashes, sortedGeneMerhashes = sort_geneMers(geneMer,
                                                                            rcGeneMer)
        # execution
        actual_canonicalGeneMer, actual_reverseCanonicalGeneMer = choose_canonical_geneMer(geneMer,
                                                                                        geneMerHashes,
                                                                                        rcGeneMer,
                                                                                        rcGeneMerHashes,
                                                                                        sortedGeneMerhashes)
        # assertion
        self.assertNotEqual(actual_canonicalGeneMer, [])
        self.assertNotEqual(actual_reverseCanonicalGeneMer, [])
        self.assertNotEqual([g.__hash__() for g in actual_canonicalGeneMer], [g.__hash__() for g in actual_reverseCanonicalGeneMer])
        self.assertTrue((actual_canonicalGeneMer == geneMer and actual_reverseCanonicalGeneMer == rcGeneMer) or (actual_canonicalGeneMer == rcGeneMer and actual_reverseCanonicalGeneMer == geneMer))

    def test___define_geneMer_size_3(self):
        # setup
        geneMer = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
        # execution
        actual_canonicalGeneMer, actual_reverseCanonicalGeneMer = define_geneMer(geneMer)
        # assertion
        expected_geneMer = geneMer
        expected_rcGeneMer = define_rc_geneMer(geneMer)
        self.assertTrue((actual_canonicalGeneMer == expected_geneMer and actual_reverseCanonicalGeneMer == expected_rcGeneMer) or (actual_canonicalGeneMer == expected_rcGeneMer and actual_reverseCanonicalGeneMer == expected_geneMer))

    def test___define_geneMer_size_1_positive(self):
        # setup
        geneMer = [Gene("+gene1")]
        # execution
        actual_canonicalGeneMer, actual_reverseCanonicalGeneMer = define_geneMer(geneMer)
        # assertion
        expected_geneMer = geneMer
        expected_rcGeneMer = define_rc_geneMer(geneMer)
        self.assertTrue((actual_canonicalGeneMer == expected_geneMer and actual_reverseCanonicalGeneMer == expected_rcGeneMer) or (actual_canonicalGeneMer == expected_rcGeneMer and actual_reverseCanonicalGeneMer == expected_geneMer))

    def test___define_geneMer_size_1_negative(self):
        # setup
        geneMer = [Gene("-gene2")]
        # execution
        actual_canonicalGeneMer, actual_reverseCanonicalGeneMer = define_geneMer(geneMer)
        # assertion
        expected_geneMer = geneMer
        expected_rcGeneMer = define_rc_geneMer(geneMer)
        self.assertTrue((actual_canonicalGeneMer == expected_geneMer and actual_reverseCanonicalGeneMer == expected_rcGeneMer) or (actual_canonicalGeneMer == expected_rcGeneMer and actual_reverseCanonicalGeneMer == expected_geneMer))

    def test__define_geneMer_empty_list(self):
        # setup
        geneMer = []
        # assertion
        self.assertRaises(AssertionError, define_geneMer, geneMer)

    def test__define_geneMer_empty_string(self):
        # setup
        geneMer = [""]
        # assertion
        self.assertRaises(AssertionError, define_geneMer, geneMer)

    def test__define_geneMer_space_string(self):
        # setup
        geneMer = [" "]
        # assertion
        self.assertRaises(AssertionError, define_geneMer, geneMer)

    def test___eq_geneMer_equal_size_3(self):
        # setup
        genes = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
        otherGenes = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
        # execution
        geneMer = GeneMer(genes)
        otherGeneMer = GeneMer(otherGenes)
        # assertion
        self.assertTrue((geneMer.__eq__(otherGeneMer) and otherGeneMer.__eq__(geneMer)))

    def test___eq_rc_geneMer_equal_size_3(self):
        # setup
        genes = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
        otherGenes = [Gene("-gene3"), Gene("+gene2"), Gene("-gene1")]
        # execution
        geneMer = GeneMer(genes)
        otherGeneMer = GeneMer(otherGenes)
        # assertion
        self.assertTrue((geneMer.__eq__(otherGeneMer) and otherGeneMer.__eq__(geneMer)))

    def test___eq_geneMer_different_size_3(self):
        # setup
        genes = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
        otherGenes = [Gene("-gene2"), Gene("+gene4"), Gene("-gene5")]
        # execution
        geneMer = GeneMer(genes)
        otherGeneMer = GeneMer(otherGenes)
        # assertion
        self.assertFalse((geneMer.__eq__(otherGeneMer) and otherGeneMer.__eq__(geneMer)))

    def test___eq_geneMer_equal_size_1(self):
        # setup
        genes = [Gene("+gene1")]
        otherGenes = [Gene("+gene1")]
        # execution
        geneMer = GeneMer(genes)
        otherGeneMer = GeneMer(otherGenes)
        # assertion
        self.assertTrue((geneMer.__eq__(otherGeneMer) and otherGeneMer.__eq__(geneMer)))

    def test___eq_rc_geneMer_equal_size_1(self):
        # setup
        genes = [Gene("+gene1")]
        otherGenes = [Gene("-gene1")]
        # execution
        geneMer = GeneMer(genes)
        otherGeneMer = GeneMer(otherGenes)
        # assertion
        self.assertTrue((geneMer.__eq__(otherGeneMer) and otherGeneMer.__eq__(geneMer)))

    def test___eq_geneMer_different_size_1(self):
        # setup
        genes = [Gene("+gene1")]
        otherGenes = [Gene("-gene2")]
        # execution
        geneMer = GeneMer(genes)
        otherGeneMer = GeneMer(otherGenes)
        # assertion
        self.assertFalse((geneMer.__eq__(otherGeneMer) and otherGeneMer.__eq__(geneMer)))

    def test___hash_geneMer_equal_size_3(self):
        # setup
        genes = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
        otherGenes = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
        # execution
        geneMer = GeneMer(genes)
        otherGeneMer = GeneMer(otherGenes)
        # assertion
        self.assertEqual(geneMer.__hash__(), otherGeneMer.__hash__())

    def test___hash_geneMer_rc_equal_size_3(self):
        # setup
        genes = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
        otherGenes = [Gene("-gene3"), Gene("+gene2"), Gene("-gene1")]
        # execution
        geneMer = GeneMer(genes)
        otherGeneMer = GeneMer(otherGenes)
        # assertion
        self.assertEqual(geneMer.__hash__(), otherGeneMer.__hash__())

    def test___hash_geneMer_different_overlapping_size_3(self):
        # setup
        genes = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
        otherGenes = [Gene("+gene2"), Gene("-gene3"), Gene("+gene4")]
        # execution
        geneMer = GeneMer(genes)
        otherGeneMer = GeneMer(otherGenes)
        # assertion
        self.assertNotEqual(geneMer.__hash__(), otherGeneMer.__hash__())

    def test___hash_geneMer_different_not_overlapping_size_3(self):
        # setup
        genes = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
        otherGenes = [Gene("+gene4"), Gene("-gene5"), Gene("+gene6")]
        # execution
        geneMer = GeneMer(genes)
        otherGeneMer = GeneMer(otherGenes)
        # assertion
        self.assertNotEqual(geneMer.__hash__(), otherGeneMer.__hash__())
