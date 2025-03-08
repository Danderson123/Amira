import os
import shutil
import unittest

from amira.result_utils import compare_reads_to_references


class TestResultUtils(unittest.TestCase):

    def test___compare_reads_to_references_1(self):
        # setup
        reference_genes = {"sul2NG_0481161": {}}
        with open("tests/test_allele_1/01.reference_alleles.fasta") as f:
            seqs = f.read().split(">")[1:]
        for row in seqs:
            header = row.split("\n")[0]
            seq = "".join(row.split("\n")[1:])
            reference_genes["sul2NG_0481161"][header] = seq
        if not os.path.exists("tests/test_allele_1/out"):
            os.mkdir("tests/test_allele_1/out")
        # execution
        actual_result = compare_reads_to_references(
            "tests/test_allele_1/sul2NG_0481161_1.fastq.gz",
            "tests/test_allele_1/out",
            reference_genes,
            "racon",
            {},
            {},
            True,
            "minimap2",
            0.9,
            0.9,
            "samtools",
        )
        # assertion
        self.assertEqual(actual_result["Determinant name"], "sul2")
        self.assertEqual(actual_result["Closest reference"], "NG_051852")
        self.assertEqual(actual_result["Identity (%)"], 100.0)
        self.assertEqual(actual_result["Coverage (%)"], 100.0)
        self.assertEqual(actual_result["Number of reads"], 90)
        # cleanup
        shutil.rmtree("tests/test_allele_1/out")

    def test___compare_reads_to_references_2(self):
        # setup
        reference_genes = {"catB8aac6IbNG_0520521": {}}
        with open("tests/test_allele_2/01.reference_alleles.fasta") as f:
            seqs = f.read().split(">")[1:]
        for row in seqs:
            header = row.split("\n")[0]
            seq = "".join(row.split("\n")[1:])
            reference_genes["catB8aac6IbNG_0520521"][header] = seq
        if not os.path.exists("tests/test_allele_2/out"):
            os.mkdir("tests/test_allele_2/out")
        # execution
        actual_result = compare_reads_to_references(
            "tests/test_allele_2/catB8aac6IbNG_0520521_1.fastq.gz",
            "tests/test_allele_2/out",
            reference_genes,
            "racon",
            {},
            {},
            True,
            "minimap2",
            0.9,
            0.9,
            "samtools",
        )
        # assertion
        self.assertEqual(actual_result["Determinant name"], "catB3")
        self.assertEqual(actual_result["Closest reference"], "NG_052455")
        self.assertEqual(actual_result["Identity (%)"], 100.0)
        self.assertEqual(actual_result["Coverage (%)"], 69.8)
        self.assertEqual(actual_result["Number of reads"], 106)
        # cleanup
        shutil.rmtree("tests/test_allele_2/out")

    def test___compare_reads_to_references_3(self):
        # setup
        reference_genes = {"sul1NG_0480981": {}}
        with open("tests/test_allele_3/01.reference_alleles.fasta") as f:
            seqs = f.read().split(">")[1:]
        for row in seqs:
            header = row.split("\n")[0]
            seq = "".join(row.split("\n")[1:])
            reference_genes["sul1NG_0480981"][header] = seq
        if not os.path.exists("tests/test_allele_3/out"):
            os.mkdir("tests/test_allele_3/out")
        # execution
        actual_result = compare_reads_to_references(
            "tests/test_allele_3/sul1NG_0480981_1.fastq.gz",
            "tests/test_allele_3/out",
            reference_genes,
            "racon",
            {},
            {},
            True,
            "minimap2",
            0.9,
            0.9,
            "samtools",
        )
        # assertion
        self.assertEqual(actual_result["Determinant name"], "sul1")
        self.assertEqual(actual_result["Closest reference"], "NG_048082")
        self.assertEqual(actual_result["Identity (%)"], 100.0)
        self.assertEqual(actual_result["Coverage (%)"], 100.0)
        self.assertEqual(actual_result["Number of reads"], 141)
        # cleanup
        shutil.rmtree("tests/test_allele_3/out")
