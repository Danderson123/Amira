import sys
sys.path.insert(0, "../amira_prototype")

from construct_gene import convert_string_strand_to_int, split_gene_and_strand, convert_int_strand_to_string, reverse_strand, Gene

def test_convert_string_strand_to_int():
    # test input correct
    intStrand = convert_string_strand_to_int("+")
    assert intStrand == 1
    intStrand = convert_string_strand_to_int("-")
    assert intStrand == -1
    # test input int
    try:
        intStrand = convert_string_strand_to_int(+1)
    except Exception as e:
        assert isinstance(e, AssertionError)
    try:
        intStrand = convert_string_strand_to_int(-1)
    except Exception as e:
        assert isinstance(e, AssertionError)
    # test input missing
    try:
        intStrand = convert_string_strand_to_int("")
    except Exception as e:
        assert isinstance(e, AssertionError)

def test_split_gene_and_strand():
    # test all inputs correct
    returnName, returnStrand = split_gene_and_strand("+gene1")
    assert returnStrand == 1
    assert returnName == "gene1"
    returnName, returnStrand = split_gene_and_strand("-gene2")
    assert returnStrand == -1
    assert returnName == "gene2"
    # test strand information missing, assert there is an assertionerror
    try:
        split_gene_and_strand("gene1")
    except Exception as e:
        assert isinstance(e, AssertionError)
    # test name information missing
    try:
        returnName, returnStrand = split_gene_and_strand("+")
    except Exception as e:
        assert isinstance(e, AssertionError)
    try:
        returnName, returnStrand = split_gene_and_strand("-")
    except Exception as e:
        assert isinstance(e, AssertionError)
    # test gene information missing
    try:
        returnName, returnStrand = split_gene_and_strand("")
    except Exception as e:
        assert isinstance(e, AssertionError)
    # test behaviour when "+" and "-" characters in gene name
    returnName, returnStrand = split_gene_and_strand("+gene+1")
    assert returnStrand == 1
    assert returnName == "gene+1"
    returnName, returnStrand = split_gene_and_strand("-gene+1")
    assert returnStrand == -1
    assert returnName == "gene+1"
    returnName, returnStrand = split_gene_and_strand("+gene-2")
    assert returnStrand == 1
    assert returnName == "gene-2"
    returnName, returnStrand = split_gene_and_strand("-gene-2")
    assert returnStrand == -1
    assert returnName == "gene-2"

def test_reverse_strand():
    # test correct input
    reverseStrand = reverse_strand(+1)
    assert reverseStrand == -1
    reverseStrand = reverse_strand(-1)
    assert reverseStrand == 1
    # test input wrong
    try:
        reverseStrand = reverse_strand(0)
    except Exception as e:
        assert isinstance(e, AssertionError)

def test_convert_int_strand_to_string():
    # test correct input
    stringStrand = convert_int_strand_to_string(-1)
    assert stringStrand == "-"
    stringStrand = convert_int_strand_to_string(1)
    assert stringStrand == "+"
    # test input wrong
    try:
        stringStrand = convert_int_strand_to_string(0)
    except Exception as e:
        assert isinstance(e, AssertionError)
    # test input string
    try:
        stringStrand = convert_int_strand_to_string("+")
    except Exception as e:
        assert isinstance(e, AssertionError)

def test_reverse_gene():
    # test correct inputs
    geneData = Gene("+gene1")
    reverseGeneData = geneData.reverse_gene()
    assert reverseGeneData.get_name() == "gene1"
    assert reverseGeneData.get_strand() == -1
    geneData = Gene("-gene2")
    reverseGeneData = geneData.reverse_gene()
    assert reverseGeneData.get_name() == "gene2"
    assert reverseGeneData.get_strand() == 1
    # test missing strand
    try:
        geneData = Gene("gene1")
        reverseGeneData = geneData.reverse_gene()
    except Exception as e:
        assert isinstance(e, AssertionError)
    # test missing gene
    try:
        geneData = Gene("+")
        reverseGeneData = geneData.reverse_gene()
    except Exception as e:
        assert isinstance(e, AssertionError)
    # tets empty input
    try:
        geneData = Gene("")
        reverseGeneData = geneData.reverse_gene()
    except Exception as e:
        assert isinstance(e, AssertionError)

def test_eq_gene():
    # test inputs the same
    geneData = Gene("+gene1")
    otherGeneData = Gene("+gene1")
    assert geneData.__eq__(otherGeneData) and otherGeneData.__eq__(geneData)
    geneData = Gene("-gene1")
    otherGeneData = Gene("-gene1")
    assert geneData.__eq__(otherGeneData) and otherGeneData.__eq__(geneData)
    # test genes the same but strands different
    geneData = Gene("+gene1")
    otherGeneData = Gene("-gene1")
    assert geneData.__eq__(otherGeneData) and otherGeneData.__eq__(geneData)
    geneData = Gene("-gene1")
    otherGeneData = Gene("+gene1")
    assert geneData.__eq__(otherGeneData) and otherGeneData.__eq__(geneData)
    # test strands the same but genes different
    geneData = Gene("+gene1")
    otherGeneData = Gene("+gene2")
    assert not (geneData.__eq__(otherGeneData) and otherGeneData.__eq__(geneData))
    geneData = Gene("-gene2")
    otherGeneData = Gene("-gene1")
    assert not (geneData.__eq__(otherGeneData) and otherGeneData.__eq__(geneData))

def test_hash_gene():
    # test inputs the same
    geneHash = Gene("+gene1").__hash__()
    otherGeneHash = Gene("+gene1").__hash__()
    assert geneHash == otherGeneHash
    geneHash = Gene("-gene1").__hash__()
    otherGeneHash = Gene("-gene1").__hash__()
    assert geneHash == otherGeneHash
    # test genes the same but strands different
    geneHash = Gene("+gene1").__hash__()
    otherGeneHash = Gene("-gene1").__hash__()
    assert not geneHash == otherGeneHash
    # test empty input
    try:
        geneHash = Gene("").__hash__()
    except Exception as e:
        assert isinstance(e, AssertionError)

sys.stderr.write("Testing: convert_string_strand_to_int\n")
test_convert_string_strand_to_int()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing: split_gene_and_strand\n")
test_split_gene_and_strand()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing: reverse_strand\n")
test_reverse_strand()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing: convert_int_strand_to_string\n")
test_convert_int_strand_to_string()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing: Gene.reverse_gene\n")
test_reverse_gene()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing: Gene.__eq__\n")
test_eq_gene()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing: Gene.__hash__\n")
test_hash_gene()
sys.stderr.write("Test passed\n")
