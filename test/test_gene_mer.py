import sys
sys.path.insert(0, "../amira_prototype")

from gene_mer import define_rc_geneMer, sort_geneMers, choose_canonical_geneMer, define_geneMer, GeneMer
from gene import Gene

def test_define_rc_geneMer():
    # test all inputs correct
    geneMer = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
    rcGeneMer = define_rc_geneMer(geneMer)
    assert rcGeneMer == [Gene("-gene3"), Gene("+gene2"), Gene("-gene1")]
    # test non-Gene object input
    try:
        geneMer = ["+gene1", "-gene2", "+gene3"]
        rcGeneMer = define_rc_geneMer(geneMer)
    except Exception as e:
        assert isinstance(e, AssertionError)

def test_sort_geneMers():
    geneMer = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
    rcGeneMer = [Gene("-gene3"), Gene("+gene2"), Gene("-gene1")]
    geneMerHashes, rcGeneMerHashes, sortedGeneMerhashes = sort_geneMers(geneMer,
                                                                        rcGeneMer)
    assert not geneMerHashes == []
    assert not rcGeneMerHashes == []
    assert (geneMerHashes == sortedGeneMerhashes[0] and rcGeneMerHashes == sortedGeneMerhashes[1]) or (geneMerHashes == sortedGeneMerhashes[1] and rcGeneMerHashes == sortedGeneMerhashes[0])

def test_choose_canonical_geneMer():
    # test k = 3
    geneMer = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
    rcGeneMer = [Gene("-gene3"), Gene("+gene2"), Gene("-gene1")]
    geneMerHashes, rcGeneMerHashes, sortedGeneMerhashes = sort_geneMers(geneMer,
                                                                        rcGeneMer)
    canonicalGeneMer, reverseCanonicalGeneMer = choose_canonical_geneMer(geneMer,
                                                                        geneMerHashes,
                                                                        rcGeneMer,
                                                                        rcGeneMerHashes,
                                                                        sortedGeneMerhashes)
    assert not (canonicalGeneMer == [] or reverseCanonicalGeneMer == [])
    assert not canonicalGeneMer == reverseCanonicalGeneMer, "Canonical and reverse compliment gene-mers are identical"
    assert (canonicalGeneMer == geneMer and reverseCanonicalGeneMer == rcGeneMer) or (canonicalGeneMer == rcGeneMer and reverseCanonicalGeneMer == geneMer)
    # tets k = 1
    geneMer = [Gene("+gene1")]
    rcGeneMer = [Gene("-gene1")]
    geneMerHashes, rcGeneMerHashes, sortedGeneMerhashes = sort_geneMers(geneMer,
                                                                        rcGeneMer)
    canonicalGeneMer, reverseCanonicalGeneMer = choose_canonical_geneMer(geneMer,
                                                                        geneMerHashes,
                                                                        rcGeneMer,
                                                                        rcGeneMerHashes,
                                                                        sortedGeneMerhashes)
    assert not (canonicalGeneMer == [] or reverseCanonicalGeneMer == [])
    assert not [g.__hash__() for g in canonicalGeneMer] == [g.__hash__() for g in reverseCanonicalGeneMer], "Canonical and reverse compliment gene-mers are identical"
    assert (canonicalGeneMer == geneMer and reverseCanonicalGeneMer == rcGeneMer) or (canonicalGeneMer == rcGeneMer and reverseCanonicalGeneMer == geneMer)

def test_define_geneMer():
    # test k = 3
    geneMer = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
    define_geneMer(geneMer)
    # test k = 1
    geneMer = [Gene("+gene1")]
    define_geneMer(geneMer)
    geneMer = [Gene("-gene2")]
    define_geneMer(geneMer)
    # test empty input
    try:
        geneMer = []
        define_geneMer(geneMer)
    except Exception as e:
        assert isinstance(e, AssertionError)
    # test empty string
    try:
        geneMer = [""]
        define_geneMer(geneMer)
    except Exception as e:
        assert isinstance(e, AssertionError)
    try:
        geneMer = [" "]
        define_geneMer(geneMer)
    except Exception as e:
        assert isinstance(e, AssertionError)

def test_eq_geneMer():
    # test k = 3
    genes = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
    otherGenes = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
    geneMer = GeneMer(genes)
    otherGeneMer = GeneMer(otherGenes)
    assert geneMer.__eq__(otherGeneMer) and otherGeneMer.__eq__(geneMer)
    # test k = 3 input gene mer is rc
    genes = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
    otherGenes = [Gene("-gene3"), Gene("+gene2"), Gene("-gene1")]
    geneMer = GeneMer(genes)
    otherGeneMer = GeneMer(otherGenes)
    assert geneMer.__eq__(otherGeneMer) and otherGeneMer.__eq__(geneMer)
    # test k = 3 and gene mers different
    genes = [Gene("+gene1"), Gene("-gene2"), Gene("+gene3")]
    otherGenes = [Gene("-gene2"), Gene("+gene4"), Gene("-gene5")]
    geneMer = GeneMer(genes)
    otherGeneMer = GeneMer(otherGenes)
    assert not (geneMer.__eq__(otherGeneMer) and otherGeneMer.__eq__(geneMer))
    # test k = 1
    genes = [Gene("+gene1")]
    otherGenes = [Gene("+gene1")]
    geneMer = GeneMer(genes)
    otherGeneMer = GeneMer(otherGenes)
    assert geneMer.__eq__(otherGeneMer) and otherGeneMer.__eq__(geneMer)
    # test k = 1 and input gene mer is rc
    genes = [Gene("+gene1")]
    otherGenes = [Gene("-gene1")]
    geneMer = GeneMer(genes)
    otherGeneMer = GeneMer(otherGenes)
    assert geneMer.__eq__(otherGeneMer) and otherGeneMer.__eq__(geneMer)

sys.stderr.write("Testing: define_rc_geneMer\n")
test_define_rc_geneMer()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing: sort_geneMers\n")
test_sort_geneMers()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing: choose_canonical_geneMer\n")
test_choose_canonical_geneMer()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing: define_geneMer\n")
test_define_geneMer()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing: GeneMer.__eq__\n")
test_eq_geneMer()
sys.stderr.write("Test passed\n")