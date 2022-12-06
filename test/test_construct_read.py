import sys
sys.path.insert(0, "../amira_prototype")

from construct_read import Read

def test_get_geneMers():
    annotatedGenes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5"]
    readData = Read("read1",
                    annotatedGenes)
    # test a gene-mer length of 1
    geneMers = readData.get_geneMers(1)
    geneObjects = [gmer.get_canonical_geneMer() for gmer in geneMers]
    geneStrings = [[g.get_name() for g in gmer] for gmer in geneObjects]
    assert len(geneMers) == len(annotatedGenes)
    assert geneStrings == [['gene1'], ['gene2'], ['gene3'], ['gene4'], ['gene5']] or geneStrings == list(reversed([['gene1'], ['gene2'], ['gene3'], ['gene4'], ['gene5']]))
    # test a gene-mer length of 2
    geneMers = readData.get_geneMers(2)
    geneObjects = [gmer.get_canonical_geneMer() for gmer in geneMers]
    geneStrings = [[g.get_name() for g in gmer] for gmer in geneObjects]
    assert len(geneMers) == len(annotatedGenes) - (2 - 1)
    assert all(any(test_gene[1:] in g for g in geneStrings) for test_gene in annotatedGenes)
    # test a gene-mer length of 3
    geneMers = readData.get_geneMers(3)
    geneObjects = [gmer.get_canonical_geneMer() for gmer in geneMers]
    geneStrings = [[g.get_name() for g in gmer] for gmer in geneObjects]
    assert len(geneMers) == len(annotatedGenes) - (3 - 1)
    assert all(any(test_gene[1:] in g for g in geneStrings) for test_gene in annotatedGenes)
    # test a gene-mer length of 4
    geneMers = readData.get_geneMers(4)
    geneObjects = [gmer.get_canonical_geneMer() for gmer in geneMers]
    geneStrings = [[g.get_name() for g in gmer] for gmer in geneObjects]
    assert len(geneMers) == len(annotatedGenes) - (4 - 1)
    assert all(any(test_gene[1:] in g for g in geneStrings) for test_gene in annotatedGenes)
    # test a gene-mer length = len(annotatedGenes)
    geneMers = readData.get_geneMers(5)
    geneObjects = [gmer.get_canonical_geneMer() for gmer in geneMers]
    geneStrings = [[g.get_name() for g in gmer] for gmer in geneObjects]
    assert len(geneMers) == 1
    assert all(any(test_gene[1:] in g for g in geneStrings) for test_gene in annotatedGenes)
    # test a gene-mer length > len(annotatedGenes)
    geneMers = readData.get_geneMers(6)
    assert len(geneMers) == 0
    assert geneMers == []

sys.stderr.write("Testing construct_read: get_geneMers\n")
test_get_geneMers()
sys.stderr.write("Test passed\n")