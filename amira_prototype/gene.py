def convert_string_strand_to_int(stringStrand: str) -> int:
    assert stringStrand == "+" or stringStrand == "-"
    if stringStrand == "+":
        intStrand = +1
    else:
        intStrand = -1
    return intStrand

def split_gene_and_strand(gene: str):
    assert gene.replace(" ", "") != "", "Gene information is missing"
    geneStringStrand = gene[0]
    geneName = gene[1:]
    assert geneStringStrand == "-" or geneStringStrand == "+", "Strand information missing for: " + gene
    assert geneName != "", "Gene name information missing for: " + gene
    geneStrand = convert_string_strand_to_int(geneStringStrand)
    return geneName, geneStrand

def reverse_strand(geneStrand: int) -> int:
    assert geneStrand == -1 or geneStrand == 1
    reverseStrandInt = geneStrand * -1
    return reverseStrandInt

def convert_int_strand_to_string(intStrand: int) -> str:
    assert intStrand == -1 or intStrand == 1
    if intStrand == -1:
        stringStrand = "-"
    else:
        stringStrand = "+"
    return stringStrand

class Gene:
    def __init__(self,
                gene: str):
        self.name, self.strand = split_gene_and_strand(gene)
    def get_name(self) -> str:
        return self.name
    def get_strand(self) -> str:
        return self.strand
    def reverse_gene(self):
        reverseStrandInt = reverse_strand(self.strand)
        reverseStrand = convert_int_strand_to_string(reverseStrandInt)
        reverseGene = reverseStrand + self.name
        GeneData = Gene(reverseGene)
        return GeneData
    def __eq__(self,
            otherGene) -> bool:
        otherGeneName, otherGeneStrand = split_gene_and_strand(otherGene)
        return self.get_name() == otherGeneName
    def __hash__(self):
        return hash(self.name)