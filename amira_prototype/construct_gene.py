import hashlib
import pickle


def hashlib_hash(value):
    # Serialize the value using pickle
    serialized_value = pickle.dumps(value)
    hex_digest = hashlib.sha256(serialized_value).hexdigest()
    integer_hash = int(hex_digest, 16)
    return integer_hash


def convert_string_strand_to_int(stringStrand: str) -> int:
    """convert a string strand to an integer strand and return the value"""
    # ensure the string strand is "+" or "-" only
    assert stringStrand == "+" or stringStrand == "-"
    # convert the string strand to the necessary integer
    if stringStrand == "+":
        intStrand = +1
    else:
        intStrand = -1
    return intStrand


def reverse_strand(geneStrand: int) -> int:
    """returns the reverse sign for an input integer strand"""
    # ensure the integer strand input is only 1 or -1
    assert geneStrand == -1 or geneStrand == 1
    # reverse the input strand by multiplying by -1
    reverseStrandInt = geneStrand * -1
    return reverseStrandInt


def convert_int_strand_to_string(intStrand: int) -> str:
    """returns a string strand consistent with the integer strand inputted"""
    # ensure the integer strand input is only 1 or -1
    assert intStrand == -1 or intStrand == 1
    # if the integer strand is -1 return "-"
    if intStrand == -1:
        stringStrand = "-"
    # if the integer strand is +1 return "+"
    else:
        stringStrand = "+"
    return stringStrand


class Gene:
    def __init__(self, gene: str):
        def split_gene_and_strand(gene: str) -> tuple:
            """returns the separated string name and integer strand for an input gene"""
            # ensure the input gene is not an empty string
            assert gene.replace(" ", "") != "", "Gene information is missing"
            # isolate the strand
            geneStringStrand = gene[0]
            # isolate the gene name and replace any spaces with an underscore
            geneName = gene[1:].replace(" ", "_")
            # ensure the strand is only either "+" or "-"
            assert geneStringStrand == "-" or geneStringStrand == "+", (
                "Strand information missing for: " + gene
            )
            # ensure that the name of the gene is not missing
            assert geneName != "", "Gene name information missing for: " + gene
            # convert the string strand to an integer, either 1 or -1
            geneStrand = convert_string_strand_to_int(geneStringStrand)
            return geneName, geneStrand

        self.name, self.strand = split_gene_and_strand(gene)

    def get_name(self) -> str:
        """return the name of the gene without strand information"""
        return self.name

    def get_strand(self) -> int:
        """return the integer strand of the gene without gene name information"""
        return self.strand

    def reverse_gene(self):
        """return a Gene object that is the reverse of this Gene object"""
        reverseStrandInt = reverse_strand(self.strand)
        reverseStrand = convert_int_strand_to_string(reverseStrandInt)
        reverseGene = reverseStrand + self.name
        GeneData = Gene(reverseGene)
        return GeneData

    def __eq__(self, otherGene) -> bool:
        """return a bool to check if the gene names of two Gene objects are the same"""
        geneStrand = convert_int_strand_to_string(self.get_strand())
        otherGeneStrand = convert_int_strand_to_string(otherGene.get_strand())
        return geneStrand == otherGeneStrand and self.get_name() == otherGene.get_name()

    def __hash__(self) -> int:
        """return a hash of the gene name multiplied by the strand"""
        return hashlib_hash(self.name) * self.strand
