from amira.construct_gene import Gene
from amira.construct_gene_mer import GeneMer


def convert_genes(annotatedGenes):
    """return a list of Gene objects for each gene in the list of genes for this read"""
    geneMerGenes = [Gene(g) for g in annotatedGenes]
    return geneMerGenes


class Read:
    def __init__(self, readId: str, annotatedGenes, annotatedGenePositions=None):
        self.readId = readId
        self.numberOfGenes = len(annotatedGenes)
        self.listOfGenes = convert_genes(annotatedGenes)
        self._annotatedGenes = annotatedGenes
        self._annotatedGenePositions = annotatedGenePositions

    def get_readId(self) -> str:
        """return a string identifier for this read"""
        return self.readId

    def get_genes(self) -> list:
        """return a list of Gene objects for this read"""
        return self.listOfGenes

    def get_number_of_genes(self) -> int:
        """return an int of the number of genes annotated for this read"""
        return self.numberOfGenes

    def get_annotatedGenes(self) -> list:
        return self._annotatedGenes

    def get_annotatedGenePositions(self) -> list:
        return self._annotatedGenePositions

    def get_geneMers(self, kmerSize: int):
        """return a generator to create GeneMer objects of length kmerSize for this read"""
        geneMers = []
        geneMerPositions = []
        # check if number of genes is greater than gene-mer size
        if self.get_number_of_genes() > kmerSize - 1:
            # iterate through the list of genes by index
            for i in range(self.get_number_of_genes() - (kmerSize - 1)):
                # take a slice of the list of Genes from index i to i + kmerSize
                geneMerGenes = self.get_genes()[i : i + kmerSize]
                # get the gene-mer position on the read
                if self.get_annotatedGenePositions():
                    genePositions = self.get_annotatedGenePositions()[i : i + kmerSize]
                    gene_mer_start = genePositions[0][0]
                    gene_mer_end = genePositions[-1][1]
                    geneMerPositions.append((gene_mer_start, gene_mer_end))
                else:
                    geneMerPositions.append(None)
                # convert the list of Gene objects to a GeneMer object
                geneMer = GeneMer(geneMerGenes)
                # add the geneMer to the list of gene mers for this read
                geneMers.append(geneMer)
        return geneMers, geneMerPositions
