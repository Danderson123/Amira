from construct_gene import Gene
from construct_gene_mer import GeneMer

def convert_genes(annotatedGenes):
    """ return a list of Gene objects for each gene in the list of genes for this read """
    geneMerGenes = [Gene(g) for g in annotatedGenes]
    return geneMerGenes

class Read:
    def __init__(self,
                readId: str,
                annotatedGenes):
        self.readId = readId
        self.numberOfGenes = len(annotatedGenes)
        self.listOfGenes = convert_genes(annotatedGenes)
    def get_readId(self) -> str:
        """ return a string identifier for this read """
        return self.readId
    def get_genes(self) -> list:
        """ return a list of Gene objects for this read """
        return self.listOfGenes
    def get_geneMers(self,
                    kmerSize: int):
        """ return a generator to create GeneMer objects of length kmerSize for this read """
        # if the number of genes on the read is equal to or more than the kmerSize, get the gene-mers
        if self.numberOfGenes > kmerSize - 1:
            # iterate through the list of genes by index
            for i in range(self.numberOfGenes - (kmerSize - 1)):
                # take a slice of the list of Genes from index i to i + kmerSize
                geneMerGenes = self.listOfGenes[i: i + kmerSize]
                # convert the list of Gene objects to a GeneMer object
                geneMer = GeneMer(geneMerGenes)
                # add the geneMer to the list of gene mers for this read
                yield geneMer