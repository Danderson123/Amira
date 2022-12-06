from gene import Gene

def define_rc_geneMer(geneMer):
    """ return a reversed list of reverse complement geneMers for the input list of geneMers """
    # ensure that all input genes are gene objects
    assert all(type(gene) == Gene for gene in geneMer)
    # reverse the list of gene objects
    reverse_geneMer = list(reversed(geneMer))
    # get the reverse complement gene object for each gene in the geneMer
    rcGeneMer = [gene.reverse_gene() for gene in reverse_geneMer]
    return rcGeneMer

def sort_geneMers(geneMer,
                rcGeneMer):
    geneMerHashes = [g.__hash__() for g in geneMer]
    rcGeneMerHashes = [rc_g.__hash__() for rc_g in rcGeneMer]
    assert not (geneMerHashes == rcGeneMerHashes), "Gene-mer and reverse complement gene-mer are identical"
    sortedGeneMerhashes = list(sorted([geneMerHashes, rcGeneMerHashes]))
    return geneMerHashes, rcGeneMerHashes, sortedGeneMerhashes

def choose_canonical_geneMer(geneMer,
                            geneMerHashes,
                            rcGeneMer,
                            rcGeneMerHashes,
                            sortedGeneMerhashes):
    if sortedGeneMerhashes[0] == geneMerHashes and sortedGeneMerhashes[1] == rcGeneMerHashes:
        return geneMer, rcGeneMer
    if sortedGeneMerhashes[0] == rcGeneMerHashes and sortedGeneMerhashes[1] == geneMerHashes:
        return rcGeneMer, geneMer

def define_geneMer(geneMer):
    # ensure the gene-mer is a list
    assert type(geneMer) == list, "Gene-mer is not a list of Gene objects"
    # ensure the gene-mer is not an empty list
    assert not geneMer == [], "Gene-mer is empty"
    rcGeneMer = define_rc_geneMer(geneMer)
    geneMerHashes, rcGeneMerHashes, sortedGeneMerhashes = sort_geneMers(geneMer,
                                                                        rcGeneMer)
    canonicalGeneMer, reverseCanonicalGeneMer = choose_canonical_geneMer(geneMer,
                                                                        geneMerHashes,
                                                                        rcGeneMer,
                                                                        rcGeneMerHashes,
                                                                        sortedGeneMerhashes)
    return canonicalGeneMer, reverseCanonicalGeneMer

class GeneMer:
    def __init__(self,
                geneMer):
        self.canonicalGeneMer, self.rcGeneMer = define_geneMer(geneMer)
        self.geneMerSize = len(self.canonicalGeneMer)

    def get_canonical_geneMer(self): # returns the canonical gene mer
        return self.canonicalGeneMer
    def get_rc_geneMer(self): # returns the reverse-complemented gene mer
        return self.rcGeneMer
    def get_geneMer_size(self) -> int: # returns the size of this gene mer
        return self.geneMerSize
    def assign_geneMerId(self,
                        geneMerId):
        self.geneMerId = geneMerId
    def __eq__ (self,
            otherGeneMer): # equality comparison
        geneMerGenes = [g.__hash__() for g in self.canonicalGeneMer]
        rcGeneMerGenes = [g.__hash__() for g in self.rcGeneMer]
        otherGeneMerGenes = [o.__hash__() for o in otherGeneMer.get_canonical_geneMer()]
        rcOtherGeneMerGenes = [o.__hash__() for o in otherGeneMer.get_rc_geneMer()]
        return (geneMerGenes == otherGeneMerGenes and rcGeneMerGenes == rcOtherGeneMerGenes)