from amira.construct_gene import Gene, hashlib_hash


def define_rc_geneMer(geneMer):
    """return a reversed list of reverse complement geneMers for the input list of geneMers"""
    # ensure that all input genes are gene objects
    assert all(isinstance(gene, Gene) for gene in geneMer)
    # reverse the list of gene objects
    reverse_geneMer = list(reversed(geneMer))
    # get the reverse complement gene object for each gene in the geneMer
    rcGeneMer = [gene.reverse_gene() for gene in reverse_geneMer]
    return rcGeneMer


def sort_geneMers(geneMer, rcGeneMer):
    """return a list of gene-mer and rc gene-mer hashes and a sorted list of the two hash lists"""
    # get the list of hashes for the genes in this gene mer
    geneMerHashes = [g.__hash__() for g in geneMer]
    # get the list of hashes for the genes in the reverse complement gene mer
    rcGeneMerHashes = [rc_g.__hash__() for rc_g in rcGeneMer]
    # ensure that the hashes of the gene mer and rc gene mer are not the same
    # this may require us to specify k must be even
    assert not (
        geneMerHashes == rcGeneMerHashes
    ), "Gene-mer and reverse complement gene-mer are identical"
    # sort the list of hash lists to ensure the canonical gene mer is consistent
    sortedGeneMerhashes = list(sorted([geneMerHashes, rcGeneMerHashes]))
    return geneMerHashes, rcGeneMerHashes, sortedGeneMerhashes


def choose_canonical_geneMer(
    geneMer, geneMerHashes, rcGeneMer, rcGeneMerHashes, sortedGeneMerhashes
):
    """returns the selected canonical gene mer and reverse complement as lists of Gene objects"""
    # select the gene mer at index 0 in the sorted list as the canonical gene mer
    if sortedGeneMerhashes[0] == geneMerHashes and sortedGeneMerhashes[1] == rcGeneMerHashes:
        return geneMer, rcGeneMer
    if sortedGeneMerhashes[0] == rcGeneMerHashes and sortedGeneMerhashes[1] == geneMerHashes:
        return rcGeneMer, geneMer


def define_geneMer(geneMer):
    """returns the canonical gene mer and the reverse complement gene mer for this gene mer"""
    # ensure the gene-mer is a list
    assert isinstance(geneMer, list), "Gene-mer is not a list of Gene objects"
    # ensure the gene-mer is not an empty list
    assert not geneMer == [], "Gene-mer is empty"
    # get the reverse complement gene mer of the seen gene mer
    rcGeneMer = define_rc_geneMer(geneMer)
    # sort the gene mer and rc gene mer so the canonical gene mer is consistent
    geneMerHashes, rcGeneMerHashes, sortedGeneMerhashes = sort_geneMers(geneMer, rcGeneMer)
    # choose the gene mer at index 0 in the sort as the canonical gene mer
    canonicalGeneMer, reverseCanonicalGeneMer = choose_canonical_geneMer(
        geneMer, geneMerHashes, rcGeneMer, rcGeneMerHashes, sortedGeneMerhashes
    )
    return canonicalGeneMer, reverseCanonicalGeneMer


class GeneMer:
    def __init__(self, geneMer: list):
        self.canonicalGeneMer, self.rcGeneMer = define_geneMer(geneMer)
        self.geneMerSize = len(self.canonicalGeneMer)

        def define_direction(canonicalGeneMer, geneMer):
            if not canonicalGeneMer == geneMer:
                return -1
            else:
                return 1

        self.geneMerDirection = define_direction(self.canonicalGeneMer, geneMer)

    def get_canonical_geneMer(self):
        """return the canonical gene mer (list of Gene objects) for this gene mer"""
        return self.canonicalGeneMer

    def get_rc_geneMer(self):
        """return the rc gene mer (list of Gene objects) for this gene mer"""
        return self.rcGeneMer

    def get_geneMerDirection(self):
        return self.geneMerDirection

    def get_geneMer_size(self) -> int:
        """return an integer for the number of gene objects in this gene mer"""
        return self.geneMerSize

    def __eq__(self, otherGeneMer):
        """check if this gene mer is identical to another gene-mer"""
        return (
            self.canonicalGeneMer == otherGeneMer.get_canonical_geneMer()
            and self.rcGeneMer == otherGeneMer.get_rc_geneMer()
        )

    def __hash__(self):
        """return a hash of the gene mer to see if two gene mers are the same"""
        geneMerGenes = tuple([g.__hash__() for g in self.get_canonical_geneMer()])
        return hashlib_hash(geneMerGenes)
