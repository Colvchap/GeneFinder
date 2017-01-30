# -*- coding: utf-8 -*-
"""
Project 1
We are finding genes in nucleotide sequences of DNA

We were given the following list of functions with words describing what how
they should perform. However, when this was given to me, the definitions
contained no code, so everything written before a #TODO is my original
work.

While some of the variables have humorous names, this simple program is
extremely applicable to the real world through modern biology and genome
mapping.

@author: Colvin Chapaman

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return'G'
    elif nucleotide == 'G':
        return 'C'
    # TODO: implement this
    pass


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    i = -1
    AND = ''        # because AND is DNA backwards -"and" is something special
    for letter in dna:
        Complement = get_complement(dna[i])
        AND = AND + Complement
        i += -1
    return AND
    # TODO: implement this
    pass


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    separate_dna = []
    for nucleotide in dna:          # constructs list of nucleotides
        separate_dna.append(nucleotide)
    p = 0
    three_prime_end = ''
    while True:
        if len(dna) - p <= 2:
                return dna
        else:
            codon = separate_dna[p] + separate_dna[p+1] + separate_dna[p+2]

            if codon == "TAG" or codon == "TAA" or codon == "TGA":
                if p >= 3:
                    return three_prime_end
            p += 3
            three_prime_end = three_prime_end + codon
    # TODO: implement this
    pass


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    p = 0                 # Position Index
    ORFs = []             # List of ORFs
    n = 0                 # number of ORFs found
    flag = 'go'           # the loop keeps going until ending procedure

    while flag == 'go':
        if p+3 >= len(dna):     # Ending Procedure
            return ORFs

        codon = dna[p:p+3]

        if codon == 'ATG':
            ORFs.append(rest_of_ORF(dna[p:len(dna)]))
            p += len(ORFs[n])
            n += 1

        p += 3          # Shift position to new codon

    # TODO: implement this
    pass


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    ORFerLord = []
    ORF_F1 = find_all_ORFs_oneframe(dna)
    ORF_F2 = find_all_ORFs_oneframe(dna[1:])
    ORF_F3 = find_all_ORFs_oneframe(dna[2:])

    i = 0
    i1 = 0
    i2 = 0

    for ORF in ORF_F1:
        ORFerLord.append(ORF_F1[i])
        i += 1

    for ORF in ORF_F2:
        ORFerLord.append(ORF_F2[i1])
        i1 += 1

    for ORF in ORF_F1:
        ORFerLord.append(ORF_F3[i2])
        i2 += 1

    return ORFerLord

    # TODO: implement this
    pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    AND = get_reverse_complement(dna)       # Backwards DNA (AKA Compliment)

    ORF_dna = find_all_ORFs(dna)
    ORF_AND = find_all_ORFs(AND)
    MasterORFerLord = [ORF_dna, ORF_AND]

    return MasterORFerLord
    # TODO: implement this
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    #doctest.run_docstring_examples(find_all_ORFs_oneframe, globals(),
    # verbose=True)
