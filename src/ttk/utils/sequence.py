import math
from Bio.Align import PairwiseAligner

aligner = PairwiseAligner()
# filter functions
nuc_alphabet = set("ATCGU".lower())


def isNuc(seq):
    return len(set(seq.seq.lower()).difference(nuc_alphabet)) == 0


def getSeqIdentity(seq1, seq2):
    align_result = aligner.align(seq1, seq2)
    align_score = align_result.score
    return align_score / (len(seq1) + len(seq2)) * 2


def analyseSeq(seq1, seq2):
    lens1, lens2 = len(seq1), len(seq2)
    identity_score = getSeqIdentity(seq1, seq2)
    # cos angle - sqrt(2)/2  which is cos45
    length_diversity = abs(lens1 / math.sqrt(lens1**2 + lens2**2) - 0.7071067811865475)
    return identity_score, length_diversity


def analyseDimer(pdbid, record, sequence_similarity_threshold=0.95):
    assert len(record) == 2
    lens1, lens2 = len(record[0]), len(record[1])
    identity, diversity = analyseSeq(*record)
    # parse structure
    return len(record), lens1, lens2, identity, diversity
