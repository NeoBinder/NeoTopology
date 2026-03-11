import os
import tempfile

from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def muscle_alignment(seqs, muscle_exe="muscle"):
    """Align 2 sequences with muscle"""
    with tempfile.NamedTemporaryFile(delete=False) as temp:
        temp_name = temp.name
    SeqIO.write(seqs, temp_name, "fasta")
    name = os.path.splitext(temp_name)[0]
    cline = MuscleCommandline(muscle_exe, input=temp_name, out=name + ".txt")
    stdout, stderr = cline()
    return AlignIO.read(name + ".txt", "fasta")


def build_seq_record_from_seqs(seqs_dict):
    """
    return sequence record list from sequence dictionary
    """
    seq_record_ls = []
    for k, each_seq in seqs_dict.items():
        seq_record_ls.append(SeqRecord(Seq(each_seq), id=k))
    return seq_record_ls
