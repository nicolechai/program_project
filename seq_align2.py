## uniprot
# first wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebasse/

from Bio import SeqIO
import gzip

import unittest

# input sequence

# make into simplified Seq.Record class
# ?class Base

class MySeqRecord():
    def __init__(self, *bases):
        self.items = []
        for i in items:
            self.append(i)

# append only if DNA sequence
    def append(self,item):
        if not isinstance(item, Base):
            raise TypeError("Not a DNA sequence")
        self.items.append(item)

# allow indexing
    def __getitem__(self,index):
        if isinstance(self,index):
            return self.items[index]
        else:
            return MySeqRecord(*self.items[index])

#
    def __setitem__(self,index,item):
        if not isinstance(item, Base):
            raise TypeError("Not a DNA sequence")
        self.items[index] = item
# print
    def __str__(self):
        return str(self.items)

# object representation
    def __repr__(self):
        s = ", "
        s = s.join( [repr(x) for x in self.items] )
        return "MySeqRecord({})".format(s)


# blast sequence

## uniprot and rhodopsin
LIMIT_RESULTS = 20
rhodopsin = []
count = 0

#with gzip.open
# with gzip.open("uniprot_sprot.xml.gz") as f:
#     for entry in SeqIO.parse(f, "uniprot-xml"):
#         found = False
#         if ("gene_name_primary" in entry.annotations.keys()):
#             value = entry.annotations["gene_name_primary"].lower()
#             if ("rho" in value):
#                 print(value)
#                 found = value
#             print("Id:{} Primary:{}".format(entry.id, value)) if found else ""
#         if (not found and "gene_name_synonym" in entry.annotations.keys()):
#             values = [y.lower() for y in entry.annotations["gene_name_synonym"]]
#             for value in values:
#                 if ("rho" in value):
#                     found = True
#                     break
#             print("Id:{}, Secondary:{}".format(entry.id, str(values))) if found else ""
#         if (found):
#             rhodopsin.append(entry)
#             count+=1
#         if (count == LIMIT_RESULTS):
#             break
#
# ## output results to file
# with open("rhodopsin.fasta", 'w') as f:
#     SeqIO.write(rhodopsin, f, "fasta")

## alignment
# from Bio.Align.Applications import ClustalwCommandline
# # create bash command
# cline = ClustalwCommandline("clustalw", infile="rhodopsin.fasta")
# print(cline)
# ##execute Clustalw
# stdout, stderr = cline()

from Bio.Align.Applications import MuscleCommandline
fasta_rho = MuscleCommandline(input="rhodopsin.fasta", out = "rhodopsin.html", html=True)
fasta_rho()
cline = MuscleCommandline(input="rhodopsin.fasta")
stdout, stderr = cline()
from io import StringIO
from Bio import AlignIO
align = AlignIO.read(StringIO(stdout), "fasta")
print(align)

# Summary info
from Bio.Align import AlignInfo

summary_align = AlignInfo.SummaryInfo(align)
# AlignInfo.SummaryInfo.replacement_dictionary(summary_align)
consensus = summary_align.dumb_consensus()
print(consensus)

my_pssm = summary_align.pos_specific_score_matrix(consensus, chars_to_ignore=["N"])
print(my_pssm)

# commandline wrapper for secondary structure prediction
# seq.structure object storing it

# Testing
# class TestStringMethods(unittest.TestCase):
#
#     def test_upper(self):
#         self.assertEqual('foo'.upper(), 'FOO')
#
#     def test_isupper(self):
#         self.assertTrue('FOO'.isupper())
#         self.assertFalse('Foo'.isupper())
#
#     def test_split(self):
#         s = 'hello world'
#         self.assertEqual(s.split(), ['hello', 'world'])
#         # check that s.split fails when the separator is not a string
#         with self.assertRaises(TypeError):
#             s.split(2)

if __name__=='__main__':
    unittest.main()
    # ?"Equal? {}".format( seq == list(range(10))))
