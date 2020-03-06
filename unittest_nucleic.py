# Testing class NucleicAcid

import unittest
from seq_align import NucleicAcid

class NucleicAcidMethods(unittest.TestCase):
    """ Test NucleicAcid methods"""

    def test_complement(self):
        foo = NucleicAcid('ACTG')
        foo2 = 'TGAC'
        self.assertEqual(foo.complement(), foo2)

    def test_codon(self):
        foo = NucleicAcid('GGTGCG')
        foo2 = ['GGT','GCG']
        self.assertEqual(foo.codon(), foo2)
#
    def test_translate(self):
        foo = NucleicAcid('GGTGCG')
        foo2 = 'GA'
        self.assertEqual(foo.translate(), foo2)

if __name__=='__main__':
    unittest.main()
