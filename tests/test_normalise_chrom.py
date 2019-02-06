
import unittest
from spliceai.normalise_chrom import normalise_chrom


class TestNormaliseChrom(unittest.TestCase):
    def test_normalise_chrom(self):
        ''' test normalise_chrom
        '''
        self.assertEqual(normalise_chrom('1', '1'), '1')
        self.assertEqual(normalise_chrom('1', '2'), '1')
        self.assertEqual(normalise_chrom('1', 'chr2'), 'chr1')
        self.assertEqual(normalise_chrom('chr1', 'chr2'), 'chr1')
        self.assertEqual(normalise_chrom('chr1', '2'), '1')
