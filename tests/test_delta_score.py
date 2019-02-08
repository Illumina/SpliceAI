from pkg_resources import resource_filename
import unittest
from spliceai.utils import annotator, get_delta_scores


class TestDeltaScore(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        fasta_path = resource_filename(__name__, 'data/test.fa')
        self.ann = annotator(fasta_path, 'grch37')

    def test_get_delta_score_acceptor(self):
        ''' test get_delta_scores for a predicted acceptor
        '''
        class Record():
            chrom, pos, ref, alts = '10', 94077, 'A', ['C']

        record = Record()
        scores = get_delta_scores(record, self.ann)
        self.assertEqual(scores , ['C|TUBB8|0.15|0.27|0.00|0.05|89|-23|-267|193'])
    
    def test_get_delta_score_donor(self):
        ''' test get_delta_scores for a predicted donor
        '''
        class Record():
            chrom, pos, ref, alts = '10', 94555, 'C', ['T']

        record = Record()
        scores = get_delta_scores(record, self.ann)
        self.assertEqual(scores , ['T|TUBB8|0.01|0.18|0.15|0.62|-2|110|-190|0'])
