from pkg_resources import resource_filename
import unittest
from collections import namedtuple
from spliceai.utils import annotator, get_delta_scores

Record = namedtuple('Record', ['chrom', 'pos', 'ref', 'alts'])

class TestDeltaScore(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        fasta_path = resource_filename(__name__, 'data/test.fa')
        self.ann = annotator(fasta_path, 'grch37')

    def test_get_delta_score_acceptor(self):
        ''' test get_delta_scores for a predicted acceptor
        '''
        record = Record('10', 94077, 'A', ['C'])
        scores = get_delta_scores(record, self.ann)
        self.assertEqual(scores , ['C|TUBB8|0.15|0.27|0.00|0.05|89|-23|-267|193'])
    
    def test_get_delta_score_donor(self):
        ''' test get_delta_scores for a predicted donor
        '''
        record = Record('10', 94555, 'C', ['T'])
        scores = get_delta_scores(record, self.ann)
        self.assertEqual(scores , ['T|TUBB8|0.01|0.18|0.15|0.62|-2|110|-190|0'])

    def test_get_delta_score_mismatched_prefix(self):
        ''' test get_delta_scores when the chromosome prefixes don't match the fasta
        '''
        record = Record('chr10', 94555, 'C', ['T'])

        # check when the fasta lacks chr prefixes, but the variant has them
        fasta_path = resource_filename(__name__, 'data/test_without_prefix.fa')
        ann = annotator(fasta_path, None)
        scores = get_delta_scores(record, ann)
        self.assertEqual(scores , ['T|TUBB8|0.01|0.18|0.15|0.62|-2|110|-190|0'])

        # check it works when fasta and variant both have chr prefix
        scores = get_delta_scores(record, self.ann)
        self.assertEqual(scores , ['T|TUBB8|0.01|0.18|0.15|0.62|-2|110|-190|0'])
