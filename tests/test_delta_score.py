from collections import namedtuple
import unittest
from pkg_resources import resource_filename
from spliceai.utils import Annotator, get_delta_scores


Record = namedtuple('Record', ['chrom', 'pos', 'ref', 'alts'])


class TestDeltaScore(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        fasta_path = resource_filename(__name__, 'data/test.fa')
        fasta_without_prefix_path = resource_filename(__name__, 'data/test_without_prefix.fa')
        cls.ann = Annotator(fasta_path, 'grch37')
        cls.ann_without_prefix = Annotator(fasta_without_prefix_path, 'grch37')

    def test_get_delta_score_acceptor(self):

        record = Record('10', 94077, 'A', ['C'])
        scores = get_delta_scores(record, self.ann, 500, 0)
        self.assertEqual(scores, ['C|TUBB8|0.15|0.27|0.00|0.05|89|-23|-267|193'])
        scores = get_delta_scores(record, self.ann_without_prefix, 500, 0)
        self.assertEqual(scores, ['C|TUBB8|0.15|0.27|0.00|0.05|89|-23|-267|193'])

        record = Record('chr10', 94077, 'A', ['C'])
        scores = get_delta_scores(record, self.ann, 500, 0)
        self.assertEqual(scores, ['C|TUBB8|0.15|0.27|0.00|0.05|89|-23|-267|193'])
        scores = get_delta_scores(record, self.ann_without_prefix, 500, 0)
        self.assertEqual(scores, ['C|TUBB8|0.15|0.27|0.00|0.05|89|-23|-267|193'])

    def test_get_delta_score_donor(self):

        record = Record('10', 94555, 'C', ['T'])
        scores = get_delta_scores(record, self.ann, 500, 0)
        self.assertEqual(scores, ['T|TUBB8|0.01|0.18|0.15|0.62|-2|110|-190|0'])
        scores = get_delta_scores(record, self.ann_without_prefix, 500, 0)
        self.assertEqual(scores, ['T|TUBB8|0.01|0.18|0.15|0.62|-2|110|-190|0'])

        record = Record('chr10', 94555, 'C', ['T'])
        scores = get_delta_scores(record, self.ann, 500, 0)
        self.assertEqual(scores, ['T|TUBB8|0.01|0.18|0.15|0.62|-2|110|-190|0'])
        scores = get_delta_scores(record, self.ann_without_prefix, 500, 0)
        self.assertEqual(scores, ['T|TUBB8|0.01|0.18|0.15|0.62|-2|110|-190|0'])
