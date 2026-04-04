import difflib
import sys
import tempfile
import unittest

from pkg_resources import resource_filename

from spliceai.__main__ import run_spliceai


class TestBatching(unittest.TestCase):

    def run_test(self, prediction_batch_size):
        # This reference genome was built extracting a small sequence set from the reference genome
        #   samtools faidx hg19.fa chr19:0-1500000 | awk '{sub(/chr19.*/,"chr19")}1' > tests/data/chr19_small.fa
        #   samtools faidx tests/data/chr19_small.fa
        ref_genome = resource_filename(__name__, "data/chr19_grch37_small.fa")
        genome_version = "grch37"

        input_file = resource_filename(__name__, "fixtures/batch_test.vcf")
        expected_file = resource_filename(__name__, "fixtures/batch_test_expected.vcf")

        with tempfile.NamedTemporaryFile() as tf:
            output_file = tf.name
            run_spliceai(input_data=input_file, output_data=output_file, reference=ref_genome,
                         annotation=genome_version, distance=50, mask=0,
                         prediction_batch_size=prediction_batch_size, tensorflow_batch_size=32)
            with open(output_file) as fh:
                batch_file_contents = fh.readlines()

        with open(expected_file) as fh:
            expected_file_contents = fh.readlines()

        # If there isn't a match, output a nice diff of the file diffs
        if expected_file_contents != batch_file_contents:
            sys.stdout.writelines(
                difflib.unified_diff(expected_file_contents, batch_file_contents)
            )
        assert expected_file_contents == batch_file_contents

    # This test will run everything without batching
    def test_batch_size_1(self):
        self.run_test(1)

    def test_batch_size_2(self):
        self.run_test(2)

    def test_batch_size_3(self):
        self.run_test(3)

    # This test should run everything in one giant batch
    def test_batch_size_30(self):
        self.run_test(30)
