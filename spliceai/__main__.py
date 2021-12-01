# Original source code modified to add prediction batching support by Invitae in 2021.
# Modifications copyright (c) 2021 Invitae Corporation.

import sys
import argparse
import logging
import pysam

from spliceai.batch.batch import VCFPredictionBatch
from spliceai.utils import Annotator, get_delta_scores

try:
    from sys.stdin import buffer as std_in
    from sys.stdout import buffer as std_out
except ImportError:
    from sys import stdin as std_in
    from sys import stdout as std_out


def get_options():

    parser = argparse.ArgumentParser(description='Version: 1.3.1')
    parser.add_argument('-I', metavar='input', nargs='?', default=std_in,
                        help='path to the input VCF file, defaults to standard in')
    parser.add_argument('-O', metavar='output', nargs='?', default=std_out,
                        help='path to the output VCF file, defaults to standard out')
    parser.add_argument('-R', metavar='reference', required=True,
                        help='path to the reference genome fasta file')
    parser.add_argument('-A', metavar='annotation', required=True,
                        help='"grch37" (GENCODE V24lift37 canonical annotation file in '
                             'package), "grch38" (GENCODE V24 canonical annotation file in '
                             'package), or path to a similar custom gene annotation file')
    parser.add_argument('-D', metavar='distance', nargs='?', default=50,
                        type=int, choices=range(0, 5000),
                        help='maximum distance between the variant and gained/lost splice '
                             'site, defaults to 50')
    parser.add_argument('-M', metavar='mask', nargs='?', default=0,
                        type=int, choices=[0, 1],
                        help='mask scores representing annotated acceptor/donor gain and '
                             'unannotated acceptor/donor loss, defaults to 0')
    parser.add_argument('-B', '--prediction-batch-size', metavar='prediction_batch_size', default=1, type=int,
                        help='number of predictions to process at a time, note a single vcf record '
                             'may have multiple predictions for overlapping genes and multiple alts')
    parser.add_argument('-T', '--tensorflow-batch-size', metavar='tensorflow_batch_size', type=int,
                        help='tensorflow batch size for model predictions')
    parser.add_argument('-V', '--verbose', action='store_true', help='enables verbose logging')
    args = parser.parse_args()

    return args


def main():
    args = get_options()

    if args.verbose:
        logging.basicConfig(
            format='%(asctime)s %(levelname)s %(name)s: - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            level=logging.DEBUG,
        )

    if None in [args.I, args.O, args.D, args.M]:
        logging.error('Usage: spliceai [-h] [-I [input]] [-O [output]] -R reference -A annotation '
                      '[-D [distance]] [-M [mask]] [-B [prediction_batch_size]] [-T [tensorflow_batch_size]]')
        exit()

    # Default the tensorflow batch size to the prediction_batch_size if it's not supplied in the args
    tensorflow_batch_size = args.tensorflow_batch_size if args.tensorflow_batch_size else args.prediction_batch_size

    run_spliceai(input_data=args.I, output_data=args.O, reference=args.R,
                 annotation=args.A, distance=args.D, mask=args.M,
                 prediction_batch_size=args.prediction_batch_size,
                 tensorflow_batch_size=tensorflow_batch_size)


def run_spliceai(input_data, output_data, reference, annotation, distance, mask, prediction_batch_size,
                 tensorflow_batch_size):

    try:
        vcf = pysam.VariantFile(input_data)
    except (IOError, ValueError) as e:
        logging.error('{}'.format(e))
        exit()

    header = vcf.header
    header.add_line('##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3.1 variant '
                    'annotation. These include delta scores (DS) and delta positions (DP) for '
                    'acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). '
                    'Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">')

    try:
        output_data = pysam.VariantFile(output_data, mode='w', header=header)
    except (IOError, ValueError) as e:
        logging.error('{}'.format(e))
        exit()

    ann = Annotator(reference, annotation)
    batch = None

    # Only use the batching code if we are batching
    if prediction_batch_size > 1:
        batch = VCFPredictionBatch(
            ann=ann,
            output=output_data,
            dist=distance,
            mask=mask,
            prediction_batch_size=prediction_batch_size,
            tensorflow_batch_size=tensorflow_batch_size,
        )

    for record in vcf:
        if batch:
            # Add record to batch, if batch fills, then they will all be processed at once
            batch.add_record(record)
        else:
            # If we're not batching, let's run the original code
            scores = get_delta_scores(record, ann, distance, mask)
            if len(scores) > 0:
                record.info['SpliceAI'] = scores
            output_data.write(record)

    if batch:
        # Ensure we process any leftover records in the batch when we finish iterating the VCF. This
        # would be a good candidate for a context manager if we removed the original non batching code above
        batch.finish()

    vcf.close()
    output_data.close()


if __name__ == '__main__':
    main()
