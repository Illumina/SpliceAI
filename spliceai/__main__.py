import argparse
import sys
import pysam
from spliceai.utils import Annotator, get_delta_scores


try:
    from sys.stdin import buffer as std_in
    from sys.stdout import buffer as std_out
except ImportError:
    from sys import stdin as std_in
    from sys import stdout as std_out


def get_options():

    parser = argparse.ArgumentParser()
    parser.add_argument('-I', nargs='?', default=std_in,
                        help='path to the input VCF file, defaults to standard in')
    parser.add_argument('-O', nargs='?', default=std_out,
                        help='path to the output VCF file, defaults to standard out')
    parser.add_argument('-R', required=True,
                        help='path to the genome fasta file')
    parser.add_argument('-A', required=True,
                        help='"grch37" (uses GENCODE canonical annotation file in package), '
                             '"grch38" (uses GENCODE canonical annotation file in package), '
                             'or path to a similarly-constructed custom gene annotation file')
    parser.add_argument('--batch_size', type=int, default=32,
                        help='batch_size: Integer. If unspecified, it will default to 32.')
    args = parser.parse_args()

    try:
        args.I = open(args.I, 'rt')
    except TypeError:
        pass

    try:
        args.O = open(args.O, 'wt')
    except TypeError:
        pass

    return args


def main():

    args = get_options()

    vcf = pysam.VariantFile(args.I)
    header = vcf.header
    header.add_line('##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.2 variant '
                    'annotation. These include delta scores (DS) and delta positions (DP) for '
                    'acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). '
                    'Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">')
    output = pysam.VariantFile(args.O, mode='w', header=header)
    ann = Annotator(args.R, args.A, args.batch_size)

    for record in vcf:

        scores = get_delta_scores(record, ann)
        if len(scores) > 0:
            record.info['SpliceAI'] = scores
        output.write(record)


if __name__ == '__main__':
    main()
