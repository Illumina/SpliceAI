import argparse
import sys
import vcf
from spliceai.utils import annotator, get_delta_scores

def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', nargs='?', default=sys.stdin.buffer,
        help='path to VCF input, defaults to standard in')
    parser.add_argument('out', nargs='?', default=sys.stdout,
        help='path to write output VCF to, defaults to standard out.')
    parser.add_argument('-R', required=True,
        help='path to genome fasta file (must be fai indexed for quick access)')
    parser.add_argument('-A',
        help='path to gencode genes file (defaults to file in package)')
    args = parser.parse_args()

    try:
        args.vcf = open(args.vcf, 'rt')
    except TypeError:
        pass

    try:
        args.out = open(args.out, 'wt')
    except TypeError:
        pass

    return args

def main():
    args = get_options()

    vcf_reader = vcf.Reader(args.vcf)
    vcf_writer = vcf.Writer(args.out, vcf_reader)
    ann = annotator(args.R, args.A)

    for record in vcf_reader:

        delta_scores = get_delta_scores(record, ann)
        record.add_info('SpliceAI', delta_scores)
        vcf_writer.write_record(record)


if __name__ == '__main__':
    main()

