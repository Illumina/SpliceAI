import argparse
import sys
import pysam
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

    vcf = pysam.VariantFile(args.vcf)
    header = vcf.header
    header.add_line('##INFO=<ID=SpliceAI,Number=.,Type=String,Description="Splice AI annotation for variant. These include delta scores (DS) for acceptor gain (AG), acceptor loss (AL), donor gain (DG) and donor loss (DL). The distance from the variant site to the splice site is also included (DP). Format: Allele|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">')
    output = pysam.VariantFile(args.out, mode='w', header=header)
    ann = annotator(args.R, args.A)

    for record in vcf:
        scores = get_delta_scores(record, ann)
        if len(scores) > 0:
            record.info['SpliceAI'] = scores
        output.write(record)


if __name__ == '__main__':
    main()

