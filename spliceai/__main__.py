import argparse
import sys
import pysam
from spliceai.utils import annotator, get_delta_scores


def get_options():

    parser = argparse.ArgumentParser()
    parser.add_argument('-I', nargs='?', default=sys.stdin,
        help='path to the input VCF file, defaults to standard in')
    parser.add_argument('-O', nargs='?', default=sys.stdout,
        help='path to the output VCF file, defaults to standard out')
    parser.add_argument('-R', required=True,
        help='path to the genome fasta file')
    parser.add_argument('-A',
        help='path to the gene annotations file, defaults to file in package')
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
    header.add_line('##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAI variant annotation. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: Allele|Gene|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">')
    output = pysam.VariantFile(args.O, mode='w', header=header)
    ann = annotator(args.R, args.A)

    for record in vcf:
        scores = get_delta_scores(record, ann)
        if len(scores) > 0:
            record.info['SpliceAI'] = scores
        output.write(record)


if __name__ == '__main__':
    main()

