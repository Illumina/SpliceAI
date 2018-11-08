import os
import argparse
import vcf
from spliceai.utils import annotator, get_delta_scores


def main():

    cwd = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser()
    parser.add_argument('-I')
    parser.add_argument('-O')
    parser.add_argument('-R')
    parser.add_argument('-A', default=cwd+'/annotations/GENCODE.v24lift37')
    args = parser.parse_args()

    vcf_reader = vcf.Reader(open(args.I, 'r'))
    vcf_writer = vcf.Writer(open(args.O, 'w'), vcf_reader)
    ann = annotator(args.R, args.A)

    for record in vcf_reader:

        delta_scores = get_delta_scores(record, ann)
        record.add_info('SpliceAI', delta_scores)
        vcf_writer.write_record(record)


if __name__ == '__main__':
    main()

