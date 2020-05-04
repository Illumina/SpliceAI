import argparse
import logging
from multiprocessing import Process, JoinableQueue, Queue, cpu_count

import pysam

from spliceai.utils import Annotator, get_delta_scores


try:
    from sys.stdin import buffer as std_in
    from sys.stdout import buffer as std_out
except ImportError:
    from sys import stdin as std_in
    from sys import stdout as std_out


def get_options():

    parser = argparse.ArgumentParser(description='Version: 1.3')
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
    args = parser.parse_args()

    return args


class Worker(Process):
    def __init__(self, records_queue, scores_queue, cli_args):
        Process.__init__(self)
        self.records_queue = records_queue
        self.scores_queue = scores_queue
        self.cli_args = cli_args

    def run(self):
        ann = Annotator(self.cli_args.R, self.cli_args.A)
        while True:
            record_data = self.records_queue.get()
            if record_data is None:
                self.records_queue.task_done()
                break

            scores = get_delta_scores(record_data, ann, self.cli_args.D, self.cli_args.M)
            self.scores_queue.put((record_data[0], scores))
            self.records_queue.task_done()


def main():
    args = get_options()

    if None in [args.I, args.O, args.D, args.M]:
        logging.error('Usage: spliceai [-h] [-I [input]] [-O [output]] -R reference -A annotation '
                      '[-D [distance]] [-M [mask]]')
        exit()

    try:
        vcf = pysam.VariantFile(args.I)
    except (IOError, ValueError) as e:
        logging.error('{}'.format(e))
        exit()

    records_queue = JoinableQueue()
    scores_queue = Queue(200)
    for i in range(cpu_count()):
        proc = Worker(records_queue, scores_queue, args)
        proc.start()

    index = 0
    for record in vcf:
        record_data = (index, record.chrom, record.pos, record.ref, record.alts, str(record))
        records_queue.put(record_data)
        index += 1
    vcf.close()

    for i in range(cpu_count()):
        records_queue.put(None)
    records_queue.join()

    scores = []
    while not scores_queue.empty():
        scores.append(scores_queue.get())
    scores = [score for index, score in sorted(scores, key=lambda item: item[0])]
    print(scores)

    vcf = pysam.VariantFile(args.I)
    header = vcf.header
    header.add_line('##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3.1 variant '
                    'annotation. These include delta scores (DS) and delta positions (DP) for '
                    'acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). '
                    'Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">')

    try:
        output = pysam.VariantFile(args.O, mode='w', header=header)
    except (IOError, ValueError) as e:
        logging.error('{}'.format(e))
        exit()

    index = 0
    for record in vcf:
        if len(scores[index]) > 0:
            record.info['SpliceAI'] = scores[index]
        output.write(record)
        index += 1

    vcf.close()
    output.close()


if __name__ == '__main__':
    main()
