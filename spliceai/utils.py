# Original source code modified to add prediction batching support by Invitae in 2021.
# Modifications copyright (c) 2021 Invitae Corporation.

import collections

from pkg_resources import resource_filename
import pandas as pd
import numpy as np
from pyfaidx import Fasta
from keras.models import load_model
import logging


GeneInfo = collections.namedtuple('GeneInfo', 'genes strands idxs')


class Annotator:

    def __init__(self, ref_fasta, annotations):

        if annotations == 'grch37':
            annotations = resource_filename(__name__, 'annotations/grch37.txt')
        elif annotations == 'grch38':
            annotations = resource_filename(__name__, 'annotations/grch38.txt')

        try:
            df = pd.read_csv(annotations, sep='\t', dtype={'CHROM': object})
            self.genes = df['#NAME'].to_numpy()
            self.chroms = df['CHROM'].to_numpy()
            self.strands = df['STRAND'].to_numpy()
            self.tx_starts = df['TX_START'].to_numpy()+1
            self.tx_ends = df['TX_END'].to_numpy()
            self.exon_starts = [np.asarray([int(i) for i in c.split(',') if i])+1
                                for c in df['EXON_START'].to_numpy()]
            self.exon_ends = [np.asarray([int(i) for i in c.split(',') if i])
                              for c in df['EXON_END'].to_numpy()]
        except IOError as e:
            logging.error('{}'.format(e))
            exit()
        except (KeyError, pd.errors.ParserError) as e:
            logging.error('Gene annotation file {} not formatted properly: {}'.format(annotations, e))
            exit()

        try:
            self.ref_fasta = Fasta(ref_fasta, rebuild=False)
        except IOError as e:
            logging.error('{}'.format(e))
            exit()

        paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
        self.models = [load_model(resource_filename(__name__, x)) for x in paths]

    def get_name_and_strand(self, chrom, pos):

        chrom = normalise_chrom(chrom, list(self.chroms)[0])
        idxs = np.intersect1d(np.nonzero(self.chroms == chrom)[0],
                              np.intersect1d(np.nonzero(self.tx_starts <= pos)[0],
                              np.nonzero(pos <= self.tx_ends)[0]))

        if len(idxs) >= 1:
            return GeneInfo(genes=self.genes[idxs], strands=self.strands[idxs], idxs=idxs)
        else:
            return GeneInfo(genes=[], strands=[], idxs=[])

    def get_pos_data(self, idx, pos):

        dist_tx_start = self.tx_starts[idx]-pos
        dist_tx_end = self.tx_ends[idx]-pos
        dist_exon_bdry = min(np.union1d(self.exon_starts[idx], self.exon_ends[idx])-pos, key=abs)
        dist_ann = (dist_tx_start, dist_tx_end, dist_exon_bdry)

        return dist_ann


def one_hot_encode(seq):

    map = np.asarray([[0, 0, 0, 0],
                      [1, 0, 0, 0],
                      [0, 1, 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])

    seq = seq.upper().replace('A', '\x01').replace('C', '\x02')
    seq = seq.replace('G', '\x03').replace('T', '\x04').replace('N', '\x00')

    return map[np.fromstring(seq, np.int8) % 5]


def normalise_chrom(source, target):

    def has_prefix(x):
        return x.startswith('chr')

    if has_prefix(source) and not has_prefix(target):
        return source.strip('chr')
    elif not has_prefix(source) and has_prefix(target):
        return 'chr'+source

    return source


def get_cov(dist_var):
    return 2 * dist_var + 1


def get_wid(cov):
    return 10000 + cov


def encode_seqs(record, seq, ann, gene_info, gene_ix, alt_ix, wid):
    dist_ann = ann.get_pos_data(gene_info.idxs[gene_ix], record.pos)
    pad_size = [max(wid // 2 + dist_ann[0], 0), max(wid // 2 - dist_ann[1], 0)]
    ref_len = len(record.ref)
    x_ref = 'N' * pad_size[0] + seq[pad_size[0]:wid - pad_size[1]] + 'N' * pad_size[1]
    x_alt = x_ref[:wid // 2] + str(record.alts[alt_ix]) + x_ref[wid // 2 + ref_len:]

    x_ref = one_hot_encode(x_ref)[None, :]
    x_alt = one_hot_encode(x_alt)[None, :]

    if gene_info.strands[gene_ix] == '-':
        x_ref = x_ref[:, ::-1, ::-1]
        x_alt = x_alt[:, ::-1, ::-1]

    return x_ref, x_alt


def is_record_valid(record):
    try:
        record.chrom, record.pos, record.ref, len(record.alts)
    except TypeError:
        logging.warning('Skipping record (bad input): {}'.format(record))
        return False
    return True


def get_seq(record, ann, wid):
    chrom = normalise_chrom(record.chrom, list(ann.ref_fasta.keys())[0])
    try:
        seq = ann.ref_fasta[chrom][
            record.pos - wid // 2 - 1: record.pos + wid // 2
        ].seq
    except (IndexError, ValueError):
        logging.warning('Skipping record (fasta issue): {}'.format(record))
        return ""
    return seq


def is_valid_alt_record(record, alt_ix):
    if '.' in record.alts[alt_ix] or '-' in record.alts[alt_ix] or '*' in record.alts[alt_ix]:
        return False
    if '<' in record.alts[alt_ix] or '>' in record.alts[alt_ix]:
        return False
    return True


def is_location_predictable(record, seq, wid, dist_var):
    var_ref_seq = seq[wid // 2: wid // 2 + len(record.ref)].upper()
    if var_ref_seq != record.ref:
        logging.warning(
            'Skipping record (ref issue) should be {} - {}'.format(var_ref_seq, record)
        )
        return False

    if len(seq) != wid:
        logging.warning('Skipping record (near chromosome end): {}'.format(record))
        return False

    if len(record.ref) > 2 * dist_var:
        logging.warning('Skipping record (ref too long): {}'.format(record))
        return False

    return True


def create_unhandled_delta_score(alt, gene):
    return '{}|{}|.|.|.|.|.|.|.|.'.format(alt, gene)


def get_delta_scores(record, ann, dist_var, mask):

    cov = get_cov(dist_var)
    wid = get_wid(cov)
    delta_scores = []

    if not is_record_valid(record):
        return delta_scores

    gene_info = ann.get_name_and_strand(record.chrom, record.pos)
    if len(gene_info.idxs) == 0:
        return delta_scores

    seq = get_seq(record, ann, wid)
    if not seq:
        return delta_scores

    if not is_location_predictable(record, seq, wid, dist_var):
        return delta_scores

    for alt_ix in range(len(record.alts)):
        for gene_ix in range(len(gene_info.idxs)):

            if not is_valid_alt_record(record, alt_ix):
                continue

            if len(record.ref) > 1 and len(record.alts[alt_ix]) > 1:
                delta_score = create_unhandled_delta_score(record.alts[alt_ix], gene_info.genes[gene_ix])
                delta_scores.append(delta_score)
                continue

            x_ref, x_alt = encode_seqs(record=record,
                                       seq=seq,
                                       ann=ann,
                                       gene_info=gene_info,
                                       gene_ix=gene_ix,
                                       alt_ix=alt_ix,
                                       wid=wid)

            y_ref = np.mean([ann.models[m].predict(x_ref) for m in range(5)], axis=0)
            y_alt = np.mean([ann.models[m].predict(x_alt) for m in range(5)], axis=0)

            delta_score = get_alt_gene_delta_score(record=record,
                                                   ann=ann,
                                                   alt_ix=alt_ix,
                                                   gene_ix=gene_ix,
                                                   y_ref=y_ref,
                                                   y_alt=y_alt,
                                                   cov=cov,
                                                   gene_info=gene_info,
                                                   mask=mask)
            delta_scores.append(delta_score)

    return delta_scores


def get_alt_gene_delta_score(record, ann, alt_ix, gene_ix, y_ref, y_alt, cov, gene_info, mask):
    ref_len = len(record.ref)
    dist_ann = ann.get_pos_data(gene_info.idxs[gene_ix], record.pos)
    alt_len = len(record.alts[alt_ix])
    del_len = max(ref_len - alt_len, 0)

    if gene_info.strands[gene_ix] == '-':
        y_ref = y_ref[:, ::-1]
        y_alt = y_alt[:, ::-1]

    if ref_len > 1 and alt_len == 1:
        y_alt = np.concatenate([
            y_alt[:, :cov // 2 + alt_len],
            np.zeros((1, del_len, 3)),
            y_alt[:, cov // 2 + alt_len:]],
            axis=1)
    elif ref_len == 1 and alt_len > 1:
        y_alt = np.concatenate([
            y_alt[:, :cov // 2],
            np.max(y_alt[:, cov // 2:cov // 2 + alt_len], axis=1)[:, None, :],
            y_alt[:, cov // 2 + alt_len:]],
            axis=1)

    y = np.concatenate([y_ref, y_alt])

    idx_pa = (y[1, :, 1] - y[0, :, 1]).argmax()
    idx_na = (y[0, :, 1] - y[1, :, 1]).argmax()
    idx_pd = (y[1, :, 2] - y[0, :, 2]).argmax()
    idx_nd = (y[0, :, 2] - y[1, :, 2]).argmax()

    mask_pa = np.logical_and((idx_pa - cov // 2 == dist_ann[2]), mask)
    mask_na = np.logical_and((idx_na - cov // 2 != dist_ann[2]), mask)
    mask_pd = np.logical_and((idx_pd - cov // 2 == dist_ann[2]), mask)
    mask_nd = np.logical_and((idx_nd - cov // 2 != dist_ann[2]), mask)

    delta_score = "{}|{}|{:.2f}|{:.2f}|{:.2f}|{:.2f}|{}|{}|{}|{}".format(
        record.alts[alt_ix],
        gene_info.genes[gene_ix],
        (y[1, idx_pa, 1] - y[0, idx_pa, 1]) * (1 - mask_pa),
        (y[0, idx_na, 1] - y[1, idx_na, 1]) * (1 - mask_na),
        (y[1, idx_pd, 2] - y[0, idx_pd, 2]) * (1 - mask_pd),
        (y[0, idx_nd, 2] - y[1, idx_nd, 2]) * (1 - mask_nd),
        idx_pa - cov // 2,
        idx_na - cov // 2,
        idx_pd - cov // 2,
        idx_nd - cov // 2)

    return delta_score
