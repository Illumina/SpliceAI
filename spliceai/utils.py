import pandas as pd
import numpy as np
import re
from pkg_resources import resource_filename
import pyfasta
from keras.models import load_model


class annotator():

    def __init__(self, ref_fasta, annotations):

        if annotations is None:
            annotations = resource_filename(__name__, 'annotations/GENCODE.v24lift37')
        df = pd.read_csv(annotations, sep='\t')

        self.genes = df['#NAME'].get_values()
        self.chroms = df['CHROM'].get_values()
        self.strands = df['STRAND'].get_values()
        self.tx_starts = df['TX_START'].get_values()+1
        self.tx_ends = df['TX_END'].get_values()

        self.ref_fasta = pyfasta.Fasta(ref_fasta)

        paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
        self.models = [load_model(resource_filename(__name__, x)) for x in paths]

    def get_name_and_strand(self, chrom, pos):

        idxs = np.intersect1d(
                   np.nonzero(self.chroms == chrom)[0],
                   np.intersect1d(np.nonzero(self.tx_starts <= pos)[0],
                   np.nonzero(pos <= self.tx_ends)[0]))

        if len(idxs) >= 1:
            return (self.genes[idxs], self.strands[idxs], idxs)
        else:
            return ([], [], [])

    def get_pos_data(self, idx, pos):

        dist_tx_start = self.tx_starts[idx]-pos
        dist_tx_end = self.tx_ends[idx]-pos

        dist = (dist_tx_start, dist_tx_end)

        return dist


def one_hot_encode(seq):

    MAP = np.asarray([[0, 0, 0, 0],
                      [1, 0, 0, 0],
                      [0, 1, 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])

    seq = seq.upper().replace('A', '\x01').replace('C', '\x02')
    seq = seq.replace('G', '\x03').replace('T', '\x04').replace('N', '\x00')

    return MAP[np.fromstring(seq, np.int8) % 5]


def get_delta_scores(record, ann, L=1001):

    W = 10000+L

    (genes, strands, idxs) = ann.get_name_and_strand(record.CHROM, record.POS)

    delta_scores = []

    for j in range(len(record.ALT)):
        for i in range(len(idxs)):

            dist = ann.get_pos_data(idxs[i], record.POS)

            if len(record.REF) > 1 and len(record.ALT[j]) > 1:
                continue
            # Ignoring complicated INDELs
            
            pad_size = [max(W//2+dist[0], 0), max(W//2-dist[1], 0)]
            ref_len = len(record.REF)
            alt_len = len(record.ALT[j])
            del_len = max(ref_len-alt_len, 0)

            seq = ann.ref_fasta['chr'+str(record.CHROM)][
                                record.POS-W//2-1:record.POS+W//2]
            x_ref = 'N'*pad_size[0]+seq[pad_size[0]:W-pad_size[1]]\
                     +'N'*pad_size[1]
            x_alt = x_ref[:W//2]+str(record.ALT[j])+x_ref[W//2+ref_len:]

            X_ref = one_hot_encode(x_ref)[None, :]
            X_alt = one_hot_encode(x_alt)[None, :]

            if strands[i] == '-':
                X_ref = X_ref[:, ::-1, ::-1]
                X_alt = X_alt[:, ::-1, ::-1]

            Y0 = np.asarray(ann.models[0].predict(X_ref))
            Y1 = np.asarray(ann.models[1].predict(X_ref))
            Y2 = np.asarray(ann.models[2].predict(X_ref))
            Y3 = np.asarray(ann.models[3].predict(X_ref))
            Y4 = np.asarray(ann.models[4].predict(X_ref))
            Y_ref = (Y0+Y1+Y2+Y3+Y4)/5

            Y0 = np.asarray(ann.models[0].predict(X_alt))
            Y1 = np.asarray(ann.models[1].predict(X_alt))
            Y2 = np.asarray(ann.models[2].predict(X_alt))
            Y3 = np.asarray(ann.models[3].predict(X_alt))
            Y4 = np.asarray(ann.models[4].predict(X_alt))
            Y_alt = (Y0+Y1+Y2+Y3+Y4)/5

            if strands[i] == '-':
                Y_ref = Y_ref[:, ::-1]
                Y_alt = Y_alt[:, ::-1]

            if ref_len > 1 and alt_len == 1:
                Y_alt = np.concatenate([Y_alt[:, :L//2+alt_len],
                                        np.zeros((1, del_len, 3)),
                                        Y_alt[:, L//2+alt_len:]], axis=1)
            elif ref_len == 1 and alt_len > 1:
                Y_alt = np.concatenate([
                    Y_alt[:, :L//2],
                    np.max(Y_alt[:, L//2:L//2+alt_len], axis=1)[:, None, :],
                    Y_alt[:, L//2+alt_len:]], axis=1)

            Y = np.concatenate([Y_ref, Y_alt])

            idx_pA = (Y[1, :, 1]-Y[0, :, 1]).argmax()
            idx_nA = (Y[0, :, 1]-Y[1, :, 1]).argmax()
            idx_pD = (Y[1, :, 2]-Y[0, :, 2]).argmax()
            idx_nD = (Y[0, :, 2]-Y[1, :, 2]).argmax()

            delta_scores.append(
                "{}|{}|{:.2f}|{:.2f}|{:.2f}|{:.2f}|{}|{}|{}|{}".format(
                record.ALT[j],
                genes[i],
                Y[1, idx_pA, 1]-Y[0, idx_pA, 1],
                Y[0, idx_nA, 1]-Y[1, idx_nA, 1],
                Y[1, idx_pD, 2]-Y[0, idx_pD, 2],
                Y[0, idx_nD, 2]-Y[1, idx_nD, 2],
                idx_pA-L//2,
                idx_nA-L//2,
                idx_pD-L//2,
                idx_nD-L//2))

    return delta_scores
