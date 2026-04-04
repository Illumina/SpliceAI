# Original source code modified to add prediction batching support by Invitae in 2021.
# Modifications copyright (c) 2021 Invitae Corporation.

import logging

from spliceai.utils import get_alt_gene_delta_score, is_record_valid, get_seq, \
    is_location_predictable, get_cov, get_wid, is_valid_alt_record, encode_seqs, create_unhandled_delta_score

logger = logging.getLogger(__name__)


def get_preds(ann, x, batch_size=32):
    logger.debug('Running get_preds with matrix size: {}'.format(x.shape))
    return [
        ann.models[m].predict(x, batch_size=batch_size, verbose=0) for m in range(5)
    ]


# Heavily based on utils.get_delta_scores but only handles the validation and encoding
# of the record, but doesn't do any of the prediction or post-processing steps
def encode_batch_records(record, ann, dist_var, gene_info):
    cov = get_cov(dist_var)
    wid = get_wid(cov)
    # If the record is not going to get a prediction, return this empty encoding
    empty_encoding = ([], [])

    if not is_record_valid(record):
        return empty_encoding

    seq = get_seq(record, ann, wid)
    if not seq:
        return empty_encoding

    if not is_location_predictable(record, seq, wid, dist_var):
        return empty_encoding

    all_x_ref = []
    all_x_alt = []
    for alt_ix in range(len(record.alts)):
        for gene_ix in range(len(gene_info.idxs)):

            if not is_valid_alt_record(record, alt_ix):
                continue

            x_ref, x_alt = encode_seqs(record=record,
                                       seq=seq,
                                       ann=ann,
                                       gene_info=gene_info,
                                       gene_ix=gene_ix,
                                       alt_ix=alt_ix,
                                       wid=wid)

            all_x_ref.append(x_ref)
            all_x_alt.append(x_alt)

    return all_x_ref, all_x_alt


# Heavily based on utils.get_delta_scores but only handles the post-processing steps after
# the models have made the predictions
def extract_delta_scores(
    all_y_ref, all_y_alt, record, ann, dist_var, mask, gene_info
):
    cov = get_cov(dist_var)
    delta_scores = []
    pred_ix = 0
    for alt_ix in range(len(record.alts)):
        for gene_ix in range(len(gene_info.idxs)):

            # Pull prediction out of batch
            y_ref = all_y_ref[pred_ix]
            y_alt = all_y_alt[pred_ix]

            # No prediction here
            if y_ref is None or y_alt is None:
                continue

            if not is_valid_alt_record(record, alt_ix):
                continue

            if len(record.ref) > 1 and len(record.alts[alt_ix]) > 1:
                pred_ix += 1
                delta_score = create_unhandled_delta_score(record.alts[alt_ix], gene_info.genes[gene_ix])
                delta_scores.append(delta_score)
                continue

            if pred_ix >= len(all_y_ref) or pred_ix >= len(all_y_alt):
                raise LookupError(
                    'Prediction index {} does not exist in prediction matrices: ref({}) alt({})'.format(
                        pred_ix, len(all_y_ref), len(all_y_alt)
                    )
                )

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

            pred_ix += 1

    return delta_scores
