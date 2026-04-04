# Original source code modified to add prediction batching support by Invitae in 2021.
# Modifications copyright (c) 2021 Invitae Corporation.

import collections
import logging
import time

import numpy as np

from spliceai.batch.batch_utils import extract_delta_scores, get_preds, encode_batch_records

logger = logging.getLogger(__name__)

SequenceType_REF = 0
SequenceType_ALT = 1


BatchLookupIndex = collections.namedtuple(
    'BatchLookupIndex', 'sequence_type tensor_size batch_index'
)

PreparedVCFRecord = collections.namedtuple(
    'PreparedVCFRecord', 'vcf_record gene_info locations'
)


class VCFPredictionBatch:
    def __init__(self, ann, output, dist, mask, prediction_batch_size, tensorflow_batch_size):
        self.ann = ann
        self.output = output
        self.dist = dist
        self.mask = mask
        # This is the maximum number of predictions to parse/encode/predict at a time
        self.prediction_batch_size = prediction_batch_size
        # This is the size of the batch tensorflow will use to make the predictions
        self.tensorflow_batch_size = tensorflow_batch_size

        # Batch vars
        self.batches = {}
        self.prepared_vcf_records = []

        # Counts
        self.batch_predictions = 0
        self.total_predictions = 0
        self.total_vcf_records = 0

    def _clear_batch(self):
        self.batch_predictions = 0
        self.batches.clear()
        del self.prepared_vcf_records[:]

    def _process_batch(self):
        start = time.time()
        total_batch_predictions = 0
        logger.debug('Starting process_batch')

        # Sanity check dump of batch sizes
        batch_sizes = ["{}:{}".format(tensor_size, len(batch)) for tensor_size, batch in self.batches.items()]
        logger.debug('Batch Sizes: {}'.format(batch_sizes))

        # Collect each batch's predictions
        batch_preds = {}
        for tensor_size, batch in self.batches.items():
            # Convert list of encodings into a proper sized numpy matrix
            prediction_batch = np.concatenate(batch, axis=0)

            # Run predictions
            batch_preds[tensor_size] = np.mean(
                get_preds(self.ann, prediction_batch, self.prediction_batch_size), axis=0
            )

        # Iterate over original list of vcf records, reconstructing record with annotations
        for prepared_record in self.prepared_vcf_records:
            record_predictions = self._write_record(prepared_record, batch_preds)
            total_batch_predictions += record_predictions

        self._clear_batch()
        logger.debug('Predictions: {}, VCF Records: {}'.format(self.total_predictions, self.total_vcf_records))
        duration = time.time() - start
        preds_per_sec = total_batch_predictions / duration
        preds_per_hour = preds_per_sec * 60 * 60
        logger.debug('Finished in {:0.2f}s, per sec: {:0.2f}, per hour: {:0.2f}'.format(duration,
                                                                                        preds_per_sec,
                                                                                        preds_per_hour))

    def _write_record(self, prepared_record, batch_preds):
        record = prepared_record.vcf_record
        gene_info = prepared_record.gene_info
        record_predictions = 0

        all_y_ref = []
        all_y_alt = []

        # Each prediction in the batch is located and put into the correct y
        for location in prepared_record.locations:
            # No prediction here
            if location.tensor_size == 0:
                if location.sequence_type == SequenceType_REF:
                    all_y_ref.append(None)
                else:
                    all_y_alt.append(None)
                continue

            # Extract the prediction from the batch into a list of predictions for this record
            batch = batch_preds[location.tensor_size]
            if location.sequence_type == SequenceType_REF:
                all_y_ref.append(batch[[location.batch_index], :, :])
            else:
                all_y_alt.append(batch[[location.batch_index], :, :])

        delta_scores = extract_delta_scores(
            all_y_ref=all_y_ref,
            all_y_alt=all_y_alt,
            record=record,
            ann=self.ann,
            dist_var=self.dist,
            mask=self.mask,
            gene_info=gene_info,
        )

        # If there are predictions, write them to the VCF INFO section
        if len(delta_scores) > 0:
            record.info['SpliceAI'] = delta_scores
            record_predictions += len(delta_scores)

        self.output.write(record)
        return record_predictions

    def add_record(self, record):
        """
        Adds a record to a batch. It'll capture the gene information for the record and
        save it for later to avoid looking it up again, then it'll encode ref and alt from
        the VCF record and place the encoded values into lists of matching sizes. Once the
        encoded values are added, a BatchLookupIndex is created so that after the predictions
        are made, it knows where to look up the corresponding prediction for the vcf record.

        Once the batch size hits it's capacity, it'll process all the predictions for the
        encoded batches.
        """

        self.total_vcf_records += 1
        # Collect gene information for this record
        gene_info = self.ann.get_name_and_strand(record.chrom, record.pos)

        # Keep track of how many predictions we're going to make
        prediction_count = len(record.alts) * len(gene_info.genes)
        self.batch_predictions += prediction_count
        self.total_predictions += prediction_count

        # Collect lists of encoded ref/alt sequences
        x_ref, x_alt = encode_batch_records(record, self.ann, self.dist, gene_info)

        # List of BatchLookupIndex's so we know how to lookup predictions for records from
        # the batches
        batch_lookup_indexes = []

        # Process the encodings into batches
        for var_type, encoded_seq in zip((SequenceType_REF, SequenceType_ALT), (x_ref, x_alt)):

            if len(encoded_seq) == 0:
                # Add BatchLookupIndex with zeros so when the batch collects the outputs
                # it knows that there is no prediction for this record
                batch_lookup_indexes.append(BatchLookupIndex(var_type, 0, 0))
                continue

            # Iterate over the encoded sequence and drop into the correct batch by size and
            # create an index to use to pull out the result after batch is processed
            for row in encoded_seq:
                # Extract the size of the sequence that was encoded to build a batch from
                tensor_size = row.shape[1]

                # Create batch for this size
                if tensor_size not in self.batches:
                    self.batches[tensor_size] = []

                # Add encoded record to batch
                self.batches[tensor_size].append(row)

                # Get the index of the record we just added in the batch
                cur_batch_record_ix = len(self.batches[tensor_size]) - 1

                # Store a reference so we can pull out the prediction for this item from the batches
                batch_lookup_indexes.append(BatchLookupIndex(var_type, tensor_size, cur_batch_record_ix))

        # Save the batch locations for this record on the composite object
        prepared_record = PreparedVCFRecord(
            vcf_record=record, gene_info=gene_info, locations=batch_lookup_indexes
        )
        self.prepared_vcf_records.append(prepared_record)

        # If we're reached our threshold for the max items to process, then process the batch
        if self.batch_predictions >= self.prediction_batch_size:
            self._process_batch()

    def finish(self):
        """
        Method to process all the remaining items that have been added to the batch.
        """
        if len(self.prepared_vcf_records) > 0:
            self._process_batch()
