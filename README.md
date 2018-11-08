# SpliceAI: A deep learning-based tool to evaluate the splice-altering potential of variants

## Table of contents

  * [Installation](#installation)
  * [Usage](#usage)
  * [Examples](#examples)
  * [Contact](#contact)

## Installation

The simplest way to install SpliceAI is through pip:
```
pip install spliceai
```

Alternately, SpliceAI can also be installed using setup.py via the following set of commands:
```
git clone https://github.com/Illumina/SpliceAI.git
cd SpliceAI
python setup.py install
```

## Usage

SpliceAI can be run from the command line: 
```
spliceai -I input.vcf -O output.vcf -R genome.fa [-A annotations.tsv]
```
| Argument | Description |
| -------- | ----------- |
|    -I    | The input VCF file containing the variants of interest. |
|    -O    | The output VCF file containing the SpliceAI predictions. The predictions are appended to the INFO column of the input file (`ALT|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`: see below for additional details). Currently, only single nucleotide variants and simple indels (either the ref or alt should be of length 1) are supported. If a variant belongs to multiple genes, predictions for each gene are provided. Variants which do not belong to any gene are not considered for scoring. |
|    -R    | Reference genome fasta file (assembly should be hg19/GRCh37 if -A argument is set to default). |
|    -A    | A tsv file where the columns correspond to gene name, chromosome, strand, transcription start, and transcription end. This argument is optional, and the default is the hg19/GRCH37 GENCODE.v24lift37 annotations file found at `spliceai/annotations/GENCODE.v24lift37` (use it as a template if providing a custom annotations file). |

For the sake of convenience, we have already calculated the predictions for all possible single nucleotide variants within the genic regions (3.4 billion variants). The results are available [here](https://basespace.illumina.com/s/5u6ThOblecrh).

## Examples

An example input file and the corresponding output file can be found at `examples/input.vcf` and `examples/output.vcf` respectively.

## Contact

Kishore Jaganathan: kishorejaganathan@gmail.com


