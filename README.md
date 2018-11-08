# SpliceAI: A deep learning-based tool to evaluate the splice-altering potential of variants

## Table of contents

  * [Installation](#installation)
  * [Usage](#usage)
  * [Examples](#examples)
  * [FAQ](#faq)
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
|    -O    | The output VCF file containing the SpliceAI predictions. The predictions are appended to the INFO column of the input file. Currently, only single nucleotide variants and simple indels (either the ref or alt should be of length 1) are supported. |
|    -R    | Reference genome fasta file (assembly should be hg19/GRCh37 if -A argument is set to default). |
|    -A    | A tsv file where the columns correspond to gene name, chromosome, strand, transcription start, transcription end, exon starts and exon ends (default is the hg19/GRCH37 GENCODE.v24lift37 annotations file found at `spliceai/annotations/GENCODE.v24lift37`). |

## Examples

An example input file and the corresponding output file can be found at `examples/input.vcf` and `examples/output.vcf` respectively.

## FAQ



## Contact

Kishore Jaganathan: kishorejaganathan@gmail.com


