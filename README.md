# SpliceAI: A deep learning-based tool to identify splice variants
[![Downloads](https://pepy.tech/badge/spliceai)](https://pepy.tech/project/spliceai) (Downloads since Nov 8th, 2018, via pypi)

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
|    -O    | The output VCF file containing the SpliceAI predictions. The predictions are appended to the INFO column of the input file (`SpliceAI=ALT\|SYMBOL\|DS_AG\|DS_AL\|DS_DG\|DS_DL\|DP_AG\|DP_AL\|DP_DG\|DP_DL`: see the table below for additional details). Currently, only single nucleotide variants and simple indels (either the ref or alt length is 1) are supported. If a variant belongs to multiple genes, then predictions for each gene are provided separately. Variants which do not belong to any gene are not evaluated. |
|    -R    | Reference genome fasta file (should be on hg19/GRCh37 if using default -A parameter). |
|    -A    | A tsv file where the columns correspond to gene name, chromosome, strand, transcription start, and transcription end. This argument is optional, and the default is the hg19/GRCh37 GENCODE.v24lift37 annotations file found at `spliceai/annotations/GENCODE.v24lift37` (use it as a template if providing a custom annotations file). |

|    ID    | Description |
| -------- | ----------- |
|  ALT     | Alternate nucleotide |
|  SYMBOL  | HGNC gene symbol |
|  DS_AG   | Delta score (acceptor gain) |
|  DS_AL   | Delta score (acceptor loss) |
|  DS_DG   | Delta score (donor gain) |
|  DS_DL   | Delta score (donor loss) |
|  DP_AG   | Delta position (acceptor gain) relative to the variant position |
|  DP_AL   | Delta position (acceptor loss) relative to the variant position |
|  DP_DG   | Delta position (donor gain) relative to the variant position |
|  DP_DL   | Delta position (donor loss) relative to the variant position |

Delta scores can have values ranging from 0 to 1, and correspond to the effect size of the variant with respect to altering splicing (values are calculated separately for acceptor gain, acceptor loss, donor gain and donor loss). In the paper, a detailed characterization is provided for threshold values 0.2 (high recall), 0.5, and 0.8 (high precision). Delta positions convey information about the location where splicing changes relative to the variant position (positive/negative values imply the location is to the right/left of the variant respectively). 

For the sake of convenience, we have already calculated the predictions for all possible single nucleotide variants within the genic regions (3.4 billion variants). The results are available [here](https://basespace.illumina.com/s/5u6ThOblecrh).

## Examples

An example input file and the corresponding output file can be found at `examples/input.vcf` and `examples/output.vcf` respectively.

## Contact

Kishore Jaganathan: kishorejaganathan@gmail.com


