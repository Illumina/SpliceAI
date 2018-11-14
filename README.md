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

SpliceAI also requires tensorflow>=1.2.0, which in most cases can be installed through pip:
```
pip install tensorflow
```

Please see [here](https://www.tensorflow.org/install/) for other installation options if this installation method does not work.

## Usage

SpliceAI can be run from the command line: 
```
spliceai -I input.vcf -O output.vcf -R genome.fa [-A annotations.tsv]
```
| Argument | Description |
| -------- | ----------- |
|    -I    | The input VCF file containing the variants of interest. Currently, only single nucleotide variants and simple indels (either the ref or alt length is 1) are scored. |
|    -O    | The output VCF file containing the SpliceAI predictions. The predictions are appended to the INFO column (`SpliceAI=ALT\|SYMBOL\|DS_AG\|DS_AL\|DS_DG\|DS_DL\|DP_AG\|DP_AL\|DP_DG\|DP_DL`: see the table below for additional details). If a variant belongs to multiple genes, then predictions for each gene are provided separately. Variants which do not belong to any gene are not scored. |
|    -R    | Reference genome fasta file (should be on hg19/GRCh37 if using default -A parameter). |
|    -A    | A tsv file where the columns correspond to gene name, chromosome, strand, transcription start, and transcription end. This argument is optional, and the default is the hg19/GRCh37 GENCODE.v24lift37 gene annotations file found at `spliceai/annotations/GENCODE.v24lift37` (use it as a template if providing a custom annotations file). |

|    ID    | Description |
| -------- | ----------- |
|  ALT     | Alternate nucleotide |
|  SYMBOL  | Gene symbol |
|  DS_AG   | Delta score (acceptor gain) |
|  DS_AL   | Delta score (acceptor loss) |
|  DS_DG   | Delta score (donor gain) |
|  DS_DL   | Delta score (donor loss) |
|  DP_AG   | Delta position (acceptor gain) |
|  DP_AL   | Delta position (acceptor loss) |
|  DP_DG   | Delta position (donor gain) |
|  DP_DL   | Delta position (donor loss) |

**Delta score** of a variant ranges from 0 to 1, and can be interpreted as the probability of the variant being splice-altering. In the paper, a detailed characterization is provided for 0.2 (high recall/likely pathogenic), 0.5 (recommended/pathogenic), and 0.8 (high precision/pathogenic) cutoffs. **Delta position** conveys information about the location where splicing changes relative to the variant position (positive/negative values imply that the location is to the right/left of the variant respectively).

For the sake of convenience, we have already calculated the outputs for all possible single nucleotide variants within the genic regions (3.4 billion variants). The results are available [here](https://basespace.illumina.com/s/5u6ThOblecrh).

## Examples

A sample input file and the corresponding output file can be found at `examples/input.vcf` and `examples/output.vcf` respectively. The output `SpliceAI=T|RYR1|0.22|0.00|0.91|0.70|-107|-46|-2|90` for the variant `19:38958362 C>T` can be interpreted as follows:
* The probability that the position `19:38958255` is used as a splice acceptor increases by `0.22`.
* The probability that the position `19:38958360` is used as a splice donor increases by `0.91`.
* The probability that the position `19:38958452` is used as a splice donor decreases by `0.70`.

## Contact

Kishore Jaganathan: kishorejaganathan@gmail.com


