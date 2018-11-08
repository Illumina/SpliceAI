# SpliceAI: A deep learning-based tool to identify splice variants

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
spliceai -I input.vcf -O output.vcf -R genome.fa
```
| Argument | Description |
| -------- | ----------- |
|    -I    | The input file (VCF format) containing the list of variants to be scored |
|    -O    | The output file (VCF format), the SpliceAI predictions are added to the INFO column |
|    -R    | Reference genome fasta file |

## Examples

An example input file and the corresponding output file can be found at `examples/input.vcf` and `examples/output.vcf` respectively.

## Contact

Kishore Jaganathan: kishorejaganathan@gmail.com


