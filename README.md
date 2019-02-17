## SpliceAI: A deep learning-based tool to identify splice variants
This package annotates genetic variants with their predicted effect on splicing, as described in [Jaganathan *et al*, Cell 2019 in press](https://doi.org/10.1016/j.cell.2018.12.015).

### Installation
The simplest way to install SpliceAI is through pip:
```sh
pip install spliceai
```

Alternately, SpliceAI can be installed from the [github repository](https://github.com/Illumina/SpliceAI.git):
```sh
git clone https://github.com/Illumina/SpliceAI.git
cd SpliceAI
python setup.py install
```

SpliceAI requires [tensorflow](https://www.tensorflow.org/install/)>=1.2.0, which is best installed separately via pip: `pip install tensorflow`. See the TensorFlow website for other installation options.

### Usage
SpliceAI can be run from the command line:
```sh
spliceai -I input.vcf -O output.vcf -R genome.fa -A annotations.txt

# or you can pipe the input and output VCFs
cat input.vcf | spliceai -R genome.fa -A annotations.txt > output.vcf
```

Options:
 - **-I**: Input VCF with variants of interest.
 - **-O**: Output VCF with SpliceAI predictions `SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL` included in the INFO column (see table below for details). Only SNVs and simple INDELs (ref or alt must be a single base) within genes are annotated. Variants in multiple genes have separate predictions for each gene.
 - **-R**: Reference genome fasta file.
 - **-A**: Gene annotation file. Can instead provide `grch37` or `grch38` to use GENCODE canonical annotation files included with the package. To create custom annotation files, use `spliceai/annotations/grch37.txt` in repository as template.

**Note**: The annotations for all possible SNVs within genes are available [here](https://basespace.illumina.com/s/5u6ThOblecrh) for download.

Details of SpliceAI INFO field:

|    ID    | Description |
| -------- | ----------- |
|  ALLELE  | Alternate allele |
|  SYMBOL  | Gene symbol |
|  DS_AG   | Delta score (acceptor gain) |
|  DS_AL   | Delta score (acceptor loss) |
|  DS_DG   | Delta score (donor gain) |
|  DS_DL   | Delta score (donor loss) |
|  DP_AG   | Delta position (acceptor gain) |
|  DP_AL   | Delta position (acceptor loss) |
|  DP_DG   | Delta position (donor gain) |
|  DP_DL   | Delta position (donor loss) |

**Delta score** of a variant ranges from 0 to 1, and can be interpreted as the probability of the variant being splice-altering. In the paper, a detailed characterization is provided for 0.2 (high recall/likely pathogenic), 0.5 (recommended/pathogenic), and 0.8 (high precision/pathogenic) cutoffs. **Delta position** conveys information about the location where splicing changes relative to the variant position (positive values are upstream of the variant, negative values are downstream).

### Examples
A sample input file and the corresponding output file can be found at `examples/input.vcf` and `examples/output.vcf` respectively (`grch37` annotation). The output `SpliceAI=T|RYR1|0.22|0.00|0.91|0.70|-107|-46|-2|90` for the variant `19:38958362 C>T` can be interpreted as follows:
* The probability that the position `19:38958255` is used as a splice acceptor increases by `0.22`.
* The probability that the position `19:38958360` is used as a splice donor increases by `0.91`.
* The probability that the position `19:38958452` is used as a splice donor decreases by `0.70`.

### Contact
Kishore Jaganathan: kishorejaganathan@gmail.com
