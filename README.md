## SpliceAI: A deep learning-based tool to identify splice variants
[![release](https://img.shields.io/badge/release-v1.3.1-orange.svg)](https://img.shields.io/badge/release-v1.3.1-orange.svg)
[![license](https://img.shields.io/badge/license-GPLv3-green.svg)](https://img.shields.io/badge/license-GPLv3-green.svg)
[![downloads](https://pepy.tech/badge/spliceai)](https://pepy.tech/badge/spliceai)

This package annotates genetic variants with their predicted effect on splicing, as described in [Jaganathan *et al*, Cell 2019 in press](https://doi.org/10.1016/j.cell.2018.12.015).

**Update**: The annotations for all possible substitutions, 1 base insertions, and 1-4 base deletions within genes are available [here](https://basespace.illumina.com/s/otSPW8hnhaZR) for download.

### Installation
The simplest way to install SpliceAI is through pip or conda:
```sh
pip install spliceai
# or
conda install -c bioconda spliceai
```

Alternately, SpliceAI can be installed from the [github repository](https://github.com/Illumina/SpliceAI.git):
```sh
git clone https://github.com/Illumina/SpliceAI.git
cd SpliceAI
python setup.py install
```

SpliceAI requires ```tensorflow>=1.2.0```, which is best installed separately via pip or conda (see the [TensorFlow](https://www.tensorflow.org/) website for other installation options):
```sh
pip install tensorflow
# or
conda install tensorflow
```

### Usage
SpliceAI can be run from the command line:
```sh
spliceai -I input.vcf -O output.vcf -R genome.fa -A grch37
# or you can pipe the input and output VCFs
cat input.vcf | spliceai -R genome.fa -A grch37 > output.vcf
```

Required parameters:
 - ```-I```: Input VCF with variants of interest.
 - ```-O```: Output VCF with SpliceAI predictions `ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL` included in the INFO column (see table below for details). Only SNVs and simple INDELs (REF or ALT is a single base) within genes are annotated. Variants in multiple genes have separate predictions for each gene.
 - ```-R```: Reference genome fasta file. Can be downloaded from [GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz) or [GRCh38/hg38](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz).
 - ```-A```: Gene annotation file. Can instead provide `grch37` or `grch38` to use GENCODE V24 canonical annotation files included with the package. To create custom gene annotation files, use `spliceai/annotations/grch37.txt` in repository as template.

Optional parameters:
 - ```-D```: Maximum distance between the variant and gained/lost splice site (default: 50).
 - ```-M```: Mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss (default: 0).

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

Delta score of a variant, defined as the maximum of (DS_AG, DS_AL, DS_DG, DS_DL), ranges from 0 to 1 and can be interpreted as the probability of the variant being splice-altering. In the paper, a detailed characterization is provided for 0.2 (high recall), 0.5 (recommended), and 0.8 (high precision) cutoffs. Delta position conveys information about the location where splicing changes relative to the variant position (positive values are downstream of the variant, negative values are upstream).

### Examples
A sample input file and the corresponding output file can be found at `examples/input.vcf` and `examples/output.vcf` respectively. The output `T|RYR1|0.00|0.00|0.91|0.08|-28|-46|-2|-31` for the variant `19:38958362 C>T` can be interpreted as follows:
* The probability that the position 19:38958360 (=38958362-2) is used as a splice donor increases by 0.91.
* The probability that the position 19:38958331 (=38958362-31) is used as a splice donor decreases by 0.08.

Similarly, the output `CA|TTN|0.07|1.00|0.00|0.00|-7|-1|35|-29` for the variant `2:179415988 C>CA` has the following interpretation:
* The probability that the position 2:179415981 (=179415988-7) is used as a splice acceptor increases by 0.07.
* The probability that the position 2:179415987 (=179415988-1) is used as a splice acceptor decreases by 1.00.

### Frequently asked questions

**1. Why are some variants not scored by SpliceAI?**

SpliceAI only annotates variants within genes defined by the gene annotation file. Additionally, SpliceAI does not annotate variants if they are close to chromosome ends (5kb on either side), deletions of length greater than twice the input parameter ```-D```, or inconsistent with the reference fasta file.

**2. What are the differences between raw (```-M 0```) and masked (```-M 1```) precomputed files?**

The raw files also include splicing changes corresponding to strengthening annotated splice sites and weakening unannotated splice sites, which are typically much less pathogenic than weakening annotated splice sites and strengthening unannotated splice sites. The delta scores of such splicing changes are set to 0 in the masked files. We recommend using raw files for alternative splicing analysis and masked files for variant interpretation.

**3. Can SpliceAI be used to score custom sequences?**

Yes, install SpliceAI and use the following script:  

```python
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
import numpy as np

input_sequence = 'CGATCTGACGTGGGTGTCATCGCATTATCGATATTGCAT'
# Replace this with your custom sequence

context = 10000
paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
models = [load_model(resource_filename('spliceai', x)) for x in paths]
x = one_hot_encode('N'*(context//2) + input_sequence + 'N'*(context//2))[None, :]
y = np.mean([models[m].predict(x) for m in range(5)], axis=0)

acceptor_prob = y[0, :, 1]
donor_prob = y[0, :, 2]
```

### Contact
Kishore Jaganathan: kishorejaganathan@gmail.com
