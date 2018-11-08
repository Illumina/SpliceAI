## SpliceAI: A deep learning-based tool to identify splice variants

### Installation

The simplest way to install SpliceAI is through pip:
```
pip install spliceai
```

Alternately, SpliceAI can also be installed using setup.py via the following commands:
```
git clone https://github.com/Illumina/SpliceAI.git
cd SpliceAI
python setup.py install
```

### Running SpliceAI

SpliceAI can be run from the command line 
```
spliceai -I input.vcf -O output.vcf -R genome.fa
```
input.vcf: The input file (VCF format) containing the list of variants to be scored.
