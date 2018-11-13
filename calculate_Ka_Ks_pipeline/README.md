# A pipeline used to compute Ka and Ks

This pipeline adopted a method named [gamma-MYN method](https://biologydirect.biomedcentral.com/articles/10.1186/1745-6150-4-20) (a Modified version of Yang-Nielsen method) to estimate Ka and Ks values. [MAFFT (L-INS-i)](https://mafft.cbrc.jp/alignment/software/) is used to perform pairwise alignment of protein sequences for each duplicate gene pair.

| | |
| --- | --- |
| Authors | Xin Qiao ([qiaoxin](https://github.com/qiao-xin)) |
| | Leiting Li ([lileiting](https://github.com/lileiting)) |
| Email   | <qiaoxinqx2011@126.com> |

## Dependencies

- [Perl](https://www.perl.org)
- [BioPerl](https://bioperl.org)
- [PAL2NAL (v14)](http://www.bork.embl.de/pal2nal/#Download)
- [MAFFT (v7.402)](https://mafft.cbrc.jp/alignment/software/)
- [KaKs_Calculator 2.0](https://sourceforge.net/projects/kakscalculator2/)

## Installation

```bash
git clone https://github.com/qiao-xin/Scripts_for_GB.git
```

## Running
Once the required dependencies have been installed, try running this pipeline on the example data:
```bash
perl computing_Ka_Ks_pipe.pl -d data/Ath.cds.example -g data/Ath.tandem.pairs.example -o Ath.tandem.pairs.axt
```
**Note:** The CDS sequences (FASTA format) can be downloaded from [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html), NCBI, Ensembl Plants, etc. The different modes of duplicated gene pairs from different plants are available on [PlantDGD]().

## Results Files
### 1 - Ath.tandem.pairs.axt
The aligned pairwise sequences with AXT format was used as input file for KaKs_Calculator.

### 2 - Ath.tandem.pairs.axt.kaks
KaKs_Calculator generates this file that contains Ka, Ks, Ka/Ks values and other informations.

**Note:** KaKs_Calculator provides comprehensive information estimated from compared sequences, including numbers of synonymous and nonsynonymous sites, numbers of synonymous and nonsynonymous substitutions, GC contents, maximum-likelihood score, and AICC, in addition to synonymous and nonsynonymous substitution rates and their ratio. Meanwhile, Fisher’s exact test for small sample is applied to justify the validity of Ka and Ks calculated by these methods.

### 3 - Ath.tandem.pairs.axt.kaks.simplified
This is a simplified version of KaKs_Calculator output file, which only contains Ka, Ks, Ka/Ks and P-value.
```
Duplicate 1	Duplicate 2	Ka	Ks	Ka/Ks	P-Value
AT1G01580.1	AT1G01590.1	0.230604	1.32053	0.17463	2.55813e-42
AT1G01660.1	AT1G01670.1	0.518339	2.09868	0.246984	8.39747e-21
AT1G01670.1	AT1G01680.1	0.365625	1.26443	0.289163	7.50956e-11
AT1G02190.1	AT1G02205.3	0.310436	3.41288	0.0909601	1.52263e-75
```

### 4 - Ath.tandem.pairs.axt.KaKs_Calculator.log
KaKs_Calculator generates this log file:
```
Method(s): GMYN 
Genetic code: 1-Standard Code
Please wait while reading sequences and calculating...
1 AT1G01580.1-AT1G01590.1	[OK]
2 AT1G01660.1-AT1G01670.1	[OK]
3 AT1G01670.1-AT1G01680.1	[OK]
4 AT1G02190.1-AT1G02205.3	[OK]
5 AT1G02220.1-AT1G02230.1	[OK]
6 AT1G02230.1-AT1G02250.1	[OK]
7 AT1G02300.1-AT1G02305.1	[OK]
8 AT1G02430.1-AT1G02440.1	[OK]
9 AT1G02470.2-AT1G02475.1	[OK]
10 AT1G02520.1-AT1G02530.1	[OK]
Mission accomplished. (Time elapsed: 0:0)
```

## Citation
The manuscript is under review.
