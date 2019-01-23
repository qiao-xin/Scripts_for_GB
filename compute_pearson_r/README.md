# A pipeline used to compute Pearson's correlation coefficient  

This pipeline can be used to compute Pearson's coorelation coefficient (*r*) between expression profiles of duplicate gene pairs derived from different modes of gene duplication. Meanwhile, this pipeline can automatically generate 10,000 random gene pairs and compute *r* values for them.

| | |
| --- | --- |
| Authors | Xin Qiao ([qiaoxin](https://github.com/qiao-xin)) |
| | Qionghou Li ([Qionghou Li](https://github.com/LQHHHHH)) |
| Email   | <qiaoxinqx2011@126.com> |

## Dependencies
- [Perl](https://www.perl.org/)
- [Python](https://www.python.org/)
  - [Pandas](http://pandas.pydata.org/)
  - [Scipy](https://www.scipy.org/)
 
 ## Installation

```bash
git clone https://github.com/qiao-xin/Scripts_for_GB.git
```

## Preparing input files

Before running this pipline, some data files should be prepared:

### 1. gene expression profiles

For example, we peformed [Kallisto](http://pachterlab.github.io/kallisto/about.html) to quantify expression abundance of the transcripts in different tissues of *Arabidopsis thaliana* using RNA-seq reads generated in a previous study [(Wang et al. 2016)](http://www.plantphysiol.org/content/172/1/427.abstract). 

The result files were provided in ```example_data```:

- example_data
  - output1
  - output2
  - output3
  - output4
  - output5

|||
| --- | --- |
| output1 | Leaf |
| output2 | Flower |
| output3 | Stem |
| output4 | Root |
| output5 | Silique |

### 2. duplicate gene pairs

The duplicate gene pairs identified from different plants are available on [PlantDGD]().

Here, we only estimated ***r*** values between expression profiles of tandem gene pairs in *A. thaliana*. Also, you can include other types of duplicate gene pairs (e.g. WGD, PD, TRD) in current directory.


## Running

Once the required files have been prepared, try running this pipeline on the example data.

USAGE:
```coding
perl compute_pearsonr_pipeline.pl
```

## Results Files


### 1 - alltpm
The integrated spatial expression for each gene in different tissues. The data in ```alltpm``` looks like this (tab separated):

```

```
This pipeline can extract expression values (measured as TPM) from different tissues ([Kallisto outputs](https://github.com/qiao-xin/Scripts_for_GB/tree/master/compute_pearson_r/example_data)) for each gene, and then merge into one file.

### 2 - alltpm.log10
Log10-transformed TPM values for each gene.

### 3 - Ath.tandem.pairs.pr
Pearson's ***r*** values for tandem gene pairs.

### 4 - random.pairs
10,000 random gene pairs were generated for *A.thaliana* when running this pipeline.

### 5 - random.pairs.pr
Pearson's ***r*** values for 10,000 random gene pairs.

## Citation
*Qiao X, Li Q, Yin H, Qi K, Li L, Wang R, Zhang S\* and Paterson AH\*: Gene duplication and evolution in recurring polyploidization-diploidization cycles in plants. Under Review.*
