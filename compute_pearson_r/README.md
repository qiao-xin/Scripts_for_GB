# A pipeline used to compute Pearson's correlation coefficient  

This pipeline can be used to compute Pearson's coorelation coefficient (*r*) between expression profiles of gene pairs derived from different modes of gene duplication. Meanwhile, this pipeline can automatically generate 10,000 random gene pairs and compute *r* values for them.

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

### 1. 

- [Kallisto results](https://github.com/qiao-xin/Scripts_for_GB/tree/master/compute_pearson_r/example_data): [Kallisto](http://pachterlab.github.io/kallisto/about.html) output files(see [example Data](https://github.com/qiao-xin/Scripts_for_GB/tree/master/compute_pearson_r/example_data))

The example transcriptome data (```example_data```) from different tissues of  *Arabidopsis thaliana*  was reported in a previous study [(Wang et al. 2016)](http://www.plantphysiol.org/content/172/1/427.abstract). 

|||
| --- | --- |
| output1 | Leaf |
| output2 | Flower |
| output3 | Stem |
| output4 | Root |
| output5 | Silique |

### 2. Duplicate gene pairs

- [Ath.tandem.pairs](https://github.com/qiao-xin/Scripts_for_GB/blob/master/compute_pearson_r/Ath.tandem.pairs): Duplicate gene pairs in *Arabidopsis thaliana* 


## Running

USAGE:
Once the required files have been prepared, try running this pipeline on the example data:

```coding
perl compute_pearsonr_pipeline.pl
```

## Results Files
This command will produce:

### 1 - Alltpm
our pipeline can merge multiple expression data sets from Kallisto into one file

2.alltpm.log10: Log10 transformation for Alltpm values

3.Ath.tandem.pairs.pr：PCC values of gene pairs

4.random.pairs:1000 gene pairs randomly generated from *A.thaliana* genome

5.random.pairs.pr: PCC values of 1000 randomly generate gene pairs

## Citation
*Qiao X, Li Q, Yin H, Qi K, Li L, Wang R, Zhang S\* and Paterson AH\*: Gene duplication and evolution in recurring polyploidization-diploidization cycles in plants. Under Review.*
