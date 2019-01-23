# A pipeline used to compute Pearson's correlation coefficient between expression profiles of duplicate gene pairs  

This pipeline can be used to compute Pearson's coorelation coefficient (r) between expression profiles of gene pairs derived from different modes of gene duplication. Meanwhile, this pipeline can automatically generate 10,000 random gene pairs and compute r values for them.

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

## Running
The example transcriptome data (```example_data```) from different tissues of  *Arabidopsis thaliana*  was reported in a previous study [Evolutionary Fates and Dynamic Functionalization of Young Duplicate Genes in Arabidopsis Genomes](http://www.plantphysiol.org/content/172/1/427.abstract). 

|||
| --- | --- |
| output1 | Leaf |
| output2 | Flower |
| output3 | Stem |
| output4 | Root |
| output5 | Silique |

### 1. Input Files Preparation
Before running this pipline, some files should be prepared:

- [Ath.tandem.pairs](https://github.com/qiao-xin/Scripts_for_GB/blob/master/compute_pearson_r/Ath.tandem.pairs): Duplicate gene pairs in *Arabidopsis thaliana* 
- [Kallisto results](https://github.com/qiao-xin/Scripts_for_GB/tree/master/compute_pearson_r/example_data): [Kallisto](http://pachterlab.github.io/kallisto/about.html) output files(see [Example Data](https://github.com/qiao-xin/Scripts_for_GB/tree/master/compute_pearson_r/example_data))


### 2.Running
When you prepared Input files, this pipeline can provide a easy way to calculate PCC of gene pairs.

A typical command could look like this:
```coding
perl compute_pearsonr_pipeline.pl
```
This command will produce:

1.Alltpm: our pipeline can merge multiple expression data sets from Kallisto into one file

2.alltpm.log10: Log10 transformation for Alltpm values

3.Ath.tandem.pairs.pr：PCC values of gene pairs

4.random.pairs:1000 gene pairs randomly generated from *A.thaliana* genome

5.random.pairs.pr: PCC values of 1000 randomly generate gene pairs
