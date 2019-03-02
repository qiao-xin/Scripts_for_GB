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

The duplicate gene pairs identified from different plants are available on [PlantDGD](http://pdgd.njau.edu.cn:8080).

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
GeneID	1	2	3	4	5
AT1G09645.1	89.2424	58.1628	73.8257	79.6055	5.44928
AT5G02280.1	18.9473	24.7855	28.309	33.0487	0.923299
AT1G70700.1	477.756	183.034	181.575	20.1834	7.17899
AT1G25097.1	0.694941	0	0	0	0
AT2G30410.1	77.7176	95.5345	127.086	87.8365	4.49006
```
This pipeline can extract expression values (measured as TPM) from different tissues ([Kallisto outputs](https://github.com/qiao-xin/Scripts_for_GB/tree/master/compute_pearson_r/example_data)) for each gene, and then merge into one file.

### 2 - alltpm.log10
Log10-transformed TPM values for each gene.
```
GeneID	1	2	3	4	5
AT1G09645.1	1.95057124127048	1.76464530561024	1.86820757354698	1.90094307448557	0.736339123802573
AT5G02280.1	1.277547331526	1.39419768438907	1.45192452842012	1.51915438079804	-0.0346576348116481
AT1G70700.1	2.67920614994061	2.26253177082796	2.25905605283343	1.30499432725574	0.856063348381457
AT1G25097.1	-0.158052065139533	-4	-4	-4	-4
AT2G30410.1	1.89051938067131	1.98016023497532	2.10409771060411	1.94367502222247	0.652252144454283
```

### 3 - Ath.tandem.pairs.pr
Pearson's *r* values for tandem gene pairs.

### 4 - random.pairs
10,000 random gene pairs were generated for *A.thaliana* when running this pipeline.

### 5 - random.pairs.pr
Pearson's *r* values for 10,000 random gene pairs.

## Citation
*[Qiao X, Li Q, Yin H, Qi K, Li L, Wang R, Zhang S, Paterson AH: Gene duplication and evolution in recurring polyploidization–diploidization cycles in plants. Genome Biology 2019, 20:38.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1650-2)*
