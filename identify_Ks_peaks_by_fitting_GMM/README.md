# A pipeline used to identify Ks peaks by fitting the Gaussian Mixture Model(GMM)
This pipline aimed to calculate the Ks peaks from inter or intra-species by using Gaussian Mixture Model(GMM) and provide visualization for Ks peaks. It needs some output files from [MCScanX](http://chibba.pgml.uga.edu/mcscan2/) and [calculate_Ka_Ks_pipeline](https://github.com/qiao-xin/Scripts_for_GB/tree/master/calculate_Ka_Ks_pipeline).

| | |
| --- | --- |
| Authors | Xin Qiao ([qiaoxin](https://github.com/qiao-xin)) |
| Email   | <qiaoxinqx2011@126.com> |
| Authors | Qionghou Li ([Qionghou Li](https://github.com/LQHHHHH)) |
| Email   | <liqionghou1996@gamil.com> |

## Dependencies

- [Perl](https://www.perl.org)
- [BioPerl](https://bioperl.org)
- [Python](https://www.python.org/)
  - [Pandas](http://pandas.pydata.org/)
  - [Scipy](https://www.scipy.org/)
  - [Sklearn](https://scikit-learn.org/stable/index.html)

## Installation

```bash
git clone https://github.com/qiao-xin/Scripts_for_GB.git
```

## Running

### 1. File Preparation

Before using this pipline, some files should be prepared:
- Ath.collinearity: an output file from [MCScanX](http://chibba.pgml.uga.edu/mcscan2/).
- Ath.segmental.kaks: an output file from [calculate_Ka_Ks_pipeline](https://github.com/qiao-xin/Scripts_for_GB/tree/master/calculate_Ka_Ks_pipeline).

### 2. Adding Ka, Ks and Ka/Ks into collinearity file 

This step will add the Ka, Ks from Ath.segmental.kaks to Ath.collinearity

~~~bash
perl add_ka_ks_colinearity_file-single.pl Ath.segmental.kaks Ath.collinearity
~~~
**Note:** The output file will be automatically named to Ath.collinearity.kaks

### 3. Extracting Ks value from syntenic blocks

~~~bash
perl extract_synteny_blocks_ks-v2-single.pl Ath.collinearity.kaks
~~~
**Note:** The output file will be automatically named to Ath.synteny.blocks.ks.info and Ks values of syntenic blocks will be contained in it.

### 4. Calculating Ks peaks by fitting the GMM

~~~bash
python plot_syntenic_blocks_ks_distri.py Ath.synteny.blocks.ks.info 2 Ath
~~~
**Note:** Three parameters should be provided:file_name, component, and abbreviation of species. The file_name is the name of files obtained from last step. The component is the number of the mixture components and abbreviation of species is the abbreviation you gave to your species.

## Results Files
### 1 - Ath.collinearity.kaks
This file is an output file from step 2. the Ka, Ks and Ka/Ks values of gene pairs extracted from Ath.segmental.kaks are appended to Ath.collinearity file.

```
############### Parameters ###############
# MATCH_SCORE: 50
# MATCH_SIZE: 5
# GAP_PENALTY: -1
# OVERLAP_WINDOW: 5
# E_VALUE: 1e-05
# MAX GAPS: 25
############### Statistics ###############
# Number of collinear genes: 7544, Percentage: 27.52
# Number of all genes: 27416
##########################################
## Alignment 0: score=8970.0 e_value=0 N=190 Ath-Chr1&Ath-Chr1 plus Ka Ks Ka/Ks
  0-  0:	AT1G17240.1	AT1G72300.1	      0	0.136575	0.717908	0.19024
  0-  1:	AT1G17290.1	AT1G72330.3	      0	0.0746912	0.749983	0.0995905
  0-  2:	AT1G17310.1	AT1G72350.1	  2e-49	0.349746	1.3647	0.25628
  0-  3:	AT1G17350.2	AT1G72420.1	 5e-132	0.121138	0.9046	0.133913
  0-  4:	AT1G17380.1	AT1G72450.1	  2e-76	0.278683	1.69693	0.164228
  0-  5:	AT1G17400.1	AT1G72490.1	 3e-101	0.176897	1.22428	0.14449
  0-  6:	AT1G17420.1	AT1G72520.1	      0	0.0914749	1.3165	0.0694835
  0-  7:	AT1G17430.1	AT1G72620.1	      0	0.125315	0.760592	0.16476
```

### 2 - Ath.synteny.blocks.ks.info
This file is an output file from step 3. In this file, the Ks value of every syntenic block are computed by using average Ks values of collinear genes contained in syntenic block. The basic information of syntenic block such as blocks ID, location and block size are also contained.
```
Blocks ID	Location	Block Size	Average Ks	e-value	Score	Orientation
Alignment166	Ath-Chr3&Ath-Chr5	7	1.42200142857143	4.6e-12	330.0	plus
Alignment68	Ath-Chr1&Ath-Chr4	8	3.6449275	6.9e-15	341.0	minus
Alignment140	Ath-Chr2&Ath-ChrM	13	0.329314675	2.6e-36	610.0	plus
Alignment169	Ath-Chr3&Ath-Chr5	40	1.267692325	4.6e-136	1867.0	minus
Alignment121	Ath-Chr2&Ath-Chr4	7	1.48141166666667	4.4e-11	305.0	minus
```

### 3 - Ath.synteny.blocks.ks.distri.pdf
![Ath.synteny.blocks.ks.distri](/identify_Ks_peaks_by_fitting_GMM/data/Ath.synteny.blocks.ks.distri.png)
*Figure 1: The distribution of Ks values of syntenic blocks within Arabidopsis genome.*

This file is generated in step 4. Ks peak corresponding to WGDs of different age were inferred by fitting Gaussian mixture models to Ks distribution. The Ks peak values will be outputed when running python script.

## Citation
In preparation...
