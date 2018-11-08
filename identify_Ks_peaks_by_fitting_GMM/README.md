# A pipeline used to identify Ks peaks by fitting the Gaussian Mixture Model(GMM)
This pipline aimed to calculate the Ks peaks from inter or intra-species by using Gaussian Mixture Model(GMM) and provide visualization for Ks peaks. It needs some output files from [MCScanX](http://chibba.pgml.uga.edu/mcscan2/) and [calculate_Ka_Ks_pipeline](https://github.com/qiao-xin/Scripts_for_GB/tree/master/calculate_Ka_Ks_pipeline).

| | |
| --- | --- |
| Authors | Xin Qiao ([qiaoxin](https://github.com/qiao-xin)) |
| Email   | <qiaoxinqx2011@126.com> |

## Dependencies

- [Perl](https://www.perl.org)
- [BioPerl](https://bioperl.org)
- [Python](https://www.python.org/)
- [pandas](http://pandas.pydata.org/)
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
**Note:** The output file will be automatically named to Ath.synteny.blocks.ks.info and Ks value will be contained in it.

### 4. Calculating Ks peaks by fitting the GMM

~~~bash
python plot_syntenic_blocks_ks_distri.py Ath.synteny.blocks.ks.info 2 Ath
~~~
**Note:** Three parameters should be provided:file_name, component, and abbreviation of species. The file_name is name of files obtained from last step

