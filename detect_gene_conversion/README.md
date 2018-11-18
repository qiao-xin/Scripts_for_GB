# A pipeline used to detect gene conversion events

This pipeline incorporated a sophisticated method that was firstly reported in a previous study ([Wang et al. 2007](http://www.genetics.org/content/177/3/1753)) to detect gene conversion events. This method was originally developed by Dr. Xiyin Wang, and was applied to subsequent a series of studies ([Wang et al. 2009](https://genome.cshlp.org/content/19/6/1026.full); [Wang et al. 2011](http://www.plantcell.org/content/23/1/27); [Liu et al. 2014](https://www.nature.com/articles/ncomms4930)). 

| | |
| --- | --- |
| Authors | Xin Qiao ([qiaoxin](https://github.com/qiao-xin)) |
| | Xiyin Wang ([PGML](http://www.plantgenome.uga.edu)) |
| | Andrew Paterson ([PGML](http://www.plantgenome.uga.edu)) |
| Email   | <qiaoxinqx2011@126.com> |

## Dependencies

- [Perl](https://www.perl.org)
- [BioPerl](https://bioperl.org)
- [Clustal W](http://www.clustal.org/clustal2/#Download)
- [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html)

## Installation

```bash
git clone https://github.com/qiao-xin/Scripts_for_GB.git
```

## Preparing input files
Here, we used this pipeline to detect gene conversion events in WGD-derived gene pairs for *Arabidopsis thaliana* (Ath) using *Arabidopsis lyrata* (Aly) as outgroup.

Before running this pipline, some data files should be prepared:
- ```Ath.wgd.pairs```: duplicate gene pairs derived from different modes of gene duplication or WGD events of different ages.
- ```Ath.cds``` and ```Aly.cds```: whole-genome CDS sequences (FASTA format) of investigated species (Ath) and outgroup species (Aly).
- ```Ath_Aly.collinearity```: the execution of [MCScanX](http://chibba.pgml.uga.edu/mcscan2/) will output this file, containing pairwise collinear blocks as follows:
```
## Alignment 0: score=64951.0 e_value=0 N=1326 Aly-scaffold_1&Ath-Chr1 plus
  0-  0:	Aly333558	AT1G02210.1	  3e-11
  0-  1:	Aly470187	AT1G02260.1	      0
  0-  2:	Aly470188	AT1G02270.1	      0
  0-  3:	Aly311338	AT1G02280.1	      0
  0-  4:	Aly470192	AT1G02310.1	      0
  0-  5:	Aly470194	AT1G02335.1	 4e-152
```

## Running

**USAGE:**
```bash
perl GeConScan.pl species_A_abbr species_B_abbr gene_order_in_collinearity paralog_pairs_in_A
```
The parameter ```species_A_abbr``` and ```species_B_abbr``` indicates scientific name abbreviation you defined for the investigated species and outgroup species respectively, for example, "Ath" and "Aly". You can set the option ```gene_order_in_collinearity``` by 'f' or 'r'. If the target genome gene_id located before the outgroup genome gene_id in file "xyz.collinearity", please use 'f'. ```paralog_pairs_in_A``` indicates the input file containing paralogous gene pairs from target species.

A typical command could look like this:
```bash
perl GeConScan.pl Ath Aly r Ath.wgd.pairs.example
```
This command will produce:
- Ath_Aly.homologous.quartets
- Ath_Aly.quartets.K.CV.BP.txt

