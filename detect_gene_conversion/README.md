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
1. Ath_Aly.homologous.quartets
```
Duplicate 1	Homolog 1	Duplicate 2	Homolog 2
Ath.AT1G17380.1	Aly.Aly312706	Ath.AT1G72450.1	Aly.Aly476379
Ath.AT1G17400.1	Aly.Aly471945	Ath.AT1G72490.1	Aly.Aly476382
Ath.AT1G17455.1	Aly.Aly471952	Ath.AT1G72630.1	Aly.Aly339481
Ath.AT1G17480.1	Aly.Aly471956	Ath.AT1G72670.1	Aly.Aly476404
Ath.AT1G17530.1	Aly.Aly471960	Ath.AT1G72750.1	Aly.Aly476412
Ath.AT1G17590.1	Aly.Aly471968	Ath.AT1G72830.2	Aly.Aly339498
```

2. Ath_Aly.quartets.K.CV.BP.txt
```
Homologous_quartets	Pa(Ath-1, Ath-2)	Ps(Ath-1, Ath-2)	Pa(Aly-1, Aly-2)	Ps(Aly-1, Aly-2)	Pa(Ath-1, Aly-1)	Ps(Ath-1, Aly-1)	Pa(Ath-2, Aly-2)	Ps(Ath-2, Aly-2)	Is_conversion_between_Ath-1_and_Ath-2	Is_conversion_between_Aly-1_and_Aly-2	Bootstrap_percentage_Ks_Ath-1_Ath-2	Bootstrap_percentage_Ka_Ath-1_Ath-2	Bootstrap_percentage_Ks_Aly-1_Aly-2	Bootstrap_percentage_Ka_Aly-1_Aly-2
Ath.AT1G17380.1_Ath.AT1G72450.1_Aly.Aly312706_Aly.Aly476379	0.2493	1.0510	0.2763	0.9870	0.0493	0.1356	0.0394	0.1723	N	N	
Ath.AT1G17400.1_Ath.AT1G72490.1_Aly.Aly471945_Aly.Aly476382	0.2136	1.0556	0.2100	1.1067	0.0122	0.1398	0.0327	0.2425	N	N	
Ath.AT1G17455.1_Ath.AT1G72630.1_Aly.Aly471952_Aly.Aly339481	0.0695	0.8170	0.0738	0.7803	0.0000	0.0449	0.0040	0.2151	N	N	
Ath.AT1G17480.1_Ath.AT1G72670.1_Aly.Aly471956_Aly.Aly476404	0.1680	1.0117	0.1634	0.9545	0.0200	0.1268	0.0112	0.1642	N	N	
Ath.AT1G17530.1_Ath.AT1G72750.1_Aly.Aly471960_Aly.Aly476412	0.0888	0.9736	0.0900	0.9477	0.0224	0.1566	0.0074	0.1205	N	N	
Ath.AT1G17590.1_Ath.AT1G72830.2_Aly.Aly471968_Aly.Aly339498	0.2331	0.5479	0.1989	0.5094	0.0859	0.1327	0.0455	0.1487	N	N
```
**Note:**
- Ath-1 and Ath-2 show two Ath paralogous genes, Aly-1 and Aly-2 show two Aly paralogous genes.
- Pa () shows nonsynonymous nucleotide substitution difference, Ps () shows synonymous nucleotide substitution difference.
- 'N' shows no gene converison found between the Ath or Aly paralogs.
- 'R' shows that the two Ath paralogs have been affected by gene covnersion.
- 'S' shows that the two Aly paralogs have been affected by gene covnersion.

## Citation
The manuscript is under review.
