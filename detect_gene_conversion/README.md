# The pipeline used to detect gene conversion events

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

## Prepeations

- Duplicate gene pairs derived from different modes of gene duplication or WGD event of different ages.
- Whole-genome CDS sequences (FASTA format) of investigated species and outgroup species.

