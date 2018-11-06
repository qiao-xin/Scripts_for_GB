# A pipeline used to compute Ka and Ks.

A method named [gamma-MYN method](https://biologydirect.biomedcentral.com/articles/10.1186/1745-6150-4-20) (a Modified version of Yang-Nielsen method) was used to estimate Ka and Ks values. [MAFFT (L-INS-i)](https://mafft.cbrc.jp/alignment/software/) was used to perform pairwise alignment of protein sequences for each duplicate gene pair.

| | |
| --- | --- |
| Authors | Xin Qiao ([qiaoxin](https://github.com/qiao-xin)) |
| | Leiting Li ([lileiting](https://github.com/lileiting)) |
| Email   | <qiaoxinqx2011@126.com> |

## Dependencies

- [Perl](https://www.perl.org)
- [PAL2NAL (v14)](http://www.bork.embl.de/pal2nal/#Download)
- [MAFFT (v7.402)](https://mafft.cbrc.jp/alignment/software/)
- [KaKs_Calculator 2.0](https://sourceforge.net/projects/kakscalculator2/)

## Installation

```bash
git clone https://github.com/qiao-xin/Scripts_for_GB.git
```

## Running

```bash
perl fa_get_set_genes_from_file.pl -d data/Ath.cds -g data/Ath.tandem.pairs -o Ath.tandem.pairs.axt
```

