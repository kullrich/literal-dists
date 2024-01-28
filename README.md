[![License: GPLv3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language: C99](https://img.shields.io/badge/Language-C99-orangered.svg)](https://en.wikipedia.org/wiki/C99)

# literal-dists

Convert a FASTA alignment with ambiguous nucleotides to literal distance matrix

Accepts `A, C, G, T, R, Y, S, W, K, M, B, D, H, V, -, X, ., N`

Forked and adapted from [snp-dists](https://github.com/tseemann/snp-dists)

Almost all source code is taken from [snp-dists](https://github.com/tseemann/snp-dists),
despite the literal-distance calculation

Literal distances are implemented as mentioned by [Chang et al., 2017](https://pubmed.ncbi.nlm.nih.gov/28819774/)

## Quick Start

```
% cat test/good.aln

>seq1
AGTCAGTC
>seq2
AGGCAGTC
>seq3
AGTGAGTA
>seq4
TGTTAGAC
>seq5
TGTCAYGR

% literal-dists test/good.aln > distances.tab

Read 4 sequences of length 8

% cat distances.tab

literal-dists 0.0.1	seq1	seq2	seq3	seq4	seq5
seq1	0.000000	0.125000	0.250000	0.375000	0.500000
seq2	0.125000	0.000000	0.375000	0.500000	0.625000
seq3	0.250000	0.375000	0.000000	0.500000	0.562500
seq4	0.375000	0.500000	0.500000	0.000000	0.500000
seq5    0.500000	0.625000	0.562500	0.500000	0.000000
```

## Installation

As `literal-dists` source code is almost completely taken from [snp-dists](https://github.com/tseemann/snp-dists)

`literal-dists` is written in C to the C99 standard and only depends on `zlib`.

### Source

```
git clone https://github.com/kullrich/literal-dists.git
cd literal-dists
make

# run tests
make check

# optionally install to a specific location (default: /usr/local)
make PREFIX=/usr/local install
```

## Options

### `literal-dists -h` (help)

```
SYNOPSIS
  Pairwise literal distance matrix from a FASTA alignment
USAGE
  literal-dists [options] alignment.fasta[.gz] > matrix.tsv
OPTIONS
  -h    Show this help
  -v    Print version and exit
  -t    Use original snp-dists distance
  -q    Quiet mode; do not print progress information
  -a    Count all differences not just [AGTC]
  -k    Keep case, don't uppercase all letters
  -m    Output MOLTEN instead of TSV
  -c    Use comma instead of tab in output
  -b    Blank top left corner cell
  -g	Skip gap sites if gap frequency is met, gap sites [.-NX]
  -z	Gap frequency [default: 0.5]
URL
  https://github.com/kullrich/literal-dists
```

### `literal-dists -v` (version)

Prints the name and version separated by a space in standard Unix fashion.

```
literal-dists 0.0.2
```

### `literal-dists -q` (quiet mode)

Don't print informational messages, only errors.

### `literal-dists -c` (CSV instead of TSV)

```
literal-dists 0.0.1,seq1,seq2,seq3,seq4,seq5
seq1,0.000000,0.125000,0.250000,0.375000,0.500000
seq2,0.125000,0.000000,0.375000,0.500000,0.625000
seq3,0.250000,0.375000,0.000000,0.500000,0.562500
seq4,0.375000,0.500000,0.500000,0.000000,0.500000
seq5,0.500000,0.625000,0.562500,0.500000,0.000000
```

### `literal-dists -b` (omit the toolname/version)

```
	seq1	seq2	seq3	seq4	seq5
seq1	0.000000	0.125000	0.250000	0.375000	0.500000
seq2	0.125000	0.000000	0.375000	0.500000	0.625000
seq3	0.250000	0.375000	0.000000	0.500000	0.562500
seq4	0.375000	0.500000	0.500000	0.000000	0.500000
seq5	0.500000	0.625000	0.562500	0.500000	0.000000
```

### `literal-dists -t` (original snp-dists v0.7.0 distance calculation)

```
	seq1	seq2	seq3	seq4	seq5
seq1	0	1	2	3	2
seq2	1	0	3	4	3
seq3	2	3	0	4	3
seq4	3	4	4	0	2
seq5	2	3	3	2	0
```

## Advanced options (for original snp-dists implementation only)

By default, all letters are (1) uppercased and (2) ignored if not A, C, G, T.

### `literal-dists -a` (don't just count AGTC)

Normally one would not want to count ambiguous letters and gaps as a "difference"
but if you desire, you can enable this option.

```
>seq1
NGTCAGTC
>seq2
AG-CAGTC
>seq3
AGTGNGTA
```

### `literal-dists -k` (don't uppercase any letters)

You may wish to preserve case, as you may wish lower-case characters
to be masked in the comparison.
```
>seq1
AgTCAgTC
>seq2
AggCAgTC
>seq3
AgTgAgTA
```

### `literal-dists -m` ("molten" output format)

`literal-dists -m test/good.aln`

```
#seq1	seq2	dist	score	nsites
seq1	seq1	0.000000	0.000000	8.000000
seq1	seq2	0.125000	1.000000	8.000000
seq1	seq3	0.250000	2.000000	8.000000
seq1	seq4	0.375000	3.000000	8.000000
seq1	seq5	0.500000	4.000000	8.000000
seq2	seq1	0.125000	1.000000	8.000000
seq2	seq2	0.000000	0.000000	8.000000
seq2	seq3	0.375000	3.000000	8.000000
seq2	seq4	0.500000	4.000000	8.000000
seq2	seq5	0.625000	5.000000	8.000000
seq3	seq1	0.250000	2.000000	8.000000
seq3	seq2	0.375000	3.000000	8.000000
seq3	seq3	0.000000	0.000000	8.000000
seq3	seq4	0.500000	4.000000	8.000000
seq3	seq5	0.562500	4.500000	8.000000
seq4	seq1	0.375000	3.000000	8.000000
seq4	seq2	0.500000	4.000000	8.000000
seq4	seq3	0.500000	4.000000	8.000000
seq4	seq4	0.000000	0.000000	8.000000
seq4	seq5	0.500000	4.000000	8.000000
seq5	seq1	0.500000	4.000000	8.000000
seq5	seq2	0.625000	5.000000	8.000000
seq5	seq3	0.562500	4.500000	8.000000
seq5	seq4	0.500000	4.000000	8.000000
seq5	seq5	0.000000	0.000000	8.000000
```

`literal-dists -m test/ambig.aln`

```
#seq1	seq2	dist	score	nsites
seq1	seq1	0.000000	0.000000	7.000000
seq1	seq2	0.000000	0.000000	6.000000
seq1	seq3	0.333333	2.000000	6.000000
seq1	seq4	0.166667	1.000000	6.000000
seq1	seq5	0.166667	1.000000	6.000000
seq2	seq1	0.000000	0.000000	6.000000
seq2	seq2	0.000000	0.000000	7.000000
seq2	seq3	0.333333	2.000000	6.000000
seq2	seq4	0.333333	2.000000	6.000000
seq2	seq5	0.333333	2.000000	6.000000
seq3	seq1	0.333333	2.000000	6.000000
seq3	seq2	0.333333	2.000000	6.000000
seq3	seq3	0.000000	0.000000	7.000000
seq3	seq4	0.500000	3.000000	6.000000
seq3	seq5	0.500000	3.000000	6.000000
seq4	seq1	0.166667	1.000000	6.000000
seq4	seq2	0.333333	2.000000	6.000000
seq4	seq3	0.500000	3.000000	6.000000
seq4	seq4	0.000000	0.000000	7.000000
seq4	seq5	0.000000	0.000000	7.000000
seq5	seq1	0.166667	1.000000	6.000000
seq5	seq2	0.333333	2.000000	6.000000
seq5	seq3	0.500000	3.000000	6.000000
seq5	seq4	0.000000	0.000000	7.000000
seq5	seq5	0.000000	0.000000	7.000000
```

### `literal-dists -g` (skip gap sites [.-NX] if gap frequency is met)

`literal-dists -m -g test/ambig.aln`

```
#seq1	seq2	dist	score	nsites	gapsites
seq1	seq1	0.000000	0.000000	7.000000	0
seq1	seq2	0.000000	0.000000	6.000000	0
seq1	seq3	0.333333	2.000000	6.000000	0
seq1	seq4	0.166667	1.000000	6.000000	0
seq1	seq5	0.166667	1.000000	6.000000	0
seq2	seq1	0.000000	0.000000	6.000000	0
seq2	seq2	0.000000	0.000000	7.000000	0
seq2	seq3	0.333333	2.000000	6.000000	0
seq2	seq4	0.333333	2.000000	6.000000	0
seq2	seq5	0.333333	2.000000	6.000000	0
seq3	seq1	0.333333	2.000000	6.000000	0
seq3	seq2	0.333333	2.000000	6.000000	0
seq3	seq3	0.000000	0.000000	7.000000	0
seq3	seq4	0.500000	3.000000	6.000000	0
seq3	seq5	0.500000	3.000000	6.000000	0
seq4	seq1	0.166667	1.000000	6.000000	0
seq4	seq2	0.333333	2.000000	6.000000	0
seq4	seq3	0.500000	3.000000	6.000000	0
seq4	seq4	0.000000	0.000000	7.000000	0
seq4	seq5	0.000000	0.000000	7.000000	0
seq5	seq1	0.166667	1.000000	6.000000	0
seq5	seq2	0.333333	2.000000	6.000000	0
seq5	seq3	0.500000	3.000000	6.000000	0
seq5	seq4	0.000000	0.000000	7.000000	0
seq5	seq5	0.000000	0.000000	7.000000	0
```

### `literal-dists -g -z 0.2` (skip gap sites [.-NX] using gap frequency of 0.2)

`literal-dists -m -g -z 0.2 test/ambig.aln`

```
#seq1	seq2	dist	score	nsites	gapsites
seq1	seq1	0.000000	0.000000	4.000000	4
seq1	seq2	0.000000	0.000000	4.000000	4
seq1	seq3	0.500000	2.000000	4.000000	4
seq1	seq4	0.250000	1.000000	4.000000	4
seq1	seq5	0.250000	1.000000	4.000000	4
seq2	seq1	0.000000	0.000000	4.000000	4
seq2	seq2	0.000000	0.000000	4.000000	4
seq2	seq3	0.500000	2.000000	4.000000	4
seq2	seq4	0.250000	1.000000	4.000000	4
seq2	seq5	0.250000	1.000000	4.000000	4
seq3	seq1	0.500000	2.000000	4.000000	4
seq3	seq2	0.500000	2.000000	4.000000	4
seq3	seq3	0.000000	0.000000	4.000000	4
seq3	seq4	0.500000	2.000000	4.000000	4
seq3	seq5	0.500000	2.000000	4.000000	4
seq4	seq1	0.250000	1.000000	4.000000	4
seq4	seq2	0.250000	1.000000	4.000000	4
seq4	seq3	0.500000	2.000000	4.000000	4
seq4	seq4	0.000000	0.000000	4.000000	4
seq4	seq5	0.000000	0.000000	4.000000	4
seq5	seq1	0.250000	1.000000	4.000000	4
seq5	seq2	0.250000	1.000000	4.000000	4
seq5	seq3	0.500000	2.000000	4.000000	4
seq5	seq4	0.000000	0.000000	4.000000	4
seq5	seq5	0.000000	0.000000	4.000000	4
```

## Issues

Report bugs and give suggesions on the
[Issues page](https://github.com/kullrich/literal-dists/issues)

## Related software

* [snp-dists](https://github.com/tseemann/snp-dists)
* [Disty](https://github.com/c2-d2/disty)
* [Panito](https://github.com/sanger-pathogens/panito)
* [pairwise_snp_differences](https://github.com/MDU-PHL/pairwise_snp_differences/blob/master/pairwise_snp_differences.Rmd)

## Licence

[GPL Version 3](https://raw.githubusercontent.com/kullrich/literal-dists/master/LICENSE)

## Authors

* [Kristian Ullrich](https://github.com/kullrich)

## References

[https://github.com/tseemann/snp-dists](https://github.com/tseemann/snp-dists)

Chang PL., Kopania E., Keeble S., Sarver BAJ., Larson E., Orth A., Belkhir K., Boursot P., Bonhomme F., Good JM., and Dean MW. (2017). **Whole exome sequencing of wild-derived inbred strains of mice improves power to link phenotype and genotype.** *Mammalian Genome*, **28**, 416-425. [https://dx.doi.org/10.1007%2Fs00335-017-9704-9](https://dx.doi.org/10.1007%2Fs00335-017-9704-9)
