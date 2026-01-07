# Translate nucleic acid sequences into codons

The function transforms `dna2codon` DNA sequences to codon sequences,
`codon2dna` transform the other way. `dna2codon` translates nucleotide
to amino acids using the function
[`trans`](https://rdrr.io/pkg/ape/man/trans.html).

## Usage

``` r
dna2codon(x, codonstart = 1, code = 1, ambiguity = "---", ...)

codon2dna(x)

dna2aa(x, codonstart = 1, code = 1)
```

## Arguments

- x:

  An object containing sequences.

- codonstart:

  an integer giving where to start the translation. This should be 1, 2,
  or 3, but larger values are accepted and have for effect to start the
  translation further within the sequence.

- code:

  The ncbi genetic code number for translation (see details). By default
  the standard genetic code is used.

- ambiguity:

  character for ambiguous character and no contrast is provided.

- ...:

  further arguments passed to or from other methods.

## Value

The functions return an object of class `phyDat`.

## Details

The following genetic codes are described here. The number preceding
each corresponds to the code argument.

|     |                                       |
|-----|---------------------------------------|
| 1   | standard                              |
| 2   | vertebrate.mitochondrial              |
| 3   | yeast.mitochondrial                   |
| 4   | protozoan.mitochondrial+mycoplasma    |
| 5   | invertebrate.mitochondrial            |
| 6   | ciliate+dasycladaceal                 |
| 9   | echinoderm+flatworm.mitochondrial     |
| 10  | euplotid                              |
| 11  | bacterial+plantplastid                |
| 12  | alternativeyeast                      |
| 13  | ascidian.mitochondrial                |
| 14  | alternativeflatworm.mitochondrial     |
| 15  | blepharism                            |
| 16  | chlorophycean.mitochondrial           |
| 21  | trematode.mitochondrial               |
| 22  | scenedesmus.mitochondrial             |
| 23  | thraustochytrium.mitochondria         |
| 24  | Pterobranchia.mitochondrial           |
| 25  | CandidateDivision.SR1+Gracilibacteria |
| 26  | Pachysolen.tannophilus                |

Alignment gaps and ambiguities are currently ignored and sites
containing these are deleted.

## References

<https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes>

## See also

[`trans`](https://rdrr.io/pkg/ape/man/trans.html),
[`phyDat`](https://klausvigo.github.io/phangorn/reference/as.phyDat.md)
and the chapter 4 in the
`vignette("phangorn-specials", package="phangorn")`

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
data(Laurasiatherian)
class(Laurasiatherian)
#> [1] "phyDat"
Laurasiatherian
#> 47 sequences with 3179 character and 1605 different site patterns.
#> The states are a c g t 
dna2codon(Laurasiatherian)
#> 47 sequences with 1059 character and 914 different site patterns.
#> The states are aaa aac aag aat aca acc acg act aga agc agg agt ata atc atg att caa cac cag cat cca ccc ccg cct cga cgc cgg cgt cta ctc ctg ctt gaa gac gag gat gca gcc gcg gct gga ggc ggg ggt gta gtc gtg gtt tac tat tca tcc tcg tct tgc tgg tgt tta ttc ttg ttt 
```
