20221216-BCoAT_genome_compare
================
Jinlong Ru
12/16/22

- <a href="#visualize-using-gggenomes"
  id="toc-visualize-using-gggenomes">Visualize using
  <code>gggenomes</code></a>
  - <a href="#bcoat-in-viral-contig-baf3_v__node_17"
    id="toc-bcoat-in-viral-contig-baf3_v__node_17">BCoAT in viral contig
    BAF3_V__NODE_17</a>
  - <a href="#bcoat-in-viral-contig-bam6_v__node_3"
    id="toc-bcoat-in-viral-contig-bam6_v__node_3">BCoAT in viral contig
    BAM6_V__NODE_3</a>

**Updated: 2023-01-30 18:47:42 CET.**

# Visualize using `gggenomes`

``` r
library(gggenomes)
```

    Loading required package: dplyr


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

    Loading required package: ggplot2

    Loading required package: vctrs


    Attaching package: 'vctrs'

    The following object is masked from 'package:dplyr':

        data_frame

    Loading required package: gggenes

    Loading required package: purrr

    Loading required package: readr

    Loading required package: stringr

    Loading required package: tidyr

    Loading required package: thacklr


    Attaching package: 'thacklr'

    The following object is masked from 'package:purrr':

        %||%

    Loading required package: tibble


    Attaching package: 'tibble'

    The following object is masked from 'package:vctrs':

        data_frame

    Loading required package: jsonlite


    Attaching package: 'jsonlite'

    The following object is masked from 'package:purrr':

        flatten

    Loading required package: snakecase

    Warning: replacing previous import 'purrr::invoke' by 'rlang::invoke' when
    loading 'gggenomes'

    Warning: replacing previous import 'purrr::flatten_raw' by 'rlang::flatten_raw'
    when loading 'gggenomes'

    Warning: replacing previous import 'purrr::flatten_dbl' by 'rlang::flatten_dbl'
    when loading 'gggenomes'

    Warning: replacing previous import 'purrr::flatten_lgl' by 'rlang::flatten_lgl'
    when loading 'gggenomes'

    Warning: replacing previous import 'purrr::flatten_int' by 'rlang::flatten_int'
    when loading 'gggenomes'

    Warning: replacing previous import 'purrr::%@%' by 'rlang::%@%' when loading
    'gggenomes'

    Warning: replacing previous import 'purrr::flatten_chr' by 'rlang::flatten_chr'
    when loading 'gggenomes'

    Warning: replacing previous import 'purrr::splice' by 'rlang::splice' when
    loading 'gggenomes'

    Warning: replacing previous import 'purrr::flatten' by 'rlang::flatten' when
    loading 'gggenomes'


    Attaching package: 'gggenomes'

    The following objects are masked from 'package:thacklr':

        read_bed, read_fai, read_paf

    The following object is masked from 'package:gggenes':

        geom_gene_label

    The following object is masked from 'package:graphics':

        layout

``` r
gggenomes(emale_genes, emale_seqs, emale_tirs, emale_ava) %>%
  add_feats(ngaros=emale_ngaros, gc=emale_gc) %>%
  add_sublinks(emale_prot_ava) %>%
  flip_by_links() +
  geom_feat(position="identity", size=6) +
  geom_seq() +
  geom_link(data=links(2)) +
  geom_bin_label() +
  geom_gene(aes(fill=name)) +
  geom_gene_tag(aes(label=name), nudge_y=0.1, check_overlap = TRUE) +
  geom_feat(data=feats(ngaros), alpha=.3, size=10, position="identity") +
  geom_feat_note(aes(label="Ngaro-transposon"), feats(ngaros),
      nudge_y=.1, vjust=0) +
  geom_ribbon(aes(x=(x+xend)/2, ymax=y+.24, ymin=y+.38-(.4*score),
      group=seq_id, linetype="GC-content"), feats(gc),
      fill="lavenderblush4", position=position_nudge(y=-.1)) +
  scale_fill_brewer("Genes", palette="Dark2", na.value="cornsilk3")
```

    Transforming sublinks with "aa2nuc". Disable with `.transform = "none"`

    Flipping: 2,3,4

    Warning: Removed 95 rows containing missing values (`geom_feat_text()`).

![](20221216-BCoAT_genome_compare_files/figure-commonmark/unnamed-chunk-1-1.png)

``` r
ggsave(path_target("emales.png"), width=8, height=4)
```

    Warning: Removed 95 rows containing missing values (`geom_feat_text()`).

------------------------------------------------------------------------

**Updated: 2022-12-16**

Two Viral BCoAT protein sequences were searched against NCBI NR
database. The top protein hit in bacterial reference genome was selected
as the best target, and the bacterial genome was compared with viral
genome as show below.

### BCoAT in viral contig BAF3_V\_\_NODE_17

Bacterial `Pseudoflavonifractor sp. AF19-9AC` contains BCoAT that most
similar to viral contig `BAF3_V__NODE_17` encoded BCoAT gene, so we
compared BCoAT region of their genomes, as show below. BCoAT gene was
colored in red, and was connected with lines according to identity by
BLASTn alignment. The high identity in BCoAT gene and no identity in
surrounding genes indicate the viral-encoded BCoAT gene was from
bacterial genome. Upper part is bacterial genome, and lower part is
viral contig.

![Figure 1. BCoAT in Pseudoflavonifractor sp. AF19-9AC and
BAF3_V\_\_NODE_17](data/00-rawdata/BAF3.svg)

### BCoAT in viral contig BAM6_V\_\_NODE_3

Bacterial `PLawsonibacter sp. OA9` contains BCoAT that most similar to
viral contig `BAM6_V__NODE_3` encoded BCoAT gene, so we compared BCoAT
region of their genomes, as show below. BCoAT gene was colored in red,
and was connected with lines according to identity by BLASTn alignment.
The high identity in BCoAT gene and no identity in surrounding genes
indicate the viral-encoded BCoAT gene was from bacterial genome. Upper
part is bacterial genome, and lower part is viral contig.

![Figure 2. BCoAT in PLawsonibacter sp. OA9 and
BAM6_V\_\_NODE_3](data/00-rawdata/BAM6.svg)
