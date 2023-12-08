


# Batch Effect Correction {#BatchEffectCorrection}

批次效应是除处理其他因素带来的影响实验结果的效应，比如在研究对照组和实验组时候，提取样本DNA时间不同、送测仪器不同等均可能引入批次效应。批次效应不能够消除，只可以降低，使得它对所研究的生物学问题有较小的影响。
评估批次效应可以通过PCA或PCoA等降纬可视化技术观察样本的聚集情况是否与非研究问题存在视觉上的相关性，也可以通过一些统计方法如PERMANOVA等计算表达谱整体水平和因素之间的关系。


## MMUPHin

通常在做多个研究或多个平台的数据整合需要考虑到消除不同研究或平台数据差异，在微生物领域又因为微生物相对丰度数据是稀疏的，所以常在转录组领域使用的校正方法如`sva::ComBat`和`limma::removeBatchEffect`等均不适用。

"MMUPHin"（Meta-analysis via Mixed Models Utilizing Public Health Information）是一个哈佛大学Huttenhover实验室开发的用于微生物组数据的统计分析包，特别是在研究与公共健康相关的多个研究的数据时使用。它在处理批次效应（batch effects）时的原理是基于混合模型（mixed models）。

在微生物组学研究中，批次效应是一个常见的问题。它指的是由于样本处理和测序过程中的技术变异而导致的非生物学差异，这些差异可能会干扰真实的生物学信号。例如，不同的实验室使用不同的样本处理方法或测序平台，可能导致数据之间的系统性差异。


MMUPHin处理批次效应的原理：

+ 混合模型：MMUPHin使用混合模型来纳入批次效应。在这种模型中，批次效应被视为随机效应，它们与研究中的固定效应（例如治疗组与对照组）分开处理。

+ 元分析方法：MMUPHin利用元分析的技术，允许来自不同研究的数据共同分析。元分析是一种统计方法，它综合并分析多个研究的结果，以获得更广泛、更全面的结论。

+ 数据整合：通过混合模型和元分析方法，MMUPHin能够在考虑批次效应的同时，整合多个研究的数据，提高分析的统计能力和结论的泛化性。

+ 校正批次效应：通过在模型中包括批次效应作为一个变量，MMUPHin可以校正这些非生物学差异，从而使研究结果更加可靠和准确。


### 加载依赖包和数据


+ 数据来源是 `curatedMetagenomicData` R包


```r
knitr::opts_chunk$set(message = FALSE, warning = FALSE)

library(tidyverse)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("MMUPHin")

library(MMUPHin)
library(ggplot2)
library(phyloseq)

data("CRC_abd", "CRC_meta")
CRC_abd[1:5, 1, drop = FALSE]
#>                                 FengQ_2015.metaphlan_bugs_list.stool:SID31004
#> s__Faecalibacterium_prausnitzii                                    0.11110668
#> s__Streptococcus_salivarius                                        0.09660736
#> s__Ruminococcus_sp_5_1_39BFAA                                      0.09115385
#> s__Subdoligranulum_unclassified                                    0.05806767
#> s__Bacteroides_stercoris                                           0.05685503
# CRC_meta[1, 1:5]
```


### 数据探索

使用`MicrobiomeAnalysis::run_ord`和`MicrobiomeAnalysis::plot_ord`可视化数据，选择PCoA的方法。

+ 安装MicrobiomeAnalysis包
```R
if (!requireNamespace(c("remotes", "devtools"), quietly=TRUE)) {
  install.packages(c("devtools", "remotes"))
}
remotes::install_github("HuaZou/MicrobiomeAnalysis")
```


```r
# run_ord需要phyloseq object
rownames(CRC_meta) <- make.names(rownames(CRC_meta))
colnames(CRC_abd) <- make.names(colnames(CRC_abd))
rownames(CRC_abd) <- make.names(rownames(CRC_abd))
tax_tab <- data.frame(Species = rownames(CRC_abd))
rownames(tax_tab) <- tax_tab$Species

CRC_phy <- phyloseq::phyloseq(
  sample_data(CRC_meta),
  otu_table(CRC_abd, taxa_are_rows = TRUE),
  tax_table(as.matrix(tax_tab))
)

# 选择两个研究
CRC_phy_new <- CRC_phy
phyloseq::sample_data(CRC_phy_new) <- phyloseq::sample_data(
  CRC_phy@sam_data %>%
    data.frame() %>%
    dplyr::filter(studyID %in% c("FengQ_2015.metaphlan_bugs_list.stool",
                                 "VogtmannE_2016.metaphlan_bugs_list.stool"))) 
# 运行
ord_result <- MicrobiomeAnalysis::run_ord(
  object = CRC_phy_new,
  variable = "study_condition",
  method = "PCoA")
MicrobiomeAnalysis::plot_ord(
  reslist = ord_result,
  variable = "study_condition",
  variable_color = c("red", "blue"),
  var_shape = "studyID",
  ellipse_type = "none")
```

<img src="12-Statistics_BatchEffect_files/figure-html/unnamed-chunk-3-1.png" width="100%" />

结果：能明显看到样本以`studyID`分散开来，研究之间的批次效应要远远大于control和CRC的生物学效应，接下来我们通过PERMANOVA做统计分析进一步评估显著性。


```r
MicrobiomeAnalysis::run_PERMANOVA(
  CRC_phy_new,
  variables = c("study_condition", "studyID"),
  mode = "one",
  method = "bray")
#>                 SumsOfSample Df SumsOfSqs   MeanSqs
#> study_condition          211  1 0.6756486 0.6756486
#> studyID                  211  1 5.6815343 5.6815343
#>                   F.Model         R2 Pr(>F) AdjustedPvalue
#> study_condition  2.295539 0.01086411  0.004          0.004
#> studyID         21.013183 0.09135643  0.001          0.002
```
结果：PERMANOVA的`Pr(>F)` < 0.05表明study_condition和studyID均与整体肠道微生物结构有显著差异，也说明studyID的批次效应(_对肠道结构解释9.1%的变异_)对后续study_condition的差异研究等具有非常大的影响，因此需要做批次校正。


### 批次效应校正

需要明确的一点是批次效应只能降低，不能完全消除，并且在做批次效应过程中可能会降低或提高所研究的生物学意义，这是因为使用线性模型校正所带来的结果。MMUPHin处理批次效应的原理：

+ 混合模型：MMUPHin使用混合模型来纳入批次效应。在这种模型中，批次效应被视为随机效应，它们与研究中的固定效应（例如治疗组与对照组）分开处理。

混合模型可以参考[GEE and MLM](https://zouhua.top/DraftNotes/GEEandMLM.html)，它提供了常用的两种混合模型GEE和MLM。MMUPHin采用的是**Zero-inflated empirical Bayes adjustment of batch effect in compositional feature abundance data**。


```r
fit_adjust_batch <- adjust_batch(
  feature_abd = CRC_abd,
  batch = "studyID",
  covariates = "study_condition",
  data = CRC_meta,
  control = list(verbose = FALSE))

CRC_abd_adj <- fit_adjust_batch$feature_abd_adj

CRC_abd_adj[1:5, 1, drop = FALSE]
#>                                 FengQ_2015.metaphlan_bugs_list.stool.SID31004
#> s__Faecalibacterium_prausnitzii                                    0.10120482
#> s__Streptococcus_salivarius                                        0.06044265
#> s__Ruminococcus_sp_5_1_39BFAA                                      0.02374596
#> s__Subdoligranulum_unclassified                                    0.03265566
#> s__Bacteroides_stercoris                                           0.31510103
```

查看校正后的结果


```r
# run_ord需要phyloseq object
rownames(CRC_meta) <- make.names(rownames(CRC_meta))
colnames(CRC_abd_adj) <- make.names(colnames(CRC_abd_adj))
rownames(CRC_abd_adj) <- make.names(rownames(CRC_abd_adj))
tax_tab_adj <- data.frame(Species = rownames(CRC_abd_adj))
rownames(tax_tab_adj) <- tax_tab_adj$Species

CRC_phy_adj <- phyloseq::phyloseq(
  sample_data(CRC_meta),
  otu_table(CRC_abd_adj, taxa_are_rows = TRUE),
  tax_table(as.matrix(tax_tab_adj))
)

# 选择两个研究
CRC_phy_adj_new <- CRC_phy_adj
phyloseq::sample_data(CRC_phy_adj_new) <- phyloseq::sample_data(
  CRC_phy_adj@sam_data %>%
    data.frame() %>%
    dplyr::filter(studyID %in% c("FengQ_2015.metaphlan_bugs_list.stool",
                                 "VogtmannE_2016.metaphlan_bugs_list.stool"))) 
# 运行
ord_result <- MicrobiomeAnalysis::run_ord(
  object = CRC_phy_adj_new,
  variable = "study_condition",
  method = "PCoA")
MicrobiomeAnalysis::plot_ord(
  reslist = ord_result,
  variable = "study_condition",
  variable_color = c("red", "blue"),
  var_shape = "studyID",
  ellipse_type = "none")
```

<img src="12-Statistics_BatchEffect_files/figure-html/unnamed-chunk-6-1.png" width="100%" />

结果：相比校正前，studyID带来的效应明显降低。进一步通过PERMANOVA结果分析。


```r
MicrobiomeAnalysis::run_PERMANOVA(
  CRC_phy_adj_new,
  variables = c("study_condition", "studyID"),
  mode = "one",
  method = "bray")
#>                 SumsOfSample Df SumsOfSqs   MeanSqs
#> study_condition          211  1 0.7256598 0.7256598
#> studyID                  211  1 1.0133205 1.0133205
#>                  F.Model         R2 Pr(>F) AdjustedPvalue
#> study_condition 2.511293 0.01187309  0.002          0.002
#> studyID         3.523584 0.01657973  0.001          0.002
```

结果：相比校正前，studyID的解释肠道结构总体变异度从9.1%降低到了1.7%。与此同时，study_condition的总体变异则几乎没有太大变化。


### 荟萃分析

荟萃分析的目的是汇总各个研究的共同结果进而获得一个共有的效果，先前也有很多工具提供类似的研究，例如meta包。MMUPHin也提供了`lm_meta`函数用于分析。

+ 先用Maaslin2计算不同研究在control和CRC之间的差异物种；

+ 再使用混合模型汇总所有的结果。

因为采用了线性回归方式，所以可以加入一些协变量如年龄、性别和BMI等人口统计变量作为校正因素。coef表示EffectSize也即是在Control和CRC组间的区别。


```r
if(!dir.exists("./InputData/MMUPHin_lm_meta")) {
  dir.create("./InputData/MMUPHin_lm_meta", recursive = TRUE)
}

fit_lm_meta <- lm_meta(
  feature_abd = CRC_abd_adj,
  batch = "studyID",
  exposure = "study_condition",
  covariates = c("gender", "age", "BMI"),
  data = CRC_meta,
  control = list(verbose = FALSE, 
                 output = "./InputData/MMUPHin_lm_meta"))

fit_lm_meta$meta_fits %>% 
  filter(qval.fdr < 0.05) %>% 
  arrange(coef) %>% 
  mutate(feature = factor(feature, levels = feature)) %>% 
  ggplot(aes(y = coef, x = feature)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw()
```

<img src="12-Statistics_BatchEffect_files/figure-html/unnamed-chunk-8-1.png" width="100%" />


### 总结

+ MMUPHin提供了可以校正多个数据来源的批次效应函数

+ MMUPHin在做荟萃分析的时提供了工具


## Systemic information

```r
devtools::session_info()
#> ─ Session info ───────────────────────────────────────────
#>  setting  value
#>  version  R version 4.1.3 (2022-03-10)
#>  os       macOS Big Sur/Monterey 10.16
#>  system   x86_64, darwin17.0
#>  ui       X11
#>  language (EN)
#>  collate  en_US.UTF-8
#>  ctype    en_US.UTF-8
#>  tz       Asia/Shanghai
#>  date     2023-12-08
#>  pandoc   3.1.3 @ /Users/zouhua/opt/anaconda3/bin/ (via rmarkdown)
#> 
#> ─ Packages ───────────────────────────────────────────────
#>  package              * version    date (UTC) lib source
#>  ade4                   1.7-22     2023-02-06 [2] CRAN (R 4.1.2)
#>  ANCOMBC                1.4.0      2021-10-26 [2] Bioconductor
#>  annotate               1.72.0     2021-10-26 [2] Bioconductor
#>  AnnotationDbi          1.60.2     2023-03-10 [2] Bioconductor
#>  ape                    5.7-1      2023-03-13 [2] CRAN (R 4.1.2)
#>  biglm                  0.9-2.1    2020-11-27 [2] CRAN (R 4.1.0)
#>  Biobase                2.54.0     2021-10-26 [2] Bioconductor
#>  BiocGenerics           0.40.0     2021-10-26 [2] Bioconductor
#>  BiocParallel           1.28.3     2021-12-09 [2] Bioconductor
#>  biomformat             1.22.0     2021-10-26 [2] Bioconductor
#>  Biostrings             2.62.0     2021-10-26 [2] Bioconductor
#>  bit                    4.0.5      2022-11-15 [2] CRAN (R 4.1.2)
#>  bit64                  4.0.5      2020-08-30 [2] CRAN (R 4.1.0)
#>  bitops                 1.0-7      2021-04-24 [2] CRAN (R 4.1.0)
#>  blob                   1.2.4      2023-03-17 [2] CRAN (R 4.1.2)
#>  bookdown               0.34       2023-05-09 [2] CRAN (R 4.1.2)
#>  bslib                  0.6.0      2023-11-21 [1] CRAN (R 4.1.3)
#>  cachem                 1.0.8      2023-05-01 [2] CRAN (R 4.1.2)
#>  callr                  3.7.3      2022-11-02 [2] CRAN (R 4.1.2)
#>  caTools                1.18.2     2021-03-28 [2] CRAN (R 4.1.0)
#>  cli                    3.6.1      2023-03-23 [2] CRAN (R 4.1.2)
#>  cluster                2.1.4      2022-08-22 [2] CRAN (R 4.1.2)
#>  codetools              0.2-19     2023-02-01 [2] CRAN (R 4.1.2)
#>  colorspace             2.1-0      2023-01-23 [2] CRAN (R 4.1.2)
#>  cowplot                1.1.1      2020-12-30 [2] CRAN (R 4.1.0)
#>  crayon                 1.5.2      2022-09-29 [2] CRAN (R 4.1.2)
#>  data.table             1.14.8     2023-02-17 [2] CRAN (R 4.1.2)
#>  DBI                    1.1.3      2022-06-18 [2] CRAN (R 4.1.2)
#>  DelayedArray           0.20.0     2021-10-26 [2] Bioconductor
#>  DEoptimR               1.0-14     2023-06-09 [2] CRAN (R 4.1.3)
#>  DESeq2                 1.34.0     2021-10-26 [2] Bioconductor
#>  devtools               2.4.5      2022-10-11 [2] CRAN (R 4.1.2)
#>  digest                 0.6.33     2023-07-07 [1] CRAN (R 4.1.3)
#>  downlit                0.4.3      2023-06-29 [2] CRAN (R 4.1.3)
#>  dplyr                * 1.1.2      2023-04-20 [2] CRAN (R 4.1.2)
#>  ellipsis               0.3.2      2021-04-29 [2] CRAN (R 4.1.0)
#>  evaluate               0.21       2023-05-05 [2] CRAN (R 4.1.2)
#>  fansi                  1.0.4      2023-01-22 [2] CRAN (R 4.1.2)
#>  farver                 2.1.1      2022-07-06 [2] CRAN (R 4.1.2)
#>  fastmap                1.1.1      2023-02-24 [2] CRAN (R 4.1.2)
#>  forcats              * 1.0.0      2023-01-29 [2] CRAN (R 4.1.2)
#>  foreach                1.5.2      2022-02-02 [2] CRAN (R 4.1.2)
#>  fs                     1.6.2      2023-04-25 [2] CRAN (R 4.1.2)
#>  genefilter             1.76.0     2021-10-26 [2] Bioconductor
#>  geneplotter            1.72.0     2021-10-26 [2] Bioconductor
#>  generics               0.1.3      2022-07-05 [2] CRAN (R 4.1.2)
#>  GenomeInfoDb           1.30.1     2022-01-30 [2] Bioconductor
#>  GenomeInfoDbData       1.2.7      2022-03-09 [2] Bioconductor
#>  GenomicRanges          1.46.1     2021-11-18 [2] Bioconductor
#>  getopt                 1.20.3     2019-03-22 [2] CRAN (R 4.1.0)
#>  ggplot2              * 3.4.2      2023-04-03 [2] CRAN (R 4.1.2)
#>  glmnet                 4.1-7      2023-03-23 [2] CRAN (R 4.1.2)
#>  glue                   1.6.2      2022-02-24 [2] CRAN (R 4.1.2)
#>  gplots                 3.1.3      2022-04-25 [2] CRAN (R 4.1.2)
#>  gtable                 0.3.3      2023-03-21 [2] CRAN (R 4.1.2)
#>  gtools                 3.9.4      2022-11-27 [2] CRAN (R 4.1.2)
#>  hash                   2.2.6.2    2022-03-22 [2] CRAN (R 4.1.2)
#>  highr                  0.10       2022-12-22 [2] CRAN (R 4.1.2)
#>  hms                    1.1.3      2023-03-21 [2] CRAN (R 4.1.2)
#>  htmltools              0.5.7      2023-11-03 [1] CRAN (R 4.1.3)
#>  htmlwidgets            1.6.2      2023-03-17 [2] CRAN (R 4.1.2)
#>  httpuv                 1.6.11     2023-05-11 [2] CRAN (R 4.1.3)
#>  httr                   1.4.6      2023-05-08 [2] CRAN (R 4.1.2)
#>  igraph                 1.5.0      2023-06-16 [1] CRAN (R 4.1.3)
#>  IRanges                2.28.0     2021-10-26 [2] Bioconductor
#>  iterators              1.0.14     2022-02-05 [2] CRAN (R 4.1.2)
#>  jquerylib              0.1.4      2021-04-26 [2] CRAN (R 4.1.0)
#>  jsonlite               1.8.7      2023-06-29 [2] CRAN (R 4.1.3)
#>  KEGGREST               1.34.0     2021-10-26 [2] Bioconductor
#>  KernSmooth             2.23-22    2023-07-10 [2] CRAN (R 4.1.3)
#>  knitr                  1.43       2023-05-25 [2] CRAN (R 4.1.3)
#>  labeling               0.4.2      2020-10-20 [2] CRAN (R 4.1.0)
#>  later                  1.3.1      2023-05-02 [2] CRAN (R 4.1.2)
#>  lattice                0.21-8     2023-04-05 [2] CRAN (R 4.1.2)
#>  lifecycle              1.0.3      2022-10-07 [2] CRAN (R 4.1.2)
#>  limma                  3.50.3     2022-04-07 [2] Bioconductor
#>  locfit                 1.5-9.8    2023-06-11 [2] CRAN (R 4.1.3)
#>  logging                0.10-108   2019-07-14 [2] CRAN (R 4.1.0)
#>  lpsymphony             1.22.0     2021-10-26 [2] Bioconductor (R 4.1.1)
#>  lubridate            * 1.9.2      2023-02-10 [2] CRAN (R 4.1.2)
#>  Maaslin2               1.8.0      2021-10-26 [2] Bioconductor
#>  magrittr               2.0.3      2022-03-30 [2] CRAN (R 4.1.2)
#>  MASS                   7.3-60     2023-05-04 [2] CRAN (R 4.1.2)
#>  mathjaxr               1.6-0      2022-02-28 [2] CRAN (R 4.1.2)
#>  Matrix                 1.6-0      2023-07-08 [2] CRAN (R 4.1.3)
#>  MatrixGenerics         1.6.0      2021-10-26 [2] Bioconductor
#>  matrixStats            1.0.0      2023-06-02 [2] CRAN (R 4.1.3)
#>  memoise                2.0.1      2021-11-26 [2] CRAN (R 4.1.0)
#>  metadat                1.2-0      2022-04-06 [2] CRAN (R 4.1.2)
#>  metafor                4.2-0      2023-05-08 [2] CRAN (R 4.1.2)
#>  metagenomeSeq          1.36.0     2021-10-26 [2] Bioconductor
#>  mgcv                   1.8-42     2023-03-02 [2] CRAN (R 4.1.2)
#>  microbiome             1.16.0     2021-10-26 [2] Bioconductor
#>  MicrobiomeAnalysis     1.0.3      2023-12-02 [1] Bioconductor
#>  mime                   0.12       2021-09-28 [2] CRAN (R 4.1.0)
#>  miniUI                 0.1.1.1    2018-05-18 [2] CRAN (R 4.1.0)
#>  MMUPHin              * 1.8.2      2022-04-03 [1] Bioconductor
#>  multtest               2.50.0     2021-10-26 [2] Bioconductor
#>  munsell                0.5.0      2018-06-12 [2] CRAN (R 4.1.0)
#>  mvtnorm                1.2-2      2023-06-08 [2] CRAN (R 4.1.3)
#>  nlme                   3.1-162    2023-01-31 [2] CRAN (R 4.1.2)
#>  nloptr                 2.0.3      2022-05-26 [2] CRAN (R 4.1.2)
#>  numDeriv               2016.8-1.1 2019-06-06 [2] CRAN (R 4.1.0)
#>  optparse               1.7.3      2022-07-20 [2] CRAN (R 4.1.2)
#>  pbapply                1.7-2      2023-06-27 [2] CRAN (R 4.1.3)
#>  pcaPP                  2.0-3      2022-10-24 [2] CRAN (R 4.1.2)
#>  permute                0.9-7      2022-01-27 [2] CRAN (R 4.1.2)
#>  phyloseq             * 1.38.0     2021-10-26 [2] Bioconductor
#>  pillar                 1.9.0      2023-03-22 [2] CRAN (R 4.1.2)
#>  pkgbuild               1.4.2      2023-06-26 [2] CRAN (R 4.1.3)
#>  pkgconfig              2.0.3      2019-09-22 [2] CRAN (R 4.1.0)
#>  pkgload                1.3.2.1    2023-07-08 [2] CRAN (R 4.1.3)
#>  plyr                   1.8.8      2022-11-11 [2] CRAN (R 4.1.2)
#>  png                    0.1-8      2022-11-29 [2] CRAN (R 4.1.2)
#>  prettyunits            1.1.1      2020-01-24 [2] CRAN (R 4.1.0)
#>  processx               3.8.2      2023-06-30 [2] CRAN (R 4.1.3)
#>  profvis                0.3.8      2023-05-02 [2] CRAN (R 4.1.2)
#>  promises               1.2.0.1    2021-02-11 [2] CRAN (R 4.1.0)
#>  ps                     1.7.5      2023-04-18 [2] CRAN (R 4.1.2)
#>  purrr                * 1.0.1      2023-01-10 [2] CRAN (R 4.1.2)
#>  R6                     2.5.1      2021-08-19 [2] CRAN (R 4.1.0)
#>  ragg                   1.2.5      2023-01-12 [2] CRAN (R 4.1.2)
#>  rbibutils              2.2.13     2023-01-13 [2] CRAN (R 4.1.2)
#>  RColorBrewer           1.1-3      2022-04-03 [2] CRAN (R 4.1.2)
#>  Rcpp                   1.0.11     2023-07-06 [1] CRAN (R 4.1.3)
#>  RCurl                  1.98-1.12  2023-03-27 [2] CRAN (R 4.1.2)
#>  Rdpack                 2.4        2022-07-20 [2] CRAN (R 4.1.2)
#>  readr                * 2.1.4      2023-02-10 [2] CRAN (R 4.1.2)
#>  remotes                2.4.2      2021-11-30 [2] CRAN (R 4.1.0)
#>  reshape2               1.4.4      2020-04-09 [2] CRAN (R 4.1.0)
#>  rhdf5                  2.38.1     2022-03-10 [2] Bioconductor
#>  rhdf5filters           1.6.0      2021-10-26 [2] Bioconductor
#>  Rhdf5lib               1.16.0     2021-10-26 [2] Bioconductor
#>  rlang                  1.1.1      2023-04-28 [1] CRAN (R 4.1.2)
#>  rmarkdown              2.23       2023-07-01 [2] CRAN (R 4.1.3)
#>  robustbase             0.99-0     2023-06-16 [2] CRAN (R 4.1.3)
#>  RSQLite                2.3.1      2023-04-03 [2] CRAN (R 4.1.2)
#>  rstudioapi             0.15.0     2023-07-07 [2] CRAN (R 4.1.3)
#>  Rtsne                  0.16       2022-04-17 [2] CRAN (R 4.1.2)
#>  S4Vectors              0.32.4     2022-03-29 [2] Bioconductor
#>  sass                   0.4.6      2023-05-03 [2] CRAN (R 4.1.2)
#>  scales                 1.2.1      2022-08-20 [2] CRAN (R 4.1.2)
#>  sessioninfo            1.2.2      2021-12-06 [2] CRAN (R 4.1.0)
#>  shape                  1.4.6      2021-05-19 [2] CRAN (R 4.1.0)
#>  shiny                  1.7.4.1    2023-07-06 [2] CRAN (R 4.1.3)
#>  stringi                1.7.12     2023-01-11 [2] CRAN (R 4.1.2)
#>  stringr              * 1.5.0      2022-12-02 [2] CRAN (R 4.1.2)
#>  SummarizedExperiment   1.24.0     2021-10-26 [2] Bioconductor
#>  survival               3.5-5      2023-03-12 [2] CRAN (R 4.1.2)
#>  systemfonts            1.0.4      2022-02-11 [2] CRAN (R 4.1.2)
#>  textshaping            0.3.6      2021-10-13 [2] CRAN (R 4.1.0)
#>  tibble               * 3.2.1      2023-03-20 [2] CRAN (R 4.1.2)
#>  tidyr                * 1.3.0      2023-01-24 [2] CRAN (R 4.1.2)
#>  tidyselect             1.2.0      2022-10-10 [2] CRAN (R 4.1.2)
#>  tidyverse            * 2.0.0      2023-02-22 [1] CRAN (R 4.1.2)
#>  timechange             0.2.0      2023-01-11 [2] CRAN (R 4.1.2)
#>  tzdb                   0.4.0      2023-05-12 [2] CRAN (R 4.1.3)
#>  urlchecker             1.0.1      2021-11-30 [2] CRAN (R 4.1.0)
#>  usethis                2.2.2      2023-07-06 [2] CRAN (R 4.1.3)
#>  utf8                   1.2.3      2023-01-31 [2] CRAN (R 4.1.2)
#>  vctrs                  0.6.3      2023-06-14 [1] CRAN (R 4.1.3)
#>  vegan                  2.6-4      2022-10-11 [2] CRAN (R 4.1.2)
#>  withr                  2.5.0      2022-03-03 [2] CRAN (R 4.1.2)
#>  Wrench                 1.12.0     2021-10-26 [2] Bioconductor
#>  xfun                   0.40       2023-08-09 [1] CRAN (R 4.1.3)
#>  XML                    3.99-0.14  2023-03-19 [2] CRAN (R 4.1.2)
#>  xml2                   1.3.5      2023-07-06 [2] CRAN (R 4.1.3)
#>  xtable                 1.8-4      2019-04-21 [2] CRAN (R 4.1.0)
#>  XVector                0.34.0     2021-10-26 [2] Bioconductor
#>  yaml                   2.3.7      2023-01-23 [2] CRAN (R 4.1.2)
#>  zlibbioc               1.40.0     2021-10-26 [2] Bioconductor
#> 
#>  [1] /Users/zouhua/Library/R/x86_64/4.1/library
#>  [2] /Library/Frameworks/R.framework/Versions/4.1/Resources/library
#> 
#> ──────────────────────────────────────────────────────────
```


## Reference

+ [MMUPHin](https://huttenhower.sph.harvard.edu/mmuphin)

+ [Doing Meta-Analysis with R: A Hands-On Guide](https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/)
