


# Dimension Reduction {#DimensionReduction}


降维分析是为了查看样本的差异或挖掘潜在标志物。方法有以下三种：

+ 主成分分析（PCA）

主成分分析PCA（Principal Component Analysis）是一种常用的数据分析方法。PCA通过线性变换将原始数据变换为一组各维度线性无关的表示，可用于提取数据的主要特征分量，常用于高维数据的降维。

+ 偏最小二乘回归分析法（PLS-DA）

偏最小二乘判别分析（PLS-DA）是一种有监督模式识别的多元统计分析方法，将多维数据在压缩前先按需要寻找的差异因素分组。

+ 正交偏最小二乘判别分析（OPLS-DA）

正交偏最小二乘判别分析（OPLS-DA）结合了正交信号矫正（OSC）和PLS-DA方法，能够将X矩阵信息分解成与Y相关和不相关的两类信息，通过去除不相关的差异来筛选差异变量。


## 加载R包

```r
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
library(tidyverse)
library(ropls)
library(SummarizedExperiment)
library(MicrobiomeAnalysis)

# rm(list = ls())
options(stringsAsFactors = F)
options(future.globals.maxSize = 1000 * 1024^2)

grp_names <- c("None", "Mild", "Moderate", "Severe")
grp_colors <- c("#7DD06F", "#844081", "#688EC1", "#C17E73")
```


## 导入数据

对数据[OmicsDataSet-Zeybel et al. - 2022.xlsx](https://github.com/HuaZou/DraftNotes/blob/main/InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx)处理后生成的，可参考数据预处理章节。

> ```R
> saveRDS(se_scale, "./InputData/result/QC/se_scale.RDS", compress = TRUE)
> ```


```r
data_meta <- readRDS("./InputData/result/QC/se_scale.RDS")
```


## 函数

```r
DR_fun <- function(
    x,
    group,
    group_names,
    group_colors,
    DRtype = c("PCA", "PLS", "OPLS"),
    occ_cutoff = 0.5) {
  
  # x = data_meta
  # group = "LiverFatClass"
  # group_names = grp_names
  # group_colors = grp_colors
  # DRtype = "PLS" # PCA
  # occ_cutoff = 0.5
  
  # dataseat
  metadata <- SummarizedExperiment::colData(x) %>%
      as.data.frame()
  profile <- SummarizedExperiment::assay(x) %>%
      as.data.frame()

  colnames(metadata)[which(colnames(metadata) == group)] <- "CompVar"
  phenotype <- metadata %>%
    dplyr::filter(CompVar %in% group_names) %>%
    dplyr::mutate(CompVar = as.character(CompVar)) %>%
    dplyr::mutate(CompVar = factor(CompVar, levels = group_names))
  
  sid <- intersect(rownames(phenotype), colnames(profile))
  phen <- phenotype[pmatch(sid, rownames(phenotype)), , ]
  prof <- profile %>%
    dplyr::select(all_of(sid))
  
  if (!all(colnames(prof) == rownames(phen))) {
    stop("Wrong Order")
  }
  
  trim_FeatureOrSample <- function(x, nRow, threshold) {
  
    df_occ <- apply(x, nRow, function(x) {
      length(x[c(which(!is.na(x) & x!=0))]) / length(x)
    }) %>%
      data.frame() %>% stats::setNames("Occ") %>%
      tibble::rownames_to_column("type")
    if(nRow == 1){
      rownames(df_occ) <- rownames(x)
    }else{
      rownames(df_occ) <- colnames(x)
    }
    df_KEEP <- apply(df_occ > threshold, 1, all) %>%
      data.frame() %>% stats::setNames("Status") %>%
      dplyr::filter(Status)
  
    res <- x %>%
      tibble::rownames_to_column("featureid") %>%
      dplyr::filter(featureid %in% rownames(df_KEEP)) %>%
      tibble::column_to_rownames("featureid")
    
    return(res)
  }

  prof_cln <- trim_FeatureOrSample(prof, 1, occ_cutoff)
  dataMatrix <- prof_cln %>% t() # row->sampleID; col->features
  sampleMetadata <- phen %>% # row->sampleID; col->features
    dplyr::mutate(CompVar = factor(CompVar, levels = group_names)) %>%
    dplyr::mutate(Color = factor(CompVar, 
                                  levels = group_names,
                                  labels = group_colors),
                  Color = as.character(Color))    
  
  if (DRtype == "PCA") {
    fit <- opls(dataMatrix)
    plot(fit,
         typeVc = "x-score",
         parAsColFcVn = sampleMetadata$CompVar,
         )
  } else if (DRtype == "PLS") {
    fit <- opls(dataMatrix, sampleMetadata$CompVar, predI = 2)
    plot(fit,
         typeVc = "x-score",
         parAsColFcVn = sampleMetadata$CompVar,
         )
  } else if (DRtype == "OPLS") {
    # only for binary classification
    fit <- opls(dataMatrix, sampleMetadata$CompVar, predI = 1, orthoI = NA)
    plot(fit,
         typeVc = "x-score",
         parAsColFcVn = sampleMetadata$CompVar,
         )
  }

  return(fit)
}
```


## 主成分分析（PCA）

```r
PCA_res <- DR_fun(
    x = data_meta,
    group = "LiverFatClass",
    group_names = grp_names,
    group_colors = grp_colors,    
    DRtype = "PCA",
    occ_cutoff = 0.5)
#> PCA
#> 55 samples x 659 variables
#> standard scaling of predictors
#>       R2X(cum) pre ort
#> Total    0.526  10   0
```

<img src="21-Metabolome_DimensionReduction_files/figure-html/unnamed-chunk-5-1.png" width="100%" /><img src="21-Metabolome_DimensionReduction_files/figure-html/unnamed-chunk-5-2.png" width="100%" />


结果：PCA将原始features重新组合成新的主成分，进而形成了基于主成分的矩阵。每个主成分是由原来features的线性组合而成。

+ **左上图**: 碎石图展示不同主成分的可解释度，也即是代表全体的features，评估主成分是否足够

+ **右上图**: 离群点图，它展示的各样本在投影平面和投影平面正交的距离，数值较高样本表明它们与其他样本间的差异较大 

+ **左下图**：每个样本的主成分得分（可简单理解为系数，此处只展示第一和第二主成分），该数值表示主成分能解释多少features: 该数值由features矩阵分解得到。各个样本在PC1和PC2轴的排序坐标差异，通过它们可评估样本在组成上的差异。

+ **右下图**：载荷值: 理解为构成每个主成分的features，可以认为一个主成分是由N个feature线性组合合成，每个features对该主成分的贡献程度作为loading的数值。边缘处的变量表示它们在各样本中的含量差别明显（如在某些样本中具有较大/较小的极端值等），即对排序空间的贡献较大，暗示它们可能为一些重要的代谢物。



## 偏最小二乘回归分析法（PLS-DA）

PLS-DA需要提供Y响应变量，它通过投影分别将预测变量（Y响应变量）和观测变量（自变量）投影到一个新空间，来寻找一个线性回归模型。通过建立组学数据与样本类别之间的关系模型，实现对样本类别的预测，为有监督的建模方式。

偏最小二乘（PLS）是一种基于预测变量和响应变量之间协方差的潜在变量回归方法，找到不同类别的分割信息最优的线性组合，PLS的每一个成分都赋予一个权值，即对应该成分的分离能力。



```r
PLS_res <- DR_fun(
    x = data_meta,
    group = "LiverFatClass",
    group_names = grp_names,
    group_colors = grp_colors,    
    DRtype = "PLS",
    occ_cutoff = 0.5)
#> PLS-DA
#> 55 samples x 659 variables and 1 response
#> standard scaling of predictors and response(s)
#>       R2X(cum) R2Y(cum) Q2(cum) RMSEE pre ort pR2Y  pQ2
#> Total    0.117    0.376  -0.211 0.351   2   0  0.1 0.65
```

<img src="21-Metabolome_DimensionReduction_files/figure-html/unnamed-chunk-6-1.png" width="100%" /><img src="21-Metabolome_DimensionReduction_files/figure-html/unnamed-chunk-6-2.png" width="100%" />

结果：R2X和R2Y分别表示所建模型对X和Y矩阵的解释率，Q2标示模型的预测能力，它们的值越接近于1表明模型的拟合度越好，训练集的样本越能够被准确划分到其原始归属中。

+ **左上图**: 展示了2个正交轴的R2Y和Q2Y

+ **右上图**: PLS-DA模型的R2Y和Q2Y与随机置换数据后获得的相应值进行比较 

+ **左下图**：展示了各样本在投影平面内以及正交投影面的距离，具有高值的样本标注出名称，表明它们与其它样本间的差异较大。颜色代表性别分组。

+ **右下图**：各样本在PLS-DA轴中的坐标，颜色代表不同分组。我们可以看到，相对于上文的PCA（仅通过方差特征值分解），PLS-DA在区分组间差异时更有效（带监督的偏最小二乘判别分析）。图的下方还提供了R2X、R2Y等值，用于评估模型优度。


此外，还可通过变量投影重要度（Variable Importance for the Projection，VIP）衡量各代谢物组分含量对样本分类判别的影响强度和解释能力，辅助标志代谢物的筛选。通常以**VIP > 1**作为筛选标准。


```r
vip_values <- getVipVn(PLS_res)

vip_select <- vip_values[vip_values > 1]

head(vip_select)
#> Chem_100002356 Chem_100000657 Chem_100009014 Chem_100009002 
#>       1.004056       1.070549       1.436139       1.091628 
#> Chem_100009009 Chem_100009069 
#>       2.010289       1.100467
```



## 正交偏最小二乘判别分析（OPLS-DA）

PLS-DA容易出现过拟合（overfitting）问题。所谓过拟合，即通过训练集建立了一个预测模型，它在训练集上表现出色，但通过测试集测试时却表现不佳。过拟合是机器学习中的一个常见问题，主要出现在具有比样本数量更多的变量数量的数据集的分析中。主要是因为PLS-DA是有监督的预测模型且偏向训练数据。

相比PLS-DA，叠加了正交分解的OPLS-DA能更好地避免过拟合现象，这是因为它会通过将与X和Y正交的变量删除，避免不差异的变量影响结果。OPLS从给定的自变量数据集中移除正交变量，并把这些正交变量和非正交变量区分开来。然后再对非正交变量做偏最小二乘法获得载荷矩阵也即是在Y情况下的最大方差。


+ None vs Severe

```r
DR_fun(
    x = data_meta,
    group = "LiverFatClass",
    group_names = grp_names[c(1, 4)],
    group_colors = grp_colors[c(1, 4)],
    DRtype = "OPLS",
    occ_cutoff = 0.5)
#> OPLS-DA
#> 22 samples x 659 variables and 1 response
#> standard scaling of predictors and response(s)
#>       R2X(cum) R2Y(cum) Q2(cum)  RMSEE pre ort pR2Y  pQ2
#> Total    0.244    0.988   0.622 0.0606   1   2 0.95 0.05
```

<img src="21-Metabolome_DimensionReduction_files/figure-html/unnamed-chunk-8-1.png" width="100%" /><img src="21-Metabolome_DimensionReduction_files/figure-html/unnamed-chunk-8-2.png" width="100%" />

```
#> OPLS-DA
#> 22 samples x 659 variables and 1 response
#> standard scaling of predictors and response(s)
#>       R2X(cum) R2Y(cum) Q2(cum)  RMSEE pre ort pR2Y  pQ2
#> Total    0.244    0.988   0.622 0.0606   1   2 0.95 0.05
```


结果和PLS-DA一致。

+ **左上图**: 展示了2个正交轴的R2Y和Q2Y

+ **右上图**: PLS-DA模型的R2Y和Q2Y与随机置换数据后获得的相应值进行比较 

+ **左下图**：展示了各样本在投影平面内以及正交投影面的距离，具有高值的样本标注出名称，表明它们与其它样本间的差异较大。颜色代表性别分组。

+ **右下图**：各样本在PLS-DA轴中的坐标，颜色代表不同分组。我们可以看到，相对于上文的PCA（仅通过方差特征值分解），PLS-DA在区分组间差异时更有效（带监督的偏最小二乘判别分析）。图的下方还提供了R2X、R2Y等值，用于评估模型优度。



## Session info

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
#>  date     2023-12-06
#>  pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
#> 
#> ─ Packages ───────────────────────────────────────────────
#>  package              * version   date (UTC) lib source
#>  ade4                   1.7-22    2023-02-06 [2] CRAN (R 4.1.2)
#>  ANCOMBC                1.4.0     2021-10-26 [2] Bioconductor
#>  annotate               1.72.0    2021-10-26 [2] Bioconductor
#>  AnnotationDbi          1.60.2    2023-03-10 [2] Bioconductor
#>  ape                    5.7-1     2023-03-13 [2] CRAN (R 4.1.2)
#>  Biobase              * 2.54.0    2021-10-26 [2] Bioconductor
#>  BiocGenerics         * 0.40.0    2021-10-26 [2] Bioconductor
#>  BiocParallel           1.28.3    2021-12-09 [2] Bioconductor
#>  biomformat             1.22.0    2021-10-26 [2] Bioconductor
#>  Biostrings             2.62.0    2021-10-26 [2] Bioconductor
#>  bit                    4.0.5     2022-11-15 [2] CRAN (R 4.1.2)
#>  bit64                  4.0.5     2020-08-30 [2] CRAN (R 4.1.0)
#>  bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.1.0)
#>  blob                   1.2.4     2023-03-17 [2] CRAN (R 4.1.2)
#>  bookdown               0.34      2023-05-09 [2] CRAN (R 4.1.2)
#>  bslib                  0.6.0     2023-11-21 [1] CRAN (R 4.1.3)
#>  cachem                 1.0.8     2023-05-01 [2] CRAN (R 4.1.2)
#>  callr                  3.7.3     2022-11-02 [2] CRAN (R 4.1.2)
#>  caTools                1.18.2    2021-03-28 [2] CRAN (R 4.1.0)
#>  cli                    3.6.1     2023-03-23 [2] CRAN (R 4.1.2)
#>  cluster                2.1.4     2022-08-22 [2] CRAN (R 4.1.2)
#>  codetools              0.2-19    2023-02-01 [2] CRAN (R 4.1.2)
#>  colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.1.2)
#>  crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.1.2)
#>  data.table             1.14.8    2023-02-17 [2] CRAN (R 4.1.2)
#>  DBI                    1.1.3     2022-06-18 [2] CRAN (R 4.1.2)
#>  DelayedArray           0.20.0    2021-10-26 [2] Bioconductor
#>  DESeq2                 1.34.0    2021-10-26 [2] Bioconductor
#>  devtools               2.4.5     2022-10-11 [2] CRAN (R 4.1.2)
#>  digest                 0.6.33    2023-07-07 [1] CRAN (R 4.1.3)
#>  downlit                0.4.3     2023-06-29 [2] CRAN (R 4.1.3)
#>  dplyr                * 1.1.2     2023-04-20 [2] CRAN (R 4.1.2)
#>  ellipsis               0.3.2     2021-04-29 [2] CRAN (R 4.1.0)
#>  evaluate               0.21      2023-05-05 [2] CRAN (R 4.1.2)
#>  fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.1.2)
#>  fastmap                1.1.1     2023-02-24 [2] CRAN (R 4.1.2)
#>  forcats              * 1.0.0     2023-01-29 [2] CRAN (R 4.1.2)
#>  foreach                1.5.2     2022-02-02 [2] CRAN (R 4.1.2)
#>  fs                     1.6.2     2023-04-25 [2] CRAN (R 4.1.2)
#>  genefilter             1.76.0    2021-10-26 [2] Bioconductor
#>  geneplotter            1.72.0    2021-10-26 [2] Bioconductor
#>  generics               0.1.3     2022-07-05 [2] CRAN (R 4.1.2)
#>  GenomeInfoDb         * 1.30.1    2022-01-30 [2] Bioconductor
#>  GenomeInfoDbData       1.2.7     2022-03-09 [2] Bioconductor
#>  GenomicRanges        * 1.46.1    2021-11-18 [2] Bioconductor
#>  ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.1.2)
#>  glmnet                 4.1-7     2023-03-23 [2] CRAN (R 4.1.2)
#>  glue                   1.6.2     2022-02-24 [2] CRAN (R 4.1.2)
#>  gplots                 3.1.3     2022-04-25 [2] CRAN (R 4.1.2)
#>  gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.1.2)
#>  gtools                 3.9.4     2022-11-27 [2] CRAN (R 4.1.2)
#>  highr                  0.10      2022-12-22 [2] CRAN (R 4.1.2)
#>  hms                    1.1.3     2023-03-21 [2] CRAN (R 4.1.2)
#>  htmltools              0.5.7     2023-11-03 [1] CRAN (R 4.1.3)
#>  htmlwidgets            1.6.2     2023-03-17 [2] CRAN (R 4.1.2)
#>  httpuv                 1.6.11    2023-05-11 [2] CRAN (R 4.1.3)
#>  httr                   1.4.6     2023-05-08 [2] CRAN (R 4.1.2)
#>  igraph                 1.5.0     2023-06-16 [1] CRAN (R 4.1.3)
#>  IRanges              * 2.28.0    2021-10-26 [2] Bioconductor
#>  iterators              1.0.14    2022-02-05 [2] CRAN (R 4.1.2)
#>  jquerylib              0.1.4     2021-04-26 [2] CRAN (R 4.1.0)
#>  jsonlite               1.8.7     2023-06-29 [2] CRAN (R 4.1.3)
#>  KEGGREST               1.34.0    2021-10-26 [2] Bioconductor
#>  KernSmooth             2.23-22   2023-07-10 [2] CRAN (R 4.1.3)
#>  knitr                  1.43      2023-05-25 [2] CRAN (R 4.1.3)
#>  later                  1.3.1     2023-05-02 [2] CRAN (R 4.1.2)
#>  lattice                0.21-8    2023-04-05 [2] CRAN (R 4.1.2)
#>  lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.1.2)
#>  limma                  3.50.3    2022-04-07 [2] Bioconductor
#>  locfit                 1.5-9.8   2023-06-11 [2] CRAN (R 4.1.3)
#>  lubridate            * 1.9.2     2023-02-10 [2] CRAN (R 4.1.2)
#>  magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.1.2)
#>  MASS                   7.3-60    2023-05-04 [2] CRAN (R 4.1.2)
#>  Matrix                 1.6-0     2023-07-08 [2] CRAN (R 4.1.3)
#>  MatrixGenerics       * 1.6.0     2021-10-26 [2] Bioconductor
#>  matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.1.3)
#>  memoise                2.0.1     2021-11-26 [2] CRAN (R 4.1.0)
#>  metagenomeSeq          1.36.0    2021-10-26 [2] Bioconductor
#>  mgcv                   1.8-42    2023-03-02 [2] CRAN (R 4.1.2)
#>  microbiome             1.16.0    2021-10-26 [2] Bioconductor
#>  MicrobiomeAnalysis   * 1.0.3     2023-12-02 [1] Bioconductor
#>  mime                   0.12      2021-09-28 [2] CRAN (R 4.1.0)
#>  miniUI                 0.1.1.1   2018-05-18 [2] CRAN (R 4.1.0)
#>  multtest               2.50.0    2021-10-26 [2] Bioconductor
#>  munsell                0.5.0     2018-06-12 [2] CRAN (R 4.1.0)
#>  nlme                   3.1-162   2023-01-31 [2] CRAN (R 4.1.2)
#>  nloptr                 2.0.3     2022-05-26 [2] CRAN (R 4.1.2)
#>  permute                0.9-7     2022-01-27 [2] CRAN (R 4.1.2)
#>  phyloseq               1.38.0    2021-10-26 [2] Bioconductor
#>  pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.1.2)
#>  pkgbuild               1.4.2     2023-06-26 [2] CRAN (R 4.1.3)
#>  pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.1.0)
#>  pkgload                1.3.2.1   2023-07-08 [2] CRAN (R 4.1.3)
#>  plyr                   1.8.8     2022-11-11 [2] CRAN (R 4.1.2)
#>  png                    0.1-8     2022-11-29 [2] CRAN (R 4.1.2)
#>  prettyunits            1.1.1     2020-01-24 [2] CRAN (R 4.1.0)
#>  processx               3.8.2     2023-06-30 [2] CRAN (R 4.1.3)
#>  profvis                0.3.8     2023-05-02 [2] CRAN (R 4.1.2)
#>  promises               1.2.0.1   2021-02-11 [2] CRAN (R 4.1.0)
#>  ps                     1.7.5     2023-04-18 [2] CRAN (R 4.1.2)
#>  purrr                * 1.0.1     2023-01-10 [2] CRAN (R 4.1.2)
#>  R6                     2.5.1     2021-08-19 [2] CRAN (R 4.1.0)
#>  rbibutils              2.2.13    2023-01-13 [2] CRAN (R 4.1.2)
#>  RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.1.2)
#>  Rcpp                   1.0.11    2023-07-06 [1] CRAN (R 4.1.3)
#>  RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.1.2)
#>  Rdpack                 2.4       2022-07-20 [2] CRAN (R 4.1.2)
#>  readr                * 2.1.4     2023-02-10 [2] CRAN (R 4.1.2)
#>  remotes                2.4.2     2021-11-30 [2] CRAN (R 4.1.0)
#>  reshape2               1.4.4     2020-04-09 [2] CRAN (R 4.1.0)
#>  rhdf5                  2.38.1    2022-03-10 [2] Bioconductor
#>  rhdf5filters           1.6.0     2021-10-26 [2] Bioconductor
#>  Rhdf5lib               1.16.0    2021-10-26 [2] Bioconductor
#>  rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.1.2)
#>  rmarkdown              2.23      2023-07-01 [2] CRAN (R 4.1.3)
#>  ropls                * 1.26.4    2022-01-11 [2] Bioconductor
#>  RSQLite                2.3.1     2023-04-03 [2] CRAN (R 4.1.2)
#>  rstudioapi             0.15.0    2023-07-07 [2] CRAN (R 4.1.3)
#>  Rtsne                  0.16      2022-04-17 [2] CRAN (R 4.1.2)
#>  S4Vectors            * 0.32.4    2022-03-29 [2] Bioconductor
#>  sass                   0.4.6     2023-05-03 [2] CRAN (R 4.1.2)
#>  scales                 1.2.1     2022-08-20 [2] CRAN (R 4.1.2)
#>  sessioninfo            1.2.2     2021-12-06 [2] CRAN (R 4.1.0)
#>  shape                  1.4.6     2021-05-19 [2] CRAN (R 4.1.0)
#>  shiny                  1.7.4.1   2023-07-06 [2] CRAN (R 4.1.3)
#>  stringi                1.7.12    2023-01-11 [2] CRAN (R 4.1.2)
#>  stringr              * 1.5.0     2022-12-02 [2] CRAN (R 4.1.2)
#>  SummarizedExperiment * 1.24.0    2021-10-26 [2] Bioconductor
#>  survival               3.5-5     2023-03-12 [2] CRAN (R 4.1.2)
#>  tibble               * 3.2.1     2023-03-20 [2] CRAN (R 4.1.2)
#>  tidyr                * 1.3.0     2023-01-24 [2] CRAN (R 4.1.2)
#>  tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.1.2)
#>  tidyverse            * 2.0.0     2023-02-22 [1] CRAN (R 4.1.2)
#>  timechange             0.2.0     2023-01-11 [2] CRAN (R 4.1.2)
#>  tzdb                   0.4.0     2023-05-12 [2] CRAN (R 4.1.3)
#>  urlchecker             1.0.1     2021-11-30 [2] CRAN (R 4.1.0)
#>  usethis                2.2.2     2023-07-06 [2] CRAN (R 4.1.3)
#>  utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.1.2)
#>  vctrs                  0.6.3     2023-06-14 [1] CRAN (R 4.1.3)
#>  vegan                  2.6-4     2022-10-11 [2] CRAN (R 4.1.2)
#>  withr                  2.5.0     2022-03-03 [2] CRAN (R 4.1.2)
#>  Wrench                 1.12.0    2021-10-26 [2] Bioconductor
#>  xfun                   0.40      2023-08-09 [1] CRAN (R 4.1.3)
#>  XML                    3.99-0.14 2023-03-19 [2] CRAN (R 4.1.2)
#>  xml2                   1.3.5     2023-07-06 [2] CRAN (R 4.1.3)
#>  xtable                 1.8-4     2019-04-21 [2] CRAN (R 4.1.0)
#>  XVector                0.34.0    2021-10-26 [2] Bioconductor
#>  yaml                   2.3.7     2023-01-23 [2] CRAN (R 4.1.2)
#>  zlibbioc               1.40.0    2021-10-26 [2] Bioconductor
#> 
#>  [1] /Users/zouhua/Library/R/x86_64/4.1/library
#>  [2] /Library/Frameworks/R.framework/Versions/4.1/Resources/library
#> 
#> ──────────────────────────────────────────────────────────
```


## Reference

+ [ropls: PCA, PLS(-DA) and OPLS(-DA) for multivariate analysis and feature selection of omics data](https://master.bioconductor.org/packages/release/bioc/vignettes/ropls/inst/doc/ropls-vignette.html)

