


# Other Analysis {#OtherAnalysis}


除了常见的功能分析，还有其他的功能分析方法或R包。本章节主要介绍其他功能分析的方法以及结果解析。

## FELLA: an R package to enrich metabolomics data

FELLA（[@picart2018fella]）发表于2018年，现谷歌引用67次，它是一个专门用于代谢组通路分析的R包。

> 基于前期分析得到的差异代谢物来构建基于网络的富集分析。结果包括代谢通路、模块、酶、反应及代谢物。那么除了能够提供通路列表，FELLA还能够生成输入代谢物相关的中间物质（如模块、酶、反应）。可以反映特定研究条件下代谢通路之间的交集以及靶向潜在的酶和代谢物。

它包含了以下三步：

1. **Block I**: local database从数据库抓取数据后，将其处理转存在本地；

2. **Block II**: enrichment analysis将关心的代谢物作为输入，做富集分析；

3. **Block III**: exporting results导出数据。


### 加载R包


```r
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
library(tidyverse)
library(FELLA)

# rm(list = ls())
options(stringsAsFactors = F)
options(future.globals.maxSize = 1000 * 1024^2)

grp_names <- c("None", "Mild", "Moderate", "Severe")
grp_colors <- c("#7DD06F", "#844081", "#688EC1", "#C17E73")
```


### 背景数据库生成

代谢组测试数据是来自病人血清，但也提供构建小鼠的通路背景数据库（每次下载数据库可能会发生变化，因为官网可能更新过）

+ 人的KEGG背景数据库

1. 从KEGG官网下载数据；
2. 构建背景数据库；
3. 导入内存环境.


```r
library(KEGGREST)
library(igraph)

tmpdir <- "./InputData/FELLA/hsa"

if (!file.exists("./InputData/FELLA/hsa/keggdata.graph.RData")) {
  set.seed(123)
  # 下载KEGG
  graph <- buildGraphFromKEGGREST(
    organism = "hsa", 
    filter.path = "hsa01100")
  
  tmpdir <- "./InputData/FELLA/hsa"
  unlink(tmpdir, recursive = TRUE)
  
  # 构建数据库
  buildDataFromGraph(
    keggdata.graph = graph,
    databaseDir = tmpdir,
    internalDir = FALSE,
    matrices = "none",
    normality = "diffusion",
    niter = 100)  
}

# 导入数据库进内存
fella.data <- loadKEGGdata(
  databaseDir = tmpdir,
  internalDir = FALSE,
  loadMatrix = "none")

fella.data
#> General data:
#> - KEGG graph:
#>   * Nodes:  11891 
#>   * Edges:  38085 
#>   * Density:  0.0002693728 
#>   * Categories:
#>     + pathway [349]
#>     + module [193]
#>     + enzyme [1195]
#>     + reaction [5856]
#>     + compound [4298]
#>   * Size:  5.9 Mb 
#> - KEGG names are ready.
#> -----------------------------
#> Hypergeometric test:
#> - Matrix not loaded.
#> -----------------------------
#> Heat diffusion:
#> - Matrix not loaded.
#> - RowSums are ready.
#> -----------------------------
#> PageRank:
#> - Matrix not loaded.
#> - RowSums not loaded.
```


+ 小鼠的KEGG背景数据库

1. 从KEGG官网下载数据；
2. 构建背景数据库；
3. 基因名字转换成entrez ID，酶转换成entrez ID；
4. 导入内存环境.

```R
library(KEGGREST)
library(igraph)
library(org.Mm.eg.db)

tmpdir <- "./InputData/FELLA/mmu"

if (!file.exists("./InputData/FELLA/mmu/keggdata.graph.RData")) {
  set.seed(123)
  # 下载KEGG
graph <- buildGraphFromKEGGREST(
  organism = "mmu",
  filter.path = c("01100", "01200", "01210", "01212", "01230"))
  
  tmpdir <- "./InputData/FELLA/mmu"
  unlink(tmpdir, recursive = TRUE)
  
  # 构建数据库
  buildDataFromGraph(
    keggdata.graph = graph,
    databaseDir = tmpdir,
    internalDir = FALSE,
    matrices = "none",
    normality = "diffusion",
    niter = 100)  
  
  alias2entrez <- as.list(org.Mm.eg.db::org.Mm.egSYMBOL2EG)
  entrez2ec <- KEGGREST::keggLink("enzyme", "mmu")
  entrez2path <- KEGGREST::keggLink("pathway", "mmu")

}

# 导入数据库进内存
fella.data <- loadKEGGdata(
  databaseDir = tmpdir,
  internalDir = FALSE,
  loadMatrix = "none")

fella.data
```


### 导入数据

对数据[OmicsDataSet-Zeybel et al. - 2022.xlsx](https://github.com/HuaZou/DraftNotes/blob/main/InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx)处理后生成的，可参考数据预处理等章节。


> ```R
> write.table(final_res, "./InputData/result/DA/Metabolites_FC_VIP_ttest.tsv", 
            row.names = F, quote = F, sep = "\t", fileEncoding = "UTF-8")
> ```


```r
datSignif <- data.table::fread("./InputData/result/DA/Metabolites_FC_VIP_ttest.tsv")

# DT::datatable(datSignif)

head(datSignif)
#>                                         FeatureID
#> 1: ceramide (d18:1/20:0, d16:1/22:0, d20:1/18:0)*
#> 2:                 cysteine-glutathione disulfide
#> 3:          1-palmitoyl-2-oleoyl-GPI (16:0/18:1)*
#> 4:                                         serine
#> 5:           1-stearoyl-2-oleoyl-GPI (18:0/18:1)*
#> 6:     palmitoyl-oleoyl-glycerol (16:0/18:1) [2]*
#>            Block2                Block FoldChange
#> 1: None vs Severe 10_None vs 12_Severe  0.6444244
#> 2: None vs Severe 10_None vs 12_Severe  1.7109000
#> 3: None vs Severe 10_None vs 12_Severe  0.5199556
#> 4: None vs Severe 10_None vs 12_Severe  1.2218596
#> 5: None vs Severe 10_None vs 12_Severe  0.5667863
#> 6: None vs Severe 10_None vs 12_Severe  0.5638085
#>    Log2FoldChange      VIP    CorPvalue Statistic
#> 1:     -0.6339170 2.680471 7.910566e-05 -4.875928
#> 2:      0.7747554 2.615284 1.427798e-04  4.879864
#> 3:     -0.9435396 2.535490 2.775095e-04 -4.333042
#> 4:      0.2890785 2.490737 3.929470e-04  4.314434
#> 5:     -0.8191231 2.489832 3.956533e-04 -4.144802
#> 6:     -0.8267228 2.381443 8.605197e-04 -3.781278
#>          Pvalue AdjustedPvalue Mean Abundance (All)
#> 1: 0.0001197367     0.03975257              3841099
#> 2: 0.0001073822     0.03975257              1246453
#> 3: 0.0004066026     0.06749603              2243154
#> 4: 0.0003379029     0.06749603             63358904
#> 5: 0.0006948992     0.09228261              1817773
#> 6: 0.0017211306     0.16326153              1192929
#>    Mean Abundance None Mean Abundance Severe  metabolitesID
#> 1:           2952496.1             4581602.1 Chem_100015755
#> 2:           1611743.8              942044.4 Chem_100001437
#> 3:           1491869.7             2869225.1 Chem_100009066
#> 4:          70323857.2            57554776.3       Chem_503
#> 5:           1282914.5             2263488.8 Chem_100009181
#> 6:            838913.8             1487941.0 Chem_100010917
#>                                       BIOCHEMICAL
#> 1: ceramide (d18:1/20:0, d16:1/22:0, d20:1/18:0)*
#> 2:                 cysteine-glutathione disulfide
#> 3:          1-palmitoyl-2-oleoyl-GPI (16:0/18:1)*
#> 4:                                         serine
#> 5:           1-stearoyl-2-oleoyl-GPI (18:0/18:1)*
#> 6:     palmitoyl-oleoyl-glycerol (16:0/18:1) [2]*
#>    SUPER.PATHWAY                              SUB.PATHWAY
#> 1:         Lipid                                Ceramides
#> 2:    Amino Acid                   Glutathione Metabolism
#> 3:         Lipid                Phosphatidylinositol (PI)
#> 4:    Amino Acid Glycine, Serine and Threonine Metabolism
#> 5:         Lipid                Phosphatidylinositol (PI)
#> 6:         Lipid                           Diacylglycerol
#>    COMPID        PLATFORM CHEMICALID   RI     MASS  PUBCHEM
#> 1:  57440  LC/MS Pos Late  100015755 3920 594.5820     <NA>
#> 2:  35159 LC/MS Pos Early  100001437 2465 427.0952  3080690
#> 3:  52669  LC/MS Pos Late  100009066 3140 854.5753 71296232
#> 4:   1648 LC/MS Pos Early        503 1239 106.0499     5951
#> 5:  52726  LC/MS Pos Late  100009181 3711 882.6066     <NA>
#> 6:  54942  LC/MS Pos Late  100010917 3695 612.5562  5282283
#>           CAS   KEGG SampleIDHMDBID
#> 1:       <NA>   <NA>           <NA>
#> 2: 13081-14-6 R00900    HMDB0000656
#> 3:       <NA>   <NA>    HMDB0009783
#> 4:    56-45-1 C00065    HMDB0000187
#> 5:       <NA>   <NA>           <NA>
#> 6:       <NA> C13861    HMDB0007102
```



### 准备输入代谢物

代谢物的ID要是KEGG ID，需要注意⚠️。随机挑选5个代谢物用于分析。


```r
set.seed(123)

datSignif$KEGG <- gsub(",\\S+", "", datSignif$KEGG)

datSignif_KEGG <- datSignif %>%
  dplyr::filter(!is.na(KEGG)) %>%
  dplyr::filter(SUPER.PATHWAY == "Amino Acid") %>%
  dplyr::select(BIOCHEMICAL, KEGG) 

target_metabolites <- datSignif_KEGG[sample(1:nrow(datSignif_KEGG), 5), ,]

head(target_metabolites)
#>           BIOCHEMICAL   KEGG
#> 1:         isoleucine C16424
#> 2:    N-acetylalanine C02847
#> 3:          aspartate C00049
#> 4: O-acetylhomoserine C01077
#> 5:           cysteine C00097
```


### 富集分析

富集分析的方法有三种

  - 超几何检验
  
  - Diffusion（有意义子网络）
  
  - PageRank（和Diffusion类似，对网络进行排序）
  
统计分析：对Diffusion和PageRank提供了两种统计方法

  - Normal approximation(approx = "normality")，基于无效假设的分析的期望值和协方差矩阵的z-score计算得到得分值
  
  - Monte Carlo trials(approx = "simulation")，随机变量的蒙特卡罗实验计算得分值


+ 可通过method选择不同富集方法，本次运行选择*diffusion*


```r
myAnalysis <- enrich(
    compounds = target_metabolites$KEGG, 
    method = "diffusion", # listMethods()
    approx = "normality", 
    data = fella.data)

show(myAnalysis)
#> Compounds in the input: 4
#> [1] "C00049" "C00097" "C01077" "C16424"
#> Background compounds: all available compounds (default)
#> -----------------------------
#> Hypergeometric test: not performed
#> -----------------------------
#> Heat diffusion: ready.
#> P-scores under 0.05:  166
#> -----------------------------
#> PageRank: not performed
```

结果：展示了diffusion方法下富集的结果, 有104个节点。


+ 可视化结果


```r
plot(
    x = myAnalysis, 
    method = "diffusion", 
    main = "diffusion analysis in FELLA", 
    threshold = 0.1, 
    data = fella.data,
    nlimit = 100)
```

<img src="25-Metabolome-OtherAnalysis_files/figure-html/unnamed-chunk-7-1.png" width="100%" />



+ 输出富集分析结果表格


```r
myTable <- generateResultsTable(
    object = myAnalysis, 
    method = "diffusion", 
    threshold = 0.1, 
    data = fella.data)

knitr::kable(head(myTable, 10), format = "html")
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> KEGG.id </th>
   <th style="text-align:left;"> Entry.type </th>
   <th style="text-align:left;"> KEGG.name </th>
   <th style="text-align:right;"> p.score </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> hsa00250 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> Alanine, aspartate and glutamate metabolism -... </td>
   <td style="text-align:right;"> 0.0000029 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hsa00270 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> Cysteine and methionine metabolism - Homo sap... </td>
   <td style="text-align:right;"> 0.0000010 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hsa01230 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> Biosynthesis of amino acids - Homo sapiens (h... </td>
   <td style="text-align:right;"> 0.0972264 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hsa04614 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> Renin-angiotensin system - Homo sapiens (huma... </td>
   <td style="text-align:right;"> 0.0216791 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hsa04621 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> NOD-like receptor signaling pathway - Homo sa... </td>
   <td style="text-align:right;"> 0.0556380 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hsa04974 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> Protein digestion and absorption - Homo sapie... </td>
   <td style="text-align:right;"> 0.0272929 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hsa05230 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> Central carbon metabolism in cancer - Homo sa... </td>
   <td style="text-align:right;"> 0.0000095 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> M00017 </td>
   <td style="text-align:left;"> module </td>
   <td style="text-align:left;"> Methionine biosynthesis, aspartate =&gt; homoser... </td>
   <td style="text-align:right;"> 0.0018276 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> M00027 </td>
   <td style="text-align:left;"> module </td>
   <td style="text-align:left;"> GABA (gamma-Aminobutyrate) shunt </td>
   <td style="text-align:right;"> 0.0769332 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> M00029 </td>
   <td style="text-align:left;"> module </td>
   <td style="text-align:left;"> Urea cycle </td>
   <td style="text-align:right;"> 0.0704110 </td>
  </tr>
</tbody>
</table>





### 结果解析

+ 筛选的5个代谢物富集在(**diffusion**)"hsa00250"和"hsa00270"等通路，并且这些通路大部分和氨基酸代谢相关；

+ 除了代谢通路外，还有代谢模块等其他更为具体的通路组成，比如酶和反应等；

+ 相比传统的富集分析，FELLA能将代谢通路各个层级混合在一起做成网络分析是其特点，比如**p53 signaling pathway - Homo sapiens (human)**相关的酶是**ribonucleoside-diphosphate reductase**，该酶又和反应**5-fluorodeoxyuridine-diphosphate**相关。


## Pathview: 代谢物数据可视化KEGG通路图

> Pathview（Pathway based data integration and visualization，https://pathview.uncc.edu/）是一个用于KEGG通路可视化的工具集，能够将多种生物的基因或代谢物映射到该物种的KEGG通路图上，例如在转录组、蛋白组或代谢组中展示差异表达的基因、蛋白或代谢物等。

本教程是基于代谢组数据，使用代谢组KEGGID (如`C00064`）等。

### 安装pathview包

```R
if (!requireNamespace("pathview", quietly=TRUE)) {
  BiocManager::install('pathview')
}
```

### 准备输入数据

筛选包含Log2FoldChange和KEGG的代谢物输入数据（pathview要求输入数据包含差异倍数以及代谢物的KEGGID）。


```r
set.seed(123)

datSignif <- data.table::fread("./InputData/result/DA/Metabolites_FC_VIP_ttest.tsv")

datSignif$KEGG <- gsub(",\\S+", "", datSignif$KEGG)

datSignif_KEGG <- datSignif %>%
  dplyr::filter(!is.na(KEGG)) %>%
  dplyr::filter(SUPER.PATHWAY == "Amino Acid") %>%
  dplyr::select(BIOCHEMICAL, KEGG, Log2FoldChange) 

# head(datSignif_KEGG)

compound_data <- datSignif_KEGG$Log2FoldChange
names(compound_data) <- datSignif_KEGG$KEGG

compound_data[1:6]
#>     R00900     C00065     C00135     C05568     C01188 
#>  0.7747554  0.2890785  0.1947889  0.3253860 -0.4819788 
#>     C00719 
#>  0.2584803
```


### 运行pathview

根据*myTable*可以看到富集在00250和00310等通路，最后选择00250通路展示。


```r
library(pathview)

pl <- pathview(
  cpd.data = compound_data,
  cpd.idtype   = "kegg",
  species      = "hsa", 
  kegg.native  = TRUE,
  pathway.id   = "hsa00250",
  out.suffix   = "compound",
  limit        = list(gene=1, cpd=max(abs(datSignif_KEGG$Log2FoldChange))))
```


```
#> [1] TRUE
```

<div class="figure" style="text-align: center">
<img src="./InputData/figures/hsa00250.compound.png" alt="KEGG of 00250 pathway" width="100%" height="100%" />
<p class="caption">(\#fig:unnamed-chunk-11)KEGG of 00250 pathway</p>
</div>

结果：

+ 蓝色和黄色表示富集方向，能看到3个颜色明显的代谢物在该通路中发挥的作用。



## aPEAR: an R package for autonomous visualisation of pathway enrichment networks

> The interpretation of pathway enrichment analysis results is frequently complicated by an overwhelming and redundant list of significantly affected pathways. Here, we present an R package aPEAR (Advanced Pathway Enrichment Analysis Representation) which leverages similarities between the pathway gene sets and represents them as a network of interconnected clusters. Each cluster is assigned a meaningful name that highlights the main biological themes in the experiment. Our approach enables an automated and objective overview of the data without manual and time-consuming parameter tweaking.

aPEAR利用通路间的相似性将通路进行聚类，方便解释通路富集分析的结果。aPEAR提供一个enrichmentNetwork函数，将结果可视化成网络。其原理：


> + The pairwise similarity between all pathway gene sets is evaluated using the Jaccard index (default), cosine similarity, or correlation similarity metrics.
> 
> + The similarity matrix is then used to detect clusters of redundant pathways using Markov (default) (Van Dongen 2008), hierarchical, or spectral (John et al. 2020) clustering algorithms.


+ 所有通路的基因计算成对相似性 (默认使用Jaccard距离)。

+ 将冗余的通路聚类在一起。

### 安装aPEAR包

```R
if (!requireNamespace("aPEAR", quietly=TRUE)) {
  BiocManager::install('aPEAR')
}
```

### 准备输入数据

筛选包含Log2FoldChange和KEGG的代谢物输入数据。


```r
set.seed(123)

datSignif <- data.table::fread("./InputData/result/DA/Metabolites_FC_VIP_ttest.tsv")
datSignif$KEGG <- gsub(",\\S+", "", datSignif$KEGG)
datSignif_KEGG <- datSignif %>%
  dplyr::filter(!is.na(KEGG)) %>%
  dplyr::filter(SUPER.PATHWAY == "Amino Acid") %>%
  dplyr::select(BIOCHEMICAL, KEGG, Log2FoldChange) 


hsa_kegg_compound <- read.csv("./InputData/result/KEGG/KEGG_COMPOUND_PATHWAY_hsa.csv") 
ref_cln <- hsa_kegg_compound %>%
    dplyr::select(PATHWAY_MAP, COMPOUND) %>%
    dplyr::rename(Pathway = PATHWAY_MAP)
```


### 富集分析

采用ORA富集分析方法


```r
ORA_fit <- clusterProfiler::enricher(
  gene = datSignif_KEGG$KEGG,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = ref_cln)
```


### aPEAR网络图

aPEAR画网络图


```r
library(aPEAR)

enrichmentNetwork(
  ORA_fit@result,
  fontSize = 3,
  outerCutoff = 0.5,
  drawEllipses = TRUE,
  repelLabels = TRUE)
```

<img src="25-Metabolome-OtherAnalysis_files/figure-html/unnamed-chunk-14-1.png" width="100%" />


结果: 

+ 节点表示显著通路，边表示相关性，颜色表示标准化后的富集得分；

+ 每个簇分类了一个具有生物学意义的名称；

+ NES值表示聚类簇通路的重要性，NES越高通路越重要。



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
#>  date     2023-12-05
#>  pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
#> 
#> ─ Packages ───────────────────────────────────────────────
#>  package          * version    date (UTC) lib source
#>  AnnotationDbi      1.60.2     2023-03-10 [2] Bioconductor
#>  ape                5.7-1      2023-03-13 [2] CRAN (R 4.1.2)
#>  aPEAR            * 1.0.0      2023-06-12 [1] CRAN (R 4.1.3)
#>  aplot              0.1.10     2023-03-08 [2] CRAN (R 4.1.2)
#>  arules             1.7-6      2023-03-23 [1] CRAN (R 4.1.2)
#>  bayesbio           1.0.0      2016-05-24 [1] CRAN (R 4.1.0)
#>  Biobase            2.54.0     2021-10-26 [2] Bioconductor
#>  BiocGenerics       0.40.0     2021-10-26 [2] Bioconductor
#>  BiocParallel       1.28.3     2021-12-09 [2] Bioconductor
#>  Biostrings         2.62.0     2021-10-26 [2] Bioconductor
#>  bit                4.0.5      2022-11-15 [2] CRAN (R 4.1.2)
#>  bit64              4.0.5      2020-08-30 [2] CRAN (R 4.1.0)
#>  bitops             1.0-7      2021-04-24 [2] CRAN (R 4.1.0)
#>  blob               1.2.4      2023-03-17 [2] CRAN (R 4.1.2)
#>  bookdown           0.34       2023-05-09 [2] CRAN (R 4.1.2)
#>  bslib              0.6.0      2023-11-21 [1] CRAN (R 4.1.3)
#>  cachem             1.0.8      2023-05-01 [2] CRAN (R 4.1.2)
#>  callr              3.7.3      2022-11-02 [2] CRAN (R 4.1.2)
#>  cli                3.6.1      2023-03-23 [2] CRAN (R 4.1.2)
#>  clusterProfiler    4.2.2      2022-01-13 [2] Bioconductor
#>  colorspace         2.1-0      2023-01-23 [2] CRAN (R 4.1.2)
#>  crayon             1.5.2      2022-09-29 [2] CRAN (R 4.1.2)
#>  data.table         1.14.8     2023-02-17 [2] CRAN (R 4.1.2)
#>  DBI                1.1.3      2022-06-18 [2] CRAN (R 4.1.2)
#>  devtools           2.4.5      2022-10-11 [2] CRAN (R 4.1.2)
#>  digest             0.6.33     2023-07-07 [1] CRAN (R 4.1.3)
#>  DO.db              2.9        2022-04-11 [2] Bioconductor
#>  DOSE               3.20.1     2021-11-18 [2] Bioconductor
#>  downlit            0.4.3      2023-06-29 [2] CRAN (R 4.1.3)
#>  downloader         0.4        2015-07-09 [2] CRAN (R 4.1.0)
#>  dplyr            * 1.1.2      2023-04-20 [2] CRAN (R 4.1.2)
#>  ellipsis           0.3.2      2021-04-29 [2] CRAN (R 4.1.0)
#>  enrichplot         1.14.2     2022-02-24 [2] Bioconductor
#>  evaluate           0.21       2023-05-05 [2] CRAN (R 4.1.2)
#>  expm               0.999-7    2023-01-09 [2] CRAN (R 4.1.2)
#>  fansi              1.0.4      2023-01-22 [2] CRAN (R 4.1.2)
#>  farver             2.1.1      2022-07-06 [2] CRAN (R 4.1.2)
#>  fastmap            1.1.1      2023-02-24 [2] CRAN (R 4.1.2)
#>  fastmatch          1.1-3      2021-07-23 [2] CRAN (R 4.1.0)
#>  FELLA            * 1.14.0     2021-10-26 [1] Bioconductor
#>  fgsea              1.20.0     2021-10-26 [2] Bioconductor
#>  forcats          * 1.0.0      2023-01-29 [2] CRAN (R 4.1.2)
#>  fs                 1.6.2      2023-04-25 [2] CRAN (R 4.1.2)
#>  generics           0.1.3      2022-07-05 [2] CRAN (R 4.1.2)
#>  GenomeInfoDb       1.30.1     2022-01-30 [2] Bioconductor
#>  GenomeInfoDbData   1.2.7      2022-03-09 [2] Bioconductor
#>  ggforce            0.4.1      2022-10-04 [2] CRAN (R 4.1.2)
#>  ggfun              0.1.1      2023-06-24 [2] CRAN (R 4.1.3)
#>  ggplot2          * 3.4.2      2023-04-03 [2] CRAN (R 4.1.2)
#>  ggplotify          0.1.1      2023-06-27 [2] CRAN (R 4.1.3)
#>  ggraph             2.1.0.9000 2023-07-11 [1] Github (thomasp85/ggraph@febab71)
#>  ggrepel            0.9.3      2023-02-03 [2] CRAN (R 4.1.2)
#>  ggtree             3.2.1      2021-11-16 [2] Bioconductor
#>  glue               1.6.2      2022-02-24 [2] CRAN (R 4.1.2)
#>  GO.db              3.14.0     2022-04-11 [2] Bioconductor
#>  GOSemSim           2.20.0     2021-10-26 [2] Bioconductor
#>  graph              1.72.0     2021-10-26 [2] Bioconductor
#>  graphlayouts       1.0.0      2023-05-01 [2] CRAN (R 4.1.2)
#>  gridExtra          2.3        2017-09-09 [2] CRAN (R 4.1.0)
#>  gridGraphics       0.5-1      2020-12-13 [2] CRAN (R 4.1.0)
#>  gtable             0.3.3      2023-03-21 [2] CRAN (R 4.1.2)
#>  highr              0.10       2022-12-22 [2] CRAN (R 4.1.2)
#>  hms                1.1.3      2023-03-21 [2] CRAN (R 4.1.2)
#>  htmltools          0.5.7      2023-11-03 [1] CRAN (R 4.1.3)
#>  htmlwidgets        1.6.2      2023-03-17 [2] CRAN (R 4.1.2)
#>  httpuv             1.6.11     2023-05-11 [2] CRAN (R 4.1.3)
#>  httr               1.4.6      2023-05-08 [2] CRAN (R 4.1.2)
#>  igraph           * 1.5.0      2023-06-16 [1] CRAN (R 4.1.3)
#>  IRanges            2.28.0     2021-10-26 [2] Bioconductor
#>  jquerylib          0.1.4      2021-04-26 [2] CRAN (R 4.1.0)
#>  jsonlite           1.8.7      2023-06-29 [2] CRAN (R 4.1.3)
#>  KEGGgraph          1.54.0     2021-10-26 [2] Bioconductor
#>  KEGGREST         * 1.34.0     2021-10-26 [2] Bioconductor
#>  knitr              1.43       2023-05-25 [2] CRAN (R 4.1.3)
#>  labeling           0.4.2      2020-10-20 [2] CRAN (R 4.1.0)
#>  later              1.3.1      2023-05-02 [2] CRAN (R 4.1.2)
#>  lattice            0.21-8     2023-04-05 [2] CRAN (R 4.1.2)
#>  lazyeval           0.2.2      2019-03-15 [2] CRAN (R 4.1.0)
#>  lifecycle          1.0.3      2022-10-07 [2] CRAN (R 4.1.2)
#>  lsa                0.73.3     2022-05-09 [1] CRAN (R 4.1.2)
#>  lubridate        * 1.9.2      2023-02-10 [2] CRAN (R 4.1.2)
#>  magrittr           2.0.3      2022-03-30 [2] CRAN (R 4.1.2)
#>  MASS               7.3-60     2023-05-04 [2] CRAN (R 4.1.2)
#>  Matrix             1.6-0      2023-07-08 [2] CRAN (R 4.1.3)
#>  MCL                1.0        2015-03-11 [1] CRAN (R 4.1.0)
#>  memoise            2.0.1      2021-11-26 [2] CRAN (R 4.1.0)
#>  mime               0.12       2021-09-28 [2] CRAN (R 4.1.0)
#>  miniUI             0.1.1.1    2018-05-18 [2] CRAN (R 4.1.0)
#>  munsell            0.5.0      2018-06-12 [2] CRAN (R 4.1.0)
#>  nlme               3.1-162    2023-01-31 [2] CRAN (R 4.1.2)
#>  org.Hs.eg.db       3.16.0     2023-03-22 [2] Bioconductor
#>  patchwork          1.1.2      2022-08-19 [2] CRAN (R 4.1.2)
#>  pathview         * 1.34.0     2021-10-26 [1] Bioconductor
#>  pillar             1.9.0      2023-03-22 [2] CRAN (R 4.1.2)
#>  pkgbuild           1.4.2      2023-06-26 [2] CRAN (R 4.1.3)
#>  pkgconfig          2.0.3      2019-09-22 [2] CRAN (R 4.1.0)
#>  pkgload            1.3.2.1    2023-07-08 [2] CRAN (R 4.1.3)
#>  plyr               1.8.8      2022-11-11 [2] CRAN (R 4.1.2)
#>  png                0.1-8      2022-11-29 [2] CRAN (R 4.1.2)
#>  polyclip           1.10-4     2022-10-20 [2] CRAN (R 4.1.2)
#>  prettyunits        1.1.1      2020-01-24 [2] CRAN (R 4.1.0)
#>  processx           3.8.2      2023-06-30 [2] CRAN (R 4.1.3)
#>  profvis            0.3.8      2023-05-02 [2] CRAN (R 4.1.2)
#>  promises           1.2.0.1    2021-02-11 [2] CRAN (R 4.1.0)
#>  ps                 1.7.5      2023-04-18 [2] CRAN (R 4.1.2)
#>  purrr            * 1.0.1      2023-01-10 [2] CRAN (R 4.1.2)
#>  qvalue             2.26.0     2021-10-26 [2] Bioconductor
#>  R6                 2.5.1      2021-08-19 [2] CRAN (R 4.1.0)
#>  RColorBrewer       1.1-3      2022-04-03 [2] CRAN (R 4.1.2)
#>  Rcpp               1.0.11     2023-07-06 [1] CRAN (R 4.1.3)
#>  RCurl              1.98-1.12  2023-03-27 [2] CRAN (R 4.1.2)
#>  readr            * 2.1.4      2023-02-10 [2] CRAN (R 4.1.2)
#>  remotes            2.4.2      2021-11-30 [2] CRAN (R 4.1.0)
#>  reshape2           1.4.4      2020-04-09 [2] CRAN (R 4.1.0)
#>  Rgraphviz          2.38.0     2021-10-26 [2] Bioconductor
#>  rlang              1.1.1      2023-04-28 [1] CRAN (R 4.1.2)
#>  rmarkdown          2.23       2023-07-01 [2] CRAN (R 4.1.3)
#>  RSQLite            2.3.1      2023-04-03 [2] CRAN (R 4.1.2)
#>  rstudioapi         0.15.0     2023-07-07 [2] CRAN (R 4.1.3)
#>  S4Vectors          0.32.4     2022-03-29 [2] Bioconductor
#>  sass               0.4.6      2023-05-03 [2] CRAN (R 4.1.2)
#>  scales             1.2.1      2022-08-20 [2] CRAN (R 4.1.2)
#>  scatterpie         0.2.1      2023-06-07 [2] CRAN (R 4.1.3)
#>  sessioninfo        1.2.2      2021-12-06 [2] CRAN (R 4.1.0)
#>  shadowtext         0.1.2      2022-04-22 [2] CRAN (R 4.1.2)
#>  shiny              1.7.4.1    2023-07-06 [2] CRAN (R 4.1.3)
#>  SnowballC          0.7.1      2023-04-25 [2] CRAN (R 4.1.2)
#>  stringi            1.7.12     2023-01-11 [2] CRAN (R 4.1.2)
#>  stringr          * 1.5.0      2022-12-02 [2] CRAN (R 4.1.2)
#>  tibble           * 3.2.1      2023-03-20 [2] CRAN (R 4.1.2)
#>  tidygraph          1.2.3      2023-02-01 [2] CRAN (R 4.1.2)
#>  tidyr            * 1.3.0      2023-01-24 [2] CRAN (R 4.1.2)
#>  tidyselect         1.2.0      2022-10-10 [2] CRAN (R 4.1.2)
#>  tidytree           0.4.2      2022-12-18 [2] CRAN (R 4.1.2)
#>  tidyverse        * 2.0.0      2023-02-22 [1] CRAN (R 4.1.2)
#>  timechange         0.2.0      2023-01-11 [2] CRAN (R 4.1.2)
#>  treeio             1.18.1     2021-11-14 [2] Bioconductor
#>  tweenr             2.0.2      2022-09-06 [2] CRAN (R 4.1.2)
#>  tzdb               0.4.0      2023-05-12 [2] CRAN (R 4.1.3)
#>  urlchecker         1.0.1      2021-11-30 [2] CRAN (R 4.1.0)
#>  usethis            2.2.2      2023-07-06 [2] CRAN (R 4.1.3)
#>  utf8               1.2.3      2023-01-31 [2] CRAN (R 4.1.2)
#>  vctrs              0.6.3      2023-06-14 [1] CRAN (R 4.1.3)
#>  viridis            0.6.3      2023-05-03 [2] CRAN (R 4.1.2)
#>  viridisLite        0.4.2      2023-05-02 [2] CRAN (R 4.1.2)
#>  withr              2.5.0      2022-03-03 [2] CRAN (R 4.1.2)
#>  xfun               0.40       2023-08-09 [1] CRAN (R 4.1.3)
#>  XML                3.99-0.14  2023-03-19 [2] CRAN (R 4.1.2)
#>  xml2               1.3.5      2023-07-06 [2] CRAN (R 4.1.3)
#>  xtable             1.8-4      2019-04-21 [2] CRAN (R 4.1.0)
#>  XVector            0.34.0     2021-10-26 [2] Bioconductor
#>  yaml               2.3.7      2023-01-23 [2] CRAN (R 4.1.2)
#>  yulab.utils        0.0.6      2022-12-20 [2] CRAN (R 4.1.2)
#>  zlibbioc           1.40.0     2021-10-26 [2] Bioconductor
#> 
#>  [1] /Users/zouhua/Library/R/x86_64/4.1/library
#>  [2] /Library/Frameworks/R.framework/Versions/4.1/Resources/library
#> 
#> ──────────────────────────────────────────────────────────
```


## Reference

+ [FELLA github](https://github.com/b2slab/FELLA)

+ [The MetaRbolomics book](https://rformassspectrometry.github.io/metaRbolomics-book/)
