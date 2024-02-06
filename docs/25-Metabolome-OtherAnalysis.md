


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
#>                                            <char>
#> 1: ceramide (d18:1/20:0, d16:1/22:0, d20:1/18:0)*
#> 2:                 cysteine-glutathione disulfide
#> 3:                                         serine
#> 4:          1-palmitoyl-2-oleoyl-GPI (16:0/18:1)*
#> 5:           1-stearoyl-2-oleoyl-GPI (18:0/18:1)*
#> 6:     palmitoyl-oleoyl-glycerol (16:0/18:1) [2]*
#>            Block2                Block FoldChange
#>            <char>               <char>      <num>
#> 1: None vs Severe 10_None vs 12_Severe  0.6444244
#> 2: None vs Severe 10_None vs 12_Severe  1.7109000
#> 3: None vs Severe 10_None vs 12_Severe  1.2218596
#> 4: None vs Severe 10_None vs 12_Severe  0.5199556
#> 5: None vs Severe 10_None vs 12_Severe  0.5667863
#> 6: None vs Severe 10_None vs 12_Severe  0.5638085
#>    Log2FoldChange      VIP    CorPvalue Statistic
#>             <num>    <num>        <num>     <num>
#> 1:     -0.6339170 2.713879 6.487126e-05 -4.923988
#> 2:      0.7747554 2.653166 1.139027e-04  4.989637
#> 3:      0.2890785 2.531054 3.153670e-04  4.409792
#> 4:     -0.9435396 2.539496 2.952154e-04 -4.294439
#> 5:     -0.8191231 2.488347 4.365061e-04 -4.098726
#> 6:     -0.8267228 2.398322 8.276438e-04 -3.781813
#>          Pvalue AdjustedPvalue Mean Abundance (All)
#>           <num>          <num>                <num>
#> 1: 1.234968e-04     0.03958071              3841099
#> 2: 8.636192e-05     0.03958071              1246453
#> 3: 2.705095e-04     0.05779885             63358904
#> 4: 4.600563e-04     0.07372402              2243154
#> 5: 7.802676e-04     0.10003031              1817773
#> 6: 1.814397e-03     0.16614697              1192929
#>    Mean Abundance None Mean Abundance Severe  metabolitesID
#>                  <num>                 <num>         <char>
#> 1:           2952496.1             4581602.1 Chem_100015755
#> 2:           1611743.8              942044.4 Chem_100001437
#> 3:          70323857.2            57554776.3       Chem_503
#> 4:           1491869.7             2869225.1 Chem_100009066
#> 5:           1282914.5             2263488.8 Chem_100009181
#> 6:            838913.8             1487941.0 Chem_100010917
#>                                       BIOCHEMICAL
#>                                            <char>
#> 1: ceramide (d18:1/20:0, d16:1/22:0, d20:1/18:0)*
#> 2:                 cysteine-glutathione disulfide
#> 3:                                         serine
#> 4:          1-palmitoyl-2-oleoyl-GPI (16:0/18:1)*
#> 5:           1-stearoyl-2-oleoyl-GPI (18:0/18:1)*
#> 6:     palmitoyl-oleoyl-glycerol (16:0/18:1) [2]*
#>    SUPER.PATHWAY                              SUB.PATHWAY
#>           <char>                                   <char>
#> 1:         Lipid                                Ceramides
#> 2:    Amino Acid                   Glutathione Metabolism
#> 3:    Amino Acid Glycine, Serine and Threonine Metabolism
#> 4:         Lipid                Phosphatidylinositol (PI)
#> 5:         Lipid                Phosphatidylinositol (PI)
#> 6:         Lipid                           Diacylglycerol
#>    COMPID        PLATFORM CHEMICALID    RI     MASS
#>     <int>          <char>      <int> <num>    <num>
#> 1:  57440  LC/MS Pos Late  100015755  3920 594.5820
#> 2:  35159 LC/MS Pos Early  100001437  2465 427.0952
#> 3:   1648 LC/MS Pos Early        503  1239 106.0499
#> 4:  52669  LC/MS Pos Late  100009066  3140 854.5753
#> 5:  52726  LC/MS Pos Late  100009181  3711 882.6066
#> 6:  54942  LC/MS Pos Late  100010917  3695 612.5562
#>     PUBCHEM        CAS   KEGG SampleIDHMDBID
#>      <char>     <char> <char>         <char>
#> 1:     <NA>       <NA>   <NA>           <NA>
#> 2:  3080690 13081-14-6 R00900    HMDB0000656
#> 3:     5951    56-45-1 C00065    HMDB0000187
#> 4: 71296232       <NA>   <NA>    HMDB0009783
#> 5:     <NA>       <NA>   <NA>           <NA>
#> 6:  5282283       <NA> C13861    HMDB0007102
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
#>         BIOCHEMICAL   KEGG
#>              <char> <char>
#> 1:        sarcosine C00213
#> 2:         creatine C00300
#> 3:       methionine C00073
#> 4: cysteinylglycine C01419
#> 5:          cystine C00491
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
#> Compounds in the input: 5
#> [1] "C00213" "C00300" "C00073" "C01419" "C00491"
#> Background compounds: all available compounds (default)
#> -----------------------------
#> Hypergeometric test: not performed
#> -----------------------------
#> Heat diffusion: ready.
#> P-scores under 0.05:  181
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
   <td style="text-align:left;"> hsa00260 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> Glycine, serine and threonine metabolism - Ho... </td>
   <td style="text-align:right;"> 0.0000010 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hsa00270 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> Cysteine and methionine metabolism - Homo sap... </td>
   <td style="text-align:right;"> 0.0000053 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hsa00330 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> Arginine and proline metabolism - Homo sapien... </td>
   <td style="text-align:right;"> 0.0000010 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hsa00430 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> Taurine and hypotaurine metabolism - Homo sap... </td>
   <td style="text-align:right;"> 0.0167926 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hsa00480 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> Glutathione metabolism - Homo sapiens (human) </td>
   <td style="text-align:right;"> 0.0000395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hsa01523 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> Antifolate resistance - Homo sapiens (human) </td>
   <td style="text-align:right;"> 0.0216049 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hsa02010 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> ABC transporters - Homo sapiens (human) </td>
   <td style="text-align:right;"> 0.0261448 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hsa04216 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> Ferroptosis - Homo sapiens (human) </td>
   <td style="text-align:right;"> 0.0000188 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hsa04614 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> Renin-angiotensin system - Homo sapiens (huma... </td>
   <td style="text-align:right;"> 0.0000010 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hsa04640 </td>
   <td style="text-align:left;"> pathway </td>
   <td style="text-align:left;"> Hematopoietic cell lineage - Homo sapiens (hu... </td>
   <td style="text-align:right;"> 0.0000010 </td>
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
#>  version  R version 4.3.1 (2023-06-16)
#>  os       macOS Monterey 12.2.1
#>  system   x86_64, darwin20
#>  ui       X11
#>  language (EN)
#>  collate  en_US.UTF-8
#>  ctype    en_US.UTF-8
#>  tz       Asia/Shanghai
#>  date     2024-02-06
#>  pandoc   3.1.3 @ /Users/zouhua/opt/anaconda3/bin/ (via rmarkdown)
#> 
#> ─ Packages ───────────────────────────────────────────────
#>  package          * version     date (UTC) lib source
#>  AnnotationDbi      1.64.1      2023-11-03 [1] Bioconductor
#>  ape                5.7-1       2023-03-13 [1] CRAN (R 4.3.0)
#>  aPEAR            * 1.0.0       2023-06-12 [1] CRAN (R 4.3.0)
#>  aplot              0.2.2       2023-10-06 [1] CRAN (R 4.3.0)
#>  arules             1.7-7       2023-11-29 [1] CRAN (R 4.3.0)
#>  bayesbio           1.0.0       2016-05-24 [1] CRAN (R 4.3.0)
#>  Biobase            2.62.0      2023-10-24 [1] Bioconductor
#>  BiocGenerics       0.48.1      2023-11-01 [1] Bioconductor
#>  BiocParallel       1.36.0      2023-10-24 [1] Bioconductor
#>  Biostrings         2.70.2      2024-01-28 [1] Bioconductor 3.18 (R 4.3.2)
#>  bit                4.0.5       2022-11-15 [1] CRAN (R 4.3.0)
#>  bit64              4.0.5       2020-08-30 [1] CRAN (R 4.3.0)
#>  bitops             1.0-7       2021-04-24 [1] CRAN (R 4.3.0)
#>  blob               1.2.4       2023-03-17 [1] CRAN (R 4.3.0)
#>  bookdown           0.37        2023-12-01 [1] CRAN (R 4.3.0)
#>  bslib              0.6.1       2023-11-28 [1] CRAN (R 4.3.0)
#>  cachem             1.0.8       2023-05-01 [1] CRAN (R 4.3.0)
#>  cli                3.6.2       2023-12-11 [1] CRAN (R 4.3.0)
#>  clusterProfiler    4.10.0      2023-10-24 [1] Bioconductor
#>  codetools          0.2-19      2023-02-01 [1] CRAN (R 4.3.1)
#>  colorspace         2.1-0       2023-01-23 [1] CRAN (R 4.3.0)
#>  cowplot            1.1.3       2024-01-22 [1] CRAN (R 4.3.2)
#>  crayon             1.5.2       2022-09-29 [1] CRAN (R 4.3.0)
#>  data.table         1.15.0      2024-01-30 [1] CRAN (R 4.3.2)
#>  DBI                1.2.1       2024-01-12 [1] CRAN (R 4.3.0)
#>  devtools           2.4.5       2022-10-11 [1] CRAN (R 4.3.0)
#>  digest             0.6.34      2024-01-11 [1] CRAN (R 4.3.0)
#>  DOSE               3.28.2      2023-12-10 [1] Bioconductor
#>  downlit            0.4.3       2023-06-29 [1] CRAN (R 4.3.0)
#>  dplyr            * 1.1.4       2023-11-17 [1] CRAN (R 4.3.0)
#>  ellipsis           0.3.2       2021-04-29 [1] CRAN (R 4.3.0)
#>  enrichplot         1.22.0      2023-10-24 [1] Bioconductor
#>  evaluate           0.23        2023-11-01 [1] CRAN (R 4.3.0)
#>  expm               0.999-9     2024-01-11 [1] CRAN (R 4.3.0)
#>  fansi              1.0.6       2023-12-08 [1] CRAN (R 4.3.0)
#>  farver             2.1.1       2022-07-06 [1] CRAN (R 4.3.0)
#>  fastmap            1.1.1       2023-02-24 [1] CRAN (R 4.3.0)
#>  fastmatch          1.1-4       2023-08-18 [1] CRAN (R 4.3.0)
#>  FELLA            * 1.22.0      2023-10-24 [1] Bioconductor
#>  fgsea              1.28.0      2023-10-24 [1] Bioconductor
#>  forcats          * 1.0.0       2023-01-29 [1] CRAN (R 4.3.0)
#>  fs                 1.6.3       2023-07-20 [1] CRAN (R 4.3.0)
#>  generics           0.1.3       2022-07-05 [1] CRAN (R 4.3.0)
#>  GenomeInfoDb       1.38.5      2023-12-28 [1] Bioconductor 3.18 (R 4.3.2)
#>  GenomeInfoDbData   1.2.11      2024-01-24 [1] Bioconductor
#>  ggforce            0.4.1       2022-10-04 [1] CRAN (R 4.3.0)
#>  ggfun              0.1.4       2024-01-19 [1] CRAN (R 4.3.0)
#>  ggplot2          * 3.4.4       2023-10-12 [1] CRAN (R 4.3.0)
#>  ggplotify          0.1.2       2023-08-09 [1] CRAN (R 4.3.0)
#>  ggraph             2.1.0       2022-10-09 [1] CRAN (R 4.3.0)
#>  ggrepel            0.9.5       2024-01-10 [1] CRAN (R 4.3.0)
#>  ggtree             3.10.0      2023-10-24 [1] Bioconductor
#>  glue               1.7.0       2024-01-09 [1] CRAN (R 4.3.0)
#>  GO.db              3.18.0      2024-02-06 [1] Bioconductor
#>  GOSemSim           2.28.1      2024-01-17 [1] Bioconductor 3.18 (R 4.3.2)
#>  graph              1.80.0      2023-10-24 [1] Bioconductor
#>  graphlayouts       1.1.0       2024-01-19 [1] CRAN (R 4.3.0)
#>  gridExtra          2.3         2017-09-09 [1] CRAN (R 4.3.0)
#>  gridGraphics       0.5-1       2020-12-13 [1] CRAN (R 4.3.0)
#>  gson               0.1.0       2023-03-07 [1] CRAN (R 4.3.0)
#>  gtable             0.3.4       2023-08-21 [1] CRAN (R 4.3.0)
#>  HDO.db             0.99.1      2024-02-06 [1] Bioconductor
#>  highr              0.10        2022-12-22 [1] CRAN (R 4.3.0)
#>  hms                1.1.3       2023-03-21 [1] CRAN (R 4.3.0)
#>  htmltools          0.5.7       2023-11-03 [1] CRAN (R 4.3.0)
#>  htmlwidgets        1.6.4       2023-12-06 [1] CRAN (R 4.3.0)
#>  httpuv             1.6.14      2024-01-26 [1] CRAN (R 4.3.2)
#>  httr               1.4.7       2023-08-15 [1] CRAN (R 4.3.0)
#>  igraph           * 2.0.1.1     2024-01-30 [1] CRAN (R 4.3.2)
#>  IRanges            2.36.0      2023-10-24 [1] Bioconductor
#>  jquerylib          0.1.4       2021-04-26 [1] CRAN (R 4.3.0)
#>  jsonlite           1.8.8       2023-12-04 [1] CRAN (R 4.3.0)
#>  KEGGgraph          1.62.0      2023-10-24 [1] Bioconductor
#>  KEGGREST         * 1.42.0      2023-10-24 [1] Bioconductor
#>  knitr              1.45        2023-10-30 [1] CRAN (R 4.3.0)
#>  labeling           0.4.3       2023-08-29 [1] CRAN (R 4.3.0)
#>  later              1.3.2       2023-12-06 [1] CRAN (R 4.3.0)
#>  lattice            0.21-8      2023-04-05 [1] CRAN (R 4.3.1)
#>  lazyeval           0.2.2       2019-03-15 [1] CRAN (R 4.3.0)
#>  lifecycle          1.0.4       2023-11-07 [1] CRAN (R 4.3.0)
#>  lsa                0.73.3      2022-05-09 [1] CRAN (R 4.3.0)
#>  lubridate        * 1.9.3       2023-09-27 [1] CRAN (R 4.3.0)
#>  magrittr           2.0.3       2022-03-30 [1] CRAN (R 4.3.0)
#>  MASS               7.3-60      2023-05-04 [1] CRAN (R 4.3.1)
#>  Matrix             1.6-5       2024-01-11 [1] CRAN (R 4.3.0)
#>  MCL                1.0         2015-03-11 [1] CRAN (R 4.3.0)
#>  memoise            2.0.1       2021-11-26 [1] CRAN (R 4.3.0)
#>  mime               0.12        2021-09-28 [1] CRAN (R 4.3.0)
#>  miniUI             0.1.1.1     2018-05-18 [1] CRAN (R 4.3.0)
#>  munsell            0.5.0       2018-06-12 [1] CRAN (R 4.3.0)
#>  nlme               3.1-162     2023-01-31 [1] CRAN (R 4.3.1)
#>  org.Hs.eg.db       3.18.0      2024-02-06 [1] Bioconductor
#>  patchwork          1.2.0       2024-01-08 [1] CRAN (R 4.3.0)
#>  pathview         * 1.42.0      2023-10-24 [1] Bioconductor
#>  pillar             1.9.0       2023-03-22 [1] CRAN (R 4.3.0)
#>  pkgbuild           1.4.3       2023-12-10 [1] CRAN (R 4.3.0)
#>  pkgconfig          2.0.3       2019-09-22 [1] CRAN (R 4.3.0)
#>  pkgload            1.3.4       2024-01-16 [1] CRAN (R 4.3.0)
#>  plyr               1.8.9       2023-10-02 [1] CRAN (R 4.3.0)
#>  png                0.1-8       2022-11-29 [1] CRAN (R 4.3.0)
#>  polyclip           1.10-6      2023-09-27 [1] CRAN (R 4.3.0)
#>  profvis            0.3.8       2023-05-02 [1] CRAN (R 4.3.0)
#>  promises           1.2.1       2023-08-10 [1] CRAN (R 4.3.0)
#>  purrr            * 1.0.2       2023-08-10 [1] CRAN (R 4.3.0)
#>  qvalue             2.34.0      2023-10-24 [1] Bioconductor
#>  R6                 2.5.1       2021-08-19 [1] CRAN (R 4.3.0)
#>  RColorBrewer       1.1-3       2022-04-03 [1] CRAN (R 4.3.0)
#>  Rcpp               1.0.12      2024-01-09 [1] CRAN (R 4.3.0)
#>  RCurl              1.98-1.14   2024-01-09 [1] CRAN (R 4.3.0)
#>  readr            * 2.1.5       2024-01-10 [1] CRAN (R 4.3.0)
#>  remotes            2.4.2.1     2023-07-18 [1] CRAN (R 4.3.0)
#>  reshape2           1.4.4       2020-04-09 [1] CRAN (R 4.3.0)
#>  Rgraphviz          2.46.0      2023-10-24 [1] Bioconductor
#>  rlang              1.1.3       2024-01-10 [1] CRAN (R 4.3.0)
#>  rmarkdown          2.25        2023-09-18 [1] CRAN (R 4.3.0)
#>  RSQLite            2.3.5       2024-01-21 [1] CRAN (R 4.3.0)
#>  rstudioapi         0.15.0      2023-07-07 [1] CRAN (R 4.3.0)
#>  S4Vectors          0.40.2      2023-11-23 [1] Bioconductor
#>  sass               0.4.8       2023-12-06 [1] CRAN (R 4.3.0)
#>  scales             1.3.0       2023-11-28 [1] CRAN (R 4.3.0)
#>  scatterpie         0.2.1       2023-06-07 [1] CRAN (R 4.3.0)
#>  sessioninfo        1.2.2       2021-12-06 [1] CRAN (R 4.3.0)
#>  shadowtext         0.1.3       2024-01-19 [1] CRAN (R 4.3.0)
#>  shiny              1.8.0       2023-11-17 [1] CRAN (R 4.3.0)
#>  SnowballC          0.7.1       2023-04-25 [1] CRAN (R 4.3.0)
#>  stringi            1.8.3       2023-12-11 [1] CRAN (R 4.3.0)
#>  stringr          * 1.5.1       2023-11-14 [1] CRAN (R 4.3.0)
#>  tibble           * 3.2.1       2023-03-20 [1] CRAN (R 4.3.0)
#>  tidygraph          1.3.0       2023-12-18 [1] CRAN (R 4.3.0)
#>  tidyr            * 1.3.1       2024-01-24 [1] CRAN (R 4.3.2)
#>  tidyselect         1.2.0       2022-10-10 [1] CRAN (R 4.3.0)
#>  tidytree           0.4.6       2023-12-12 [1] CRAN (R 4.3.0)
#>  tidyverse        * 2.0.0       2023-02-22 [1] CRAN (R 4.3.0)
#>  timechange         0.3.0       2024-01-18 [1] CRAN (R 4.3.0)
#>  treeio             1.26.0      2023-10-24 [1] Bioconductor
#>  tweenr             2.0.2       2022-09-06 [1] CRAN (R 4.3.0)
#>  tzdb               0.4.0       2023-05-12 [1] CRAN (R 4.3.0)
#>  urlchecker         1.0.1       2021-11-30 [1] CRAN (R 4.3.0)
#>  usethis            2.2.2       2023-07-06 [1] CRAN (R 4.3.0)
#>  utf8               1.2.4       2023-10-22 [1] CRAN (R 4.3.0)
#>  vctrs              0.6.5       2023-12-01 [1] CRAN (R 4.3.0)
#>  viridis            0.6.5       2024-01-29 [1] CRAN (R 4.3.2)
#>  viridisLite        0.4.2       2023-05-02 [1] CRAN (R 4.3.0)
#>  withr              3.0.0       2024-01-16 [1] CRAN (R 4.3.0)
#>  xfun               0.41        2023-11-01 [1] CRAN (R 4.3.0)
#>  XML                3.99-0.16.1 2024-01-22 [1] CRAN (R 4.3.2)
#>  xml2               1.3.6       2023-12-04 [1] CRAN (R 4.3.0)
#>  xtable             1.8-4       2019-04-21 [1] CRAN (R 4.3.0)
#>  XVector            0.42.0      2023-10-24 [1] Bioconductor
#>  yaml               2.3.8       2023-12-11 [1] CRAN (R 4.3.0)
#>  yulab.utils        0.1.4       2024-01-28 [1] CRAN (R 4.3.2)
#>  zlibbioc           1.48.0      2023-10-24 [1] Bioconductor
#> 
#>  [1] /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library
#> 
#> ──────────────────────────────────────────────────────────
```


## Reference

+ [FELLA github](https://github.com/b2slab/FELLA)

+ [The MetaRbolomics book](https://rformassspectrometry.github.io/metaRbolomics-book/)
