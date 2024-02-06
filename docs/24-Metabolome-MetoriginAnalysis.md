


# MetOrigin Analysis {#MetOriginAnalysis}


微生物群及其代谢产物与人类健康和疾病密切相关。然而，理解微生物组和代谢物之间复杂的相互作用是具有挑战性的。

在研究肠道代谢物时，代谢物的来源是一个无法避免的问题即代谢物到底是来自肠道微生物的代谢还是宿主本身代谢产生的。2022年Yu, Gang et.al.发表的[@yu2022metorigin]提供了一个可以区分为微生物还是宿主代谢物的工具。

该工具的核心是一个人工检验过具有微生物和宿主等标签的代谢物数据库。

+ 微生物独有的代谢物

+ 宿主独有代谢物

+ 两者都可以代谢的代谢物

+ 其他类型（未知）


该工具没有提供R版本，需要到其网站使用。本教程只提供如何准备网站所需的输入文件。


## 加载R包

```r
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
library(tidyverse)


# rm(list = ls())
options(stringsAsFactors = F)
options(future.globals.maxSize = 1000 * 1024^2)

grp_names <- c("None", "Mild", "Moderate", "Severe")
grp_colors <- c("#7DD06F", "#844081", "#688EC1", "#C17E73")
```


## 导入数据

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



## 准备输入文件

网站需要的输入文件必须要以代谢物的名称如**HMDBID**等作为代谢物唯一标识，本文也将用**HMDBID**。除此之外，还需要设置表明代谢物差异富集方向的列**Diff**。

> The metabolite table must contain at least one column of “HMDBID”, “KEGGID” or “Name”, and a column of 0/1 values indicating statistical significance (1-significant, 0-nonsignificant). If the “Diff” column is missing, all metabolites will be considered as differential metabolites.


```r
get_metabolites <- function(
  dat,
  group_names,
  index_names = c("FoldChange", "Log2FoldChange", "VIP", "CorPvalue", "Pvalue", "AdjustedPvalue"),
  index_cutoff = c(1, 1, 1, 0.05, 0.05, 0.2)) {
  
  
  colnames(dat)[which(colnames(dat) == "SampleIDHMDBID")] <- "HMDB"
  colnames(dat)[which(colnames(dat) == "KEGG")] <- "cpd_ID"
  colnames(dat)[which(colnames(dat) == "BIOCHEMICAL")] <- "Compounds"
  dat$HMDB <- gsub(",\\S+", "", dat$HMDB)
  
  temp_dat <- dat %>%
    dplyr::filter(Block2 %in% group_names) %>%
    dplyr::filter(HMDB != "-")
  
  colnames(temp_dat)[which(colnames(temp_dat) == index_names[1])] <- "DA_index1"
  colnames(temp_dat)[which(colnames(temp_dat) == index_names[2])] <- "DA_index2"
  
  temp_dat_diff <- temp_dat %>%
    dplyr::filter(abs(DA_index1) > index_cutoff[1]) %>%
    dplyr::filter(DA_index2 < index_cutoff[2]) %>%
    dplyr::mutate(Diff = 1)
  
  if (nrow(temp_dat_diff) == 0) {
    stop("Beyond these thresholds, no significant metabolites were selected")
  }
  
  temp_dat_nodiff <- temp_dat %>%
    dplyr::filter(!HMDB %in% temp_dat_diff$HMDB) %>%
    dplyr::mutate(Diff = 0)
  
  res <- rbind(temp_dat_diff, temp_dat_nodiff) %>%
    dplyr::select(HMDB, cpd_ID, Compounds, DA_index1, DA_index2, Diff) %>%
    dplyr::rename(HMDBID = HMDB,
                  KEGGID = cpd_ID,
                  Name = Compounds) %>%
    dplyr::select(HMDBID, KEGGID, Name, Diff)

  return(res)
}

pre_data <- get_metabolites(
  dat = datSignif,
  group_names = "None vs Severe",
  index_names = c("Log2FoldChange", "AdjustedPvalue"),
  index_cutoff = c(0, 0.9))
  
head(pre_data)
#>         HMDBID KEGGID
#>         <char> <char>
#> 1: HMDB0000656 R00900
#> 2: HMDB0000187 C00065
#> 3: HMDB0009783   <NA>
#> 4: HMDB0007102 C13861
#> 5: HMDB0004950   <NA>
#> 6: HMDB0000177 C00135
#>                                          Name  Diff
#>                                        <char> <num>
#> 1:             cysteine-glutathione disulfide     1
#> 2:                                     serine     1
#> 3:      1-palmitoyl-2-oleoyl-GPI (16:0/18:1)*     1
#> 4: palmitoyl-oleoyl-glycerol (16:0/18:1) [2]*     1
#> 5:       N-stearoyl-sphingosine (d18:1/18:0)*     1
#> 6:                                  histidine     1
```


## 输出结果文件

+ the csv file could be used as inputs in the MetOrigin website.


```r
if(!dir.exists("./InputData/result/MetOrigin")) {
  dir.create("./InputData/result/MetOrigin", recursive = TRUE)
}

write.csv(pre_data, "./InputData/result/MetOrigin/MetOrigin_inputs.csv", row.names = F)
```


## MetOrigin User Tutorial

MetOrigin comprises five parts:

+ **Load data**

> Please select the analysis mode at first "Simple MetOrigin Analysis" or "Deep MetOrigin Analysis". Then you can try "Load Example Data" or upload your own data. Host information needs to be confirmed before moving to the next step.

+ **Origin Analysis**

> This step is to identify where metabolites come from: host, bacteria, both, or unknown?

+ **Function Analysis**

> This step is to perform the metabolic pathway enrichment analysis according to different categories of metabolites: metabolites belonging to the host, bacteria, or both.

+ **Sankey Network**

> This step is to link all the possible bacteria that may participate in a specific metabolic reaction, helping you to identify important interplay between bacteria and metabolites.

+ **Download Results**

> All the analysis results can be downloaded for your further data exploration. 

**See more details please go to the below tutorial ([MetOrigin website](https://metorigin.met-bioinformatics.cn/))**:

![MetOrigin User Tutorial](./InputData/Tutorial/MetOrigin_User_Tutorial.pdf){width=100% height=800}



## 结果解析

+ 溯源分析：条形图和韦恩图展示代谢物的来源

+ 功能分析：分别用不同来源的差异代谢物做超几何富集分析，判断它们富集在哪些通路

+ 桑基图网络分析：某个特定代谢通路的微生物对代谢物的贡献程度，不同颜色表明微生物对代谢物的上下调关系


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
#>  package     * version date (UTC) lib source
#>  bookdown      0.37    2023-12-01 [1] CRAN (R 4.3.0)
#>  bslib         0.6.1   2023-11-28 [1] CRAN (R 4.3.0)
#>  cachem        1.0.8   2023-05-01 [1] CRAN (R 4.3.0)
#>  cli           3.6.2   2023-12-11 [1] CRAN (R 4.3.0)
#>  colorspace    2.1-0   2023-01-23 [1] CRAN (R 4.3.0)
#>  data.table    1.15.0  2024-01-30 [1] CRAN (R 4.3.2)
#>  devtools      2.4.5   2022-10-11 [1] CRAN (R 4.3.0)
#>  digest        0.6.34  2024-01-11 [1] CRAN (R 4.3.0)
#>  downlit       0.4.3   2023-06-29 [1] CRAN (R 4.3.0)
#>  dplyr       * 1.1.4   2023-11-17 [1] CRAN (R 4.3.0)
#>  ellipsis      0.3.2   2021-04-29 [1] CRAN (R 4.3.0)
#>  evaluate      0.23    2023-11-01 [1] CRAN (R 4.3.0)
#>  fansi         1.0.6   2023-12-08 [1] CRAN (R 4.3.0)
#>  fastmap       1.1.1   2023-02-24 [1] CRAN (R 4.3.0)
#>  forcats     * 1.0.0   2023-01-29 [1] CRAN (R 4.3.0)
#>  fs            1.6.3   2023-07-20 [1] CRAN (R 4.3.0)
#>  generics      0.1.3   2022-07-05 [1] CRAN (R 4.3.0)
#>  ggplot2     * 3.4.4   2023-10-12 [1] CRAN (R 4.3.0)
#>  glue          1.7.0   2024-01-09 [1] CRAN (R 4.3.0)
#>  gtable        0.3.4   2023-08-21 [1] CRAN (R 4.3.0)
#>  hms           1.1.3   2023-03-21 [1] CRAN (R 4.3.0)
#>  htmltools     0.5.7   2023-11-03 [1] CRAN (R 4.3.0)
#>  htmlwidgets   1.6.4   2023-12-06 [1] CRAN (R 4.3.0)
#>  httpuv        1.6.14  2024-01-26 [1] CRAN (R 4.3.2)
#>  jquerylib     0.1.4   2021-04-26 [1] CRAN (R 4.3.0)
#>  jsonlite      1.8.8   2023-12-04 [1] CRAN (R 4.3.0)
#>  knitr         1.45    2023-10-30 [1] CRAN (R 4.3.0)
#>  later         1.3.2   2023-12-06 [1] CRAN (R 4.3.0)
#>  lifecycle     1.0.4   2023-11-07 [1] CRAN (R 4.3.0)
#>  lubridate   * 1.9.3   2023-09-27 [1] CRAN (R 4.3.0)
#>  magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.3.0)
#>  memoise       2.0.1   2021-11-26 [1] CRAN (R 4.3.0)
#>  mime          0.12    2021-09-28 [1] CRAN (R 4.3.0)
#>  miniUI        0.1.1.1 2018-05-18 [1] CRAN (R 4.3.0)
#>  munsell       0.5.0   2018-06-12 [1] CRAN (R 4.3.0)
#>  pillar        1.9.0   2023-03-22 [1] CRAN (R 4.3.0)
#>  pkgbuild      1.4.3   2023-12-10 [1] CRAN (R 4.3.0)
#>  pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.3.0)
#>  pkgload       1.3.4   2024-01-16 [1] CRAN (R 4.3.0)
#>  profvis       0.3.8   2023-05-02 [1] CRAN (R 4.3.0)
#>  promises      1.2.1   2023-08-10 [1] CRAN (R 4.3.0)
#>  purrr       * 1.0.2   2023-08-10 [1] CRAN (R 4.3.0)
#>  R6            2.5.1   2021-08-19 [1] CRAN (R 4.3.0)
#>  Rcpp          1.0.12  2024-01-09 [1] CRAN (R 4.3.0)
#>  readr       * 2.1.5   2024-01-10 [1] CRAN (R 4.3.0)
#>  remotes       2.4.2.1 2023-07-18 [1] CRAN (R 4.3.0)
#>  rlang         1.1.3   2024-01-10 [1] CRAN (R 4.3.0)
#>  rmarkdown     2.25    2023-09-18 [1] CRAN (R 4.3.0)
#>  rstudioapi    0.15.0  2023-07-07 [1] CRAN (R 4.3.0)
#>  sass          0.4.8   2023-12-06 [1] CRAN (R 4.3.0)
#>  scales        1.3.0   2023-11-28 [1] CRAN (R 4.3.0)
#>  sessioninfo   1.2.2   2021-12-06 [1] CRAN (R 4.3.0)
#>  shiny         1.8.0   2023-11-17 [1] CRAN (R 4.3.0)
#>  stringi       1.8.3   2023-12-11 [1] CRAN (R 4.3.0)
#>  stringr     * 1.5.1   2023-11-14 [1] CRAN (R 4.3.0)
#>  tibble      * 3.2.1   2023-03-20 [1] CRAN (R 4.3.0)
#>  tidyr       * 1.3.1   2024-01-24 [1] CRAN (R 4.3.2)
#>  tidyselect    1.2.0   2022-10-10 [1] CRAN (R 4.3.0)
#>  tidyverse   * 2.0.0   2023-02-22 [1] CRAN (R 4.3.0)
#>  timechange    0.3.0   2024-01-18 [1] CRAN (R 4.3.0)
#>  tzdb          0.4.0   2023-05-12 [1] CRAN (R 4.3.0)
#>  urlchecker    1.0.1   2021-11-30 [1] CRAN (R 4.3.0)
#>  usethis       2.2.2   2023-07-06 [1] CRAN (R 4.3.0)
#>  utf8          1.2.4   2023-10-22 [1] CRAN (R 4.3.0)
#>  vctrs         0.6.5   2023-12-01 [1] CRAN (R 4.3.0)
#>  withr         3.0.0   2024-01-16 [1] CRAN (R 4.3.0)
#>  xfun          0.41    2023-11-01 [1] CRAN (R 4.3.0)
#>  xml2          1.3.6   2023-12-04 [1] CRAN (R 4.3.0)
#>  xtable        1.8-4   2019-04-21 [1] CRAN (R 4.3.0)
#>  yaml          2.3.8   2023-12-11 [1] CRAN (R 4.3.0)
#> 
#>  [1] /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library
#> 
#> ──────────────────────────────────────────────────────────
```


## Reference

+ [metorigin website](https://metorigin.met-bioinformatics.cn/home/) 

