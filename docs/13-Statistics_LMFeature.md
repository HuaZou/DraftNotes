


# Linear Model on Microbial Community {#LinearModelonMicrobialCommunity}


在微生物研究中，通常会使用基于距离矩阵的置换检验（如PERMANOVA）判断环境因素和微生物群落结构的相关性关系。本文使用另一套分析方法计算环境因素对微生物群落差异和微生物相对丰度差异的贡献，它使用到的是：

+ 相关性分析：计算环境因素和各个物种相对丰度的相关性大小；

+ 多元线性回归分析：计算主要的环境因素对各个物种相对丰度的贡献（$R^2$）；

+ 拆分线性回归的总体贡献：计算每一个环境因素对各个物种相对丰度的重要程度（类似方差解分析）。


本文为了确定影响微生物（肠道和口腔）群落差异和微生物相对丰度差异的人体特征变量。


**声明**：本文参考了小白鱼的《仿一篇文献的相关性分析和线性模型评估影响群落组成的重要环境变量》。


## 加载R包


```r
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
library(tidyverse)
library(phyloseq)

library(caret)
library(leaps)
library(relaimpo)

# rm(list = ls())
options(stringsAsFactors = F)
options(future.globals.maxSize = 1000 * 1024^2)
```


## 导入数据

对数据[OmicsDataSet-Zeybel et al. - 2022.xlsx](https://github.com/HuaZou/DraftNotes/blob/main/InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx)处理后生成的，可参考数据生成和预处理章节。


+ ps_gut：来自宿主肠道微生物组 (见[Zeybel Dataset](https://zouhua.top/DraftNotes/ZeybelDataset.html)章节)

+ ps_oral：来自宿主肠道微生物组 (见[Zeybel Dataset](https://zouhua.top/DraftNotes/ZeybelDataset.html)章节)

> ```R
> saveRDS(ps_gut, "./InputData/result/Zeybel_2022_gut_MGS_ps.RDS", compress = TRUE)
> saveRDS(ps_ora, "./InputData/result/Zeybel_2022_oral_MGS_ps.RDS", compress = TRUE)
> ```


```r
ps_gut <- readRDS("./InputData/result/Zeybel_2022_gut_MGS_ps.RDS")
ps_oral <- readRDS("./InputData/result/Zeybel_2022_oral_MGS_ps.RDS")
```


## 数据准备

+ 人体特征指标：如Liver_fat, Creatinine, Body_mass_index等

+ 门水平微生物表达谱：phylum levels


```r
meta_gut <- ps_gut@sam_data %>%
  data.frame() %>%
  na.omit() %>%
  dplyr::select(6:14)

meta_oral <- ps_oral@sam_data %>%
  data.frame() %>%
  na.omit() %>%
  dplyr::select(6:14)

sid <- intersect(rownames(meta_gut), rownames(meta_oral))

# remotes::install_github("HuaZou/MicrobiomeAnalysis")
prof_gut <- MicrobiomeAnalysis::aggregate_taxa(
  x = ps_gut,
  level = "Phylum")
prof_gut <- MicrobiomeAnalysis::trim_prevalence(
  object = prof_gut,
  cutoff = 0.1,
  trim = "feature")
prof_gut <- prof_gut@otu_table %>%
  data.frame()
prof_gut_final <- prof_gut[, pmatch(sid, colnames(prof_gut))]

prof_oral <- MicrobiomeAnalysis::aggregate_taxa(
  x = ps_oral,
  level = "Phylum")
prof_oral <- MicrobiomeAnalysis::trim_prevalence(
  object = prof_oral,
  cutoff = 0.1,
  trim = "feature")
prof_oral <- prof_oral@otu_table %>%
  data.frame()
prof_oral_final <- prof_oral[, pmatch(sid, colnames(prof_oral))]

meta_final <- meta_oral[pmatch(sid, rownames(meta_oral)), ,]

rm(ps_gut, ps_oral, meta_gut, meta_oral, sid, prof_gut, prof_oral)
```


## 相关性分析

人体特征和物种相关性分析，分别计算每个人体特征和各个物种丰度的相关系数（采用spearman相关系数）。


```r
run_cor <- function(
    data_sam,
    data_otu,
    columns  = NULL,
    method = c("spearman", "pearson", "kendall"),
    p_adjust = c("none", "fdr", "bonferroni", "holm",
                 "hochberg", "hommel", "BH", "BY")) {

  # data_sam = meta_final
  # data_otu = prof_gut_final
  # columns = NULL
  # method = "spearman"
  # p_adjust = "BH"

  if (is.null(method)) {
    method <- "spearman"
  } else {
    method <- match.arg(
      method,
      c("spearman", "pearson", "kendall")
    )
  }
  # p_adjust
  p_adjust <- match.arg(p_adjust,
                        c("none", "fdr", "bonferroni", "holm",
                          "hochberg", "hommel", "BH", "BY")
  )

  interset_sampleid <- dplyr::intersect(
      colnames(data_otu), rownames(data_sam))

  data_otu_interset <- data_otu %>% 
    as.data.frame() %>%
    dplyr::select(dplyr::all_of(interset_sampleid))
  data_sam_interset <- data_sam %>% 
    as.data.frame()

  data_sam_interset_final <- data_sam_interset[pmatch(interset_sampleid, 
                                                      rownames(data_sam_interset)), , F]

  if (!all(colnames(data_otu_interset) == rownames(data_sam_interset_final))) {
    stop("The order of SampleID is wrong, please check your inputdata")
  }

  sam_tab <- data_sam_interset_final %>%
      as.data.frame() %>%
      tibble::rownames_to_column("TempRowNames")
  otu_tab_t <- data_otu_interset %>%
      as.data.frame() %>%
      base::t() %>%
      as.data.frame()

  # columns for test
  if (!is.null(columns)) {
    sam_tab <- sam_tab %>%
      tibble::column_to_rownames("TempRowNames") %>%
      dplyr::select(dplyr::all_of(columns))
  } else {
    sam_tab <- sam_tab %>%
      tibble::column_to_rownames("TempRowNames")    
  }

  if (!all(rownames(otu_tab_t) == rownames(sam_tab))) {
    stop("The order of SampleID is wrong, please check your input")
  }

  # calculate the association between individual taxa and factors
  res <- data.frame()
  for (i in 1:ncol(sam_tab)) {
    # whether all the elements are numeric
    if (!is.numeric(sam_tab[, i])) {
      stop("Please check your input, values in ", colnames(sam_tab)[i]," are not numeric")
    }
    for (j in 1:ncol(otu_tab_t)) {
      mdat <- data.frame(x=sam_tab[, i],
                         y=otu_tab_t[, j]) %>%
        na.omit()
      fit <- stats::cor.test(mdat$x, mdat$y, method = method)
      temp <- data.frame(Phenotype=colnames(sam_tab)[i],
                         FeatureID=colnames(otu_tab_t)[j],
                         Statistic=fit$statistic,
                         Rho=fit$estimate,
                         Pvalue=fit$p.value)
      res <- rbind(res, temp)
    }
  }

  res$AdjustedPvalue <- p.adjust(as.numeric(res$Pvalue), method = p_adjust)
  rownames(res) <- NULL

  return(res)
}

gut_corres <- run_cor(
  data_sam = meta_final,
  data_otu = prof_gut_final,
  p_adjust = "BH")

oral_corres <- run_cor(
  data_sam = meta_final,
  data_otu = prof_oral_final,
  p_adjust = "BH")

head(gut_corres[, 1:3], 2)
#>   Phenotype         FeatureID Statistic
#> 1 Liver_fat p__Actinobacteria      7914
#> 2 Liver_fat  p__Bacteroidetes      4858
```



## 线性回归分析

人体特征指标对各个物种的贡献度可通过多元线性回归分析，即物种作为Y响应变量，人体特征作为X变量，通过回归的$R^2$解析总解释度。

人体指标可能存在共线性的情况，可以在线性回归计算时候选择前向或后向回归方式筛选重要的变量。两种方法对结果影响较大，可能会得到大相径庭的结果。


### 选择变量分析思路


+ 通过前向选择变量的分析思路

> 1. 通过前向逐步回归，在所有人体特征变量筛选重要的变量，尽可能减少人体特征共线性以获取简约模型，又同时尽可能保证模型的总解释率不要损失很多；
>
> 2. 使用选择的人体特征变量，分别拟合与每个物种丰度的多元线性回归，以期通过人体特征变量来解释物种丰度组成的差异；
>
> 3. 对于在第（2）步中获得的最优模型，提取它们的总解释率（即线性模型的R2或校正后的R2）等信息；
>
> 4. 对于在第（2）步中获得的最优模型，通过类似方差分解分析的方法，评估主要的人体特征对于物种丰度总方差的贡献，实现定量分析主要的人体特征相对重要性的目的。


+ 通过后向选择变量的分析思路

> 1. 使用所有人体特征变量，分别拟合与每个物种丰度的多元线性回归，以期通过人体特征变量来解释物种丰度组成的差异；
>
> 2. 通过后向逐步回归，在构建好的每个线性模型中筛选重要的变量，尽可能减少人体特征变量共线性以获取简约模型，又同时尽可能保证模型的总解释率不要损失很多；
>
> 3. 对于在第（2）步中获得的最优模型，提取它们的总解释率（即线性模型的R2或校正后的R2）等信息；
>
> 4. 对于在第（2）步中获得的最优模型，通过类似方差分解分析的方法，评估主要的人体特征对于物种丰度总方差的贡献，实现定量分析主要的人体特征相对重要性的目的。


### 函数

method参数控制筛选变量的方向：

+ "leapForward": 前向逐步回归 (forward selection)

+ "leapBackward": 后向逐步回归 (backward selection)

+ "leapSeq": 逐步回归 (stepwise selection)


```r
run_lm <- function(
    data_sam,
    data_otu,
    columns  = NULL,
    method = c("leapForward", "leapBackward", "leapSeq")) {

  # data_sam = meta_final
  # data_otu = prof_gut_final
  # columns = NULL
  # method = "leapForward"

  if (is.null(method)) {
    method <- "leapForward"
  } else {
    method <- match.arg(
      method,
      c("leapForward", "leapBackward", "leapSeq")
    )
  }

  interset_sampleid <- dplyr::intersect(
      colnames(data_otu), rownames(data_sam))

  data_otu_interset <- data_otu %>% 
    as.data.frame() %>%
    dplyr::select(dplyr::all_of(interset_sampleid))
  data_sam_interset <- data_sam %>% 
    as.data.frame()

  data_sam_interset_final <- data_sam_interset[pmatch(interset_sampleid, 
                                                      rownames(data_sam_interset)), , F]

  if (!all(colnames(data_otu_interset) == rownames(data_sam_interset_final))) {
    stop("The order of SampleID is wrong, please check your inputdata")
  }

  sam_tab <- data_sam_interset_final %>%
      as.data.frame() %>%
      tibble::rownames_to_column("TempRowNames")
  otu_tab_t <- data_otu_interset %>%
      as.data.frame() %>%
      base::t() %>%
      as.data.frame()

  # columns for test
  if (!is.null(columns)) {
    sam_tab <- sam_tab %>%
      tibble::column_to_rownames("TempRowNames") %>%
      dplyr::select(dplyr::all_of(columns))
  } else {
    sam_tab <- sam_tab %>%
      tibble::column_to_rownames("TempRowNames")    
  }

  if (!all(rownames(otu_tab_t) == rownames(sam_tab))) {
    stop("The order of SampleID is wrong, please check your input")
  }
  
  LM_func <- function(dat, feature_name) {
    
    feature <- feature_name
    colnames(dat)[1] <- "Taxa"
    
    set.seed(123)
    
    # 变量选择
    train.control <- trainControl(method = "cv", number = 3)
    step.model <- train(Taxa ~., data = dat,
                        method = method, 
                        tuneGrid = data.frame(nvmax = 1:5),
                        trControl = train.control)
    
    model.coef <- coef(step.model$finalModel, as.numeric(step.model$bestTune)) %>%
      as.data.frame()
    final_index <- rownames(model.coef)[-1]
    
    # 基于上述选择的变量，使用 lm()拟合变量与各物种丰度的多元线性回归，获取最优模型
    if (length(final_index) < 1) {
      result_model <- data.frame(FeatureID = feature,
                                 R2 = NA,
                                 AdjustedR2 = NA,
                                 Pvalue = NA)
      result_import <- data.frame(TempID = colnames(dat)[-1],
                                 Value = NA) %>%
        tibble::column_to_rownames("TempID") %>%
        setNames(feature)
    } else {
      dat_new <- dat %>%
        dplyr::select(dplyr::all_of(c("Taxa", final_index)))
      lm_fit <- lm(Taxa ~., data = dat_new)
      lm_stat <- summary(lm_fit)
      R2 <- lm_stat$r.squared 
      adjR2 <- lm_stat$adj.r.squared 
      fvalue <- lm_stat$fstatistic 
      Pvalue <- pf(fvalue[1], fvalue[2], fvalue[3])
      
      result_model <- data.frame(FeatureID = feature,
                                 R2 = R2,
                                 AdjustedR2 = adjR2,
                                 Pvalue = Pvalue)
      
      # 多元线性回归（已通过变量选择后的最优模型）中各环境变量的相对重要性
      error <- tryCatch(
        expr = {
          crf <- relaimpo::calc.relimp(lm_fit, rela = FALSE)
        },
        error = function(e){
          message('calc.relimp Caught an error!')
          print(e)
        }
      )
      
      if (length(error) > 1) {
        result_import <- data.frame(TempID = colnames(dat)[-1],
                                     Value = NA) %>%
            tibble::column_to_rownames("TempID") %>%
            setNames(feature)  
      } else {
        temp_import <- crf@lmg %>%
          as.data.frame() %>%
          setNames(feature)
        temp_import_NA <- data.frame(TempID = c(setdiff(colnames(dat)[-1], final_index)),
                                     Value = NA) %>%
          tibble::column_to_rownames("TempID") %>%
          setNames(feature)
        
        result_import <- rbind(temp_import, temp_import_NA)        
      }
    }
    
    res <- list(model = result_model,
                import = result_import)
    
    return(res)
  }  
  

  res_model <- data.frame()
  res_import <- data.frame(matrix(NA, nrow = ncol(sam_tab), ncol = 0))
  
  for (i in 1:ncol(otu_tab_t)) {
    
    mdat <- cbind(otu_tab_t[, i, ], sam_tab)
    featureID <- colnames(otu_tab_t)[i]
    
    temp_list <- LM_func(dat = mdat, feature_name = featureID)
    
    res_model <- rbind(res_model, temp_list$model)
    res_import <- cbind(res_import, temp_list$import)
  }
  
  res <- list(model = res_model,
              import = res_import)
  
  return(res)
}

gut_LM <- run_lm(
  data_sam = meta_final,
  data_otu = prof_gut_final,
  method = "leapBackward")
#> <simpleError in solve.default(covg[diese, diese], matrix(covg[diese, andere],     length(diese), p + 1 - length(diese))): 'a' is 0-diml>
#> <simpleError in solve.default(covg[diese, diese], matrix(covg[diese, andere],     length(diese), p + 1 - length(diese))): 'a' is 0-diml>
#> <simpleError in solve.default(covg[diese, diese], matrix(covg[diese, andere],     length(diese), p + 1 - length(diese))): 'a' is 0-diml>
#> <simpleError in solve.default(covg[diese, diese], matrix(covg[diese, andere],     length(diese), p + 1 - length(diese))): 'a' is 0-diml>
#> <simpleError in solve.default(covg[diese, diese], matrix(covg[diese, andere],     length(diese), p + 1 - length(diese))): 'a' is 0-diml>

oral_LM <- run_lm(
  data_sam = meta_final,
  data_otu = prof_oral_final,
  method = "leapBackward")
#> <simpleError in solve.default(covg[diese, diese], matrix(covg[diese, andere],     length(diese), p + 1 - length(diese))): 'a' is 0-diml>
#> <simpleError in solve.default(covg[diese, diese], matrix(covg[diese, andere],     length(diese), p + 1 - length(diese))): 'a' is 0-diml>
#> <simpleError in solve.default(covg[diese, diese], matrix(covg[diese, andere],     length(diese), p + 1 - length(diese))): 'a' is 0-diml>
#> <simpleError in solve.default(covg[diese, diese], matrix(covg[diese, andere],     length(diese), p + 1 - length(diese))): 'a' is 0-diml>

head(oral_LM$model[, 1:3], 2)
#>                FeatureID        R2 AdjustedR2
#> value  p__Actinobacteria 0.2355881  0.2109297
#> value1  p__Bacteroidetes 0.3249810  0.2285498
```


## 可视化

+ 热图的颜色表示相关系数大小；

+ 柱状图表示人体特征变量对物种丰度差异的变异度解释；

+ 热图的圆圈大小表示人体特征变量对物种的贡献程度大小。


```r
get_LM_plot <- function(datCor, datLM)  {
  
  dat_cor <- datCor
  dat_cor[which(dat_cor$Pvalue < 0.001), "siglabel"] <- "***"
  dat_cor[which(dat_cor$Pvalue < 0.01 & dat_cor$Pvalue > 0.001), "siglabel"] <- "**"
  dat_cor[which(dat_cor$Pvalue < 0.05 & dat_cor$Pvalue > 0.01), "siglabel"] <- "*"
  
  dat_LM <- datLM
  dat_LM_R2 <- dat_LM$model
  dat_LM_import <- dat_LM$import
  
  dat_LM_R2[which(dat_LM_R2$Pvalue < 0.001), "siglabel"] <- "***"
  dat_LM_R2[which(dat_LM_R2$Pvalue < 0.01 & dat_LM_R2$Pvalue > 0.001), "siglabel"] <- "**"
  dat_LM_R2[which(dat_LM_R2$Pvalue < 0.05 & dat_LM_R2$Pvalue > 0.01), "siglabel"] <- "*"
  
  dat_LM_import_final <- dat_LM_import %>%
    tibble::rownames_to_column("Phenotype") %>%
    tidyr::gather(key = "FeatureID", value = "importance", -Phenotype)
  
  p1 <- ggplot() +
    geom_tile(data = dat_cor, aes(x = FeatureID, y = Phenotype, fill = Rho)) +
    scale_fill_gradientn(colors = c("#2D6DB1", "white", "#DC1623"), limit = c(-1, 1)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    geom_text(data = dat_cor, aes(x = FeatureID, y = Phenotype, label = siglabel), size = 3) +
    geom_point(data = dat_LM_import_final, aes(x = FeatureID, y = Phenotype, size = importance * 100), shape = 1) +
    scale_size_continuous(range = c(0, 5)) +
    labs(y = "", x = "", fill = "Correlation", size = "Importance (%)") +  
    guides(fill = guide_legend(order = 1), 
           size = guide_legend(order = 2)) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          panel.background = element_rect(color = "black"), 
          legend.key = element_blank(), 
          axis.text.x = element_text(color = "black", angle = 45, hjust = 1, vjust = 1), 
          axis.text.y = element_text(color = "black"), 
          axis.ticks = element_line(color = "black"))
  
  
  p2 <- ggplot() +
    geom_col(data = dat_LM_R2, aes(x = FeatureID, y = R2 * 100), fill = "#4882B2", width = 0.6) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
    labs(y = "Explained variation (%)", x = "") +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(color = "black", angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black"),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"))
  
  
  pl <- cowplot::plot_grid(p2, p1, ncol = 2, align = "h", rel_widths = c(0.7, 1))
  
  return(pl)  
}

cowplot::plot_grid(
  get_LM_plot(datCor = gut_corres, datLM = gut_LM),
  get_LM_plot(datCor = oral_corres, datLM = oral_LM),
  ncol = 1, align = "hv", labels = LETTERS[1:2])
```

<img src="13-Statistics_LMFeature_files/figure-html/unnamed-chunk-7-1.png" width="100%" />


结果：*基于相关性和最优多元回归模型的人体特征指标差异对微生物群落差异和微生物门相对丰度差异的贡献。热图表示了人体特征指标和微生物的Spearman相关系数，柱形图表示了人体指标对解释微生物变异的总贡献（通过多元线性回归获得），圆圈大小表示人体指标的重要性（通过多元线性回归和方差分解分析获得）*。

+ A图是肠道微生物和人体特征变量的结果；B图是口腔微生物和人体特征变量的结果；

+ 肠道微生物和口腔微生物与人体特征变量相关性差异较大，人体特征变量对口腔微生物相对丰度差异共享度较大，Liver fat和Sodium对口腔的Bacterodies, Firmicutes和Proteobacteria物种差异和相对丰度差异有较大贡献，而肠道仅Urea_BUN对Actinobacteria物种差异和相对丰度有较大贡献；

+ 人体特征变量对微生物群落解释度差异在肠道和口腔也存在较大差异。在口腔解释度最高的是Proteobacteria，而在肠道则是Actinobacteria；

+ 相比PERMANOVA分析（各个微生物对整体人体特征变量的总变异度，也即单个微生物对人体变量总体的扰动程度），该组合方法不仅计算了各个微生物对核心人体特征变量（前后向逐步回归筛选）的解释度，而且也计算了这些核心变量对该微生物的重要程度（$R^2$分解成重要性打分）。



## Systemic information

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
#>  package                  * version    date (UTC) lib source
#>  abind                      1.4-5      2016-07-21 [1] CRAN (R 4.3.0)
#>  ade4                       1.7-22     2023-02-06 [1] CRAN (R 4.3.0)
#>  ANCOMBC                    2.4.0      2023-10-24 [1] Bioconductor
#>  ape                        5.7-1      2023-03-13 [1] CRAN (R 4.3.0)
#>  backports                  1.4.1      2021-12-13 [1] CRAN (R 4.3.0)
#>  base64enc                  0.1-3      2015-07-28 [1] CRAN (R 4.3.0)
#>  beachmat                   2.18.0     2023-10-24 [1] Bioconductor
#>  beeswarm                   0.4.0      2021-06-01 [1] CRAN (R 4.3.0)
#>  Biobase                    2.62.0     2023-10-24 [1] Bioconductor
#>  BiocGenerics               0.48.1     2023-11-01 [1] Bioconductor
#>  BiocNeighbors              1.20.2     2024-01-07 [1] Bioconductor 3.18 (R 4.3.2)
#>  BiocParallel               1.36.0     2023-10-24 [1] Bioconductor
#>  BiocSingular               1.18.0     2023-10-24 [1] Bioconductor
#>  biomformat                 1.30.0     2023-10-24 [1] Bioconductor
#>  Biostrings                 2.70.2     2024-01-28 [1] Bioconductor 3.18 (R 4.3.2)
#>  bit                        4.0.5      2022-11-15 [1] CRAN (R 4.3.0)
#>  bit64                      4.0.5      2020-08-30 [1] CRAN (R 4.3.0)
#>  bitops                     1.0-7      2021-04-24 [1] CRAN (R 4.3.0)
#>  blob                       1.2.4      2023-03-17 [1] CRAN (R 4.3.0)
#>  bluster                    1.12.0     2023-10-24 [1] Bioconductor
#>  bookdown                   0.37       2023-12-01 [1] CRAN (R 4.3.0)
#>  boot                     * 1.3-28.1   2022-11-22 [1] CRAN (R 4.3.1)
#>  bslib                      0.6.1      2023-11-28 [1] CRAN (R 4.3.0)
#>  cachem                     1.0.8      2023-05-01 [1] CRAN (R 4.3.0)
#>  caret                    * 6.0-94     2023-03-21 [1] CRAN (R 4.3.0)
#>  caTools                    1.18.2     2021-03-28 [1] CRAN (R 4.3.0)
#>  cellranger                 1.1.0      2016-07-27 [1] CRAN (R 4.3.0)
#>  checkmate                  2.3.1      2023-12-04 [1] CRAN (R 4.3.0)
#>  class                      7.3-22     2023-05-03 [1] CRAN (R 4.3.1)
#>  cli                        3.6.2      2023-12-11 [1] CRAN (R 4.3.0)
#>  cluster                    2.1.4      2022-08-22 [1] CRAN (R 4.3.1)
#>  codetools                  0.2-19     2023-02-01 [1] CRAN (R 4.3.1)
#>  colorspace                 2.1-0      2023-01-23 [1] CRAN (R 4.3.0)
#>  corpcor                    1.6.10     2021-09-16 [1] CRAN (R 4.3.0)
#>  cowplot                    1.1.3      2024-01-22 [1] CRAN (R 4.3.2)
#>  crayon                     1.5.2      2022-09-29 [1] CRAN (R 4.3.0)
#>  CVXR                       1.0-12     2024-02-02 [1] CRAN (R 4.3.2)
#>  data.table                 1.15.0     2024-01-30 [1] CRAN (R 4.3.2)
#>  DBI                        1.2.1      2024-01-12 [1] CRAN (R 4.3.0)
#>  DECIPHER                   2.30.0     2023-10-24 [1] Bioconductor
#>  decontam                   1.22.0     2023-10-24 [1] Bioconductor
#>  DelayedArray               0.28.0     2023-10-24 [1] Bioconductor
#>  DelayedMatrixStats         1.24.0     2023-10-24 [1] Bioconductor
#>  DescTools                  0.99.54    2024-02-03 [1] CRAN (R 4.3.2)
#>  DESeq2                     1.42.0     2023-10-24 [1] Bioconductor
#>  devtools                   2.4.5      2022-10-11 [1] CRAN (R 4.3.0)
#>  digest                     0.6.34     2024-01-11 [1] CRAN (R 4.3.0)
#>  DirichletMultinomial       1.44.0     2023-10-24 [1] Bioconductor
#>  doParallel                 1.0.17     2022-02-07 [1] CRAN (R 4.3.0)
#>  doRNG                      1.8.6      2023-01-16 [1] CRAN (R 4.3.0)
#>  downlit                    0.4.3      2023-06-29 [1] CRAN (R 4.3.0)
#>  dplyr                    * 1.1.4      2023-11-17 [1] CRAN (R 4.3.0)
#>  e1071                      1.7-14     2023-12-06 [1] CRAN (R 4.3.0)
#>  ellipsis                   0.3.2      2021-04-29 [1] CRAN (R 4.3.0)
#>  energy                     1.7-11     2022-12-22 [1] CRAN (R 4.3.0)
#>  evaluate                   0.23       2023-11-01 [1] CRAN (R 4.3.0)
#>  Exact                      3.2        2022-09-25 [1] CRAN (R 4.3.0)
#>  expm                       0.999-9    2024-01-11 [1] CRAN (R 4.3.0)
#>  fansi                      1.0.6      2023-12-08 [1] CRAN (R 4.3.0)
#>  farver                     2.1.1      2022-07-06 [1] CRAN (R 4.3.0)
#>  fastmap                    1.1.1      2023-02-24 [1] CRAN (R 4.3.0)
#>  forcats                  * 1.0.0      2023-01-29 [1] CRAN (R 4.3.0)
#>  foreach                    1.5.2      2022-02-02 [1] CRAN (R 4.3.0)
#>  foreign                    0.8-84     2022-12-06 [1] CRAN (R 4.3.1)
#>  Formula                    1.2-5      2023-02-24 [1] CRAN (R 4.3.0)
#>  fs                         1.6.3      2023-07-20 [1] CRAN (R 4.3.0)
#>  future                     1.33.1     2023-12-22 [1] CRAN (R 4.3.0)
#>  future.apply               1.11.1     2023-12-21 [1] CRAN (R 4.3.0)
#>  generics                   0.1.3      2022-07-05 [1] CRAN (R 4.3.0)
#>  GenomeInfoDb               1.38.5     2023-12-28 [1] Bioconductor 3.18 (R 4.3.2)
#>  GenomeInfoDbData           1.2.11     2024-01-24 [1] Bioconductor
#>  GenomicRanges              1.54.1     2023-10-29 [1] Bioconductor
#>  ggbeeswarm                 0.7.2      2023-04-29 [1] CRAN (R 4.3.0)
#>  ggplot2                  * 3.4.4      2023-10-12 [1] CRAN (R 4.3.0)
#>  ggrepel                    0.9.5      2024-01-10 [1] CRAN (R 4.3.0)
#>  gld                        2.6.6      2022-10-23 [1] CRAN (R 4.3.0)
#>  glmnet                     4.1-8      2023-08-22 [1] CRAN (R 4.3.0)
#>  globals                    0.16.2     2022-11-21 [1] CRAN (R 4.3.0)
#>  glue                       1.7.0      2024-01-09 [1] CRAN (R 4.3.0)
#>  gmp                        0.7-4      2024-01-15 [1] CRAN (R 4.3.0)
#>  gower                      1.0.1      2022-12-22 [1] CRAN (R 4.3.0)
#>  gplots                     3.1.3.1    2024-02-02 [1] CRAN (R 4.3.2)
#>  gridExtra                  2.3        2017-09-09 [1] CRAN (R 4.3.0)
#>  gsl                        2.1-8      2023-01-24 [1] CRAN (R 4.3.0)
#>  gtable                     0.3.4      2023-08-21 [1] CRAN (R 4.3.0)
#>  gtools                     3.9.5      2023-11-20 [1] CRAN (R 4.3.0)
#>  hardhat                    1.3.1      2024-02-02 [1] CRAN (R 4.3.2)
#>  highr                      0.10       2022-12-22 [1] CRAN (R 4.3.0)
#>  Hmisc                      5.1-1      2023-09-12 [1] CRAN (R 4.3.0)
#>  hms                        1.1.3      2023-03-21 [1] CRAN (R 4.3.0)
#>  htmlTable                  2.4.2      2023-10-29 [1] CRAN (R 4.3.0)
#>  htmltools                  0.5.7      2023-11-03 [1] CRAN (R 4.3.0)
#>  htmlwidgets                1.6.4      2023-12-06 [1] CRAN (R 4.3.0)
#>  httpuv                     1.6.14     2024-01-26 [1] CRAN (R 4.3.2)
#>  httr                       1.4.7      2023-08-15 [1] CRAN (R 4.3.0)
#>  igraph                     2.0.1.1    2024-01-30 [1] CRAN (R 4.3.2)
#>  ipred                      0.9-14     2023-03-09 [1] CRAN (R 4.3.0)
#>  IRanges                    2.36.0     2023-10-24 [1] Bioconductor
#>  irlba                      2.3.5.1    2022-10-03 [1] CRAN (R 4.3.0)
#>  iterators                  1.0.14     2022-02-05 [1] CRAN (R 4.3.0)
#>  jquerylib                  0.1.4      2021-04-26 [1] CRAN (R 4.3.0)
#>  jsonlite                   1.8.8      2023-12-04 [1] CRAN (R 4.3.0)
#>  KernSmooth                 2.23-21    2023-05-03 [1] CRAN (R 4.3.1)
#>  knitr                      1.45       2023-10-30 [1] CRAN (R 4.3.0)
#>  labeling                   0.4.3      2023-08-29 [1] CRAN (R 4.3.0)
#>  later                      1.3.2      2023-12-06 [1] CRAN (R 4.3.0)
#>  lattice                  * 0.21-8     2023-04-05 [1] CRAN (R 4.3.1)
#>  lava                       1.7.3      2023-11-04 [1] CRAN (R 4.3.0)
#>  lazyeval                   0.2.2      2019-03-15 [1] CRAN (R 4.3.0)
#>  leaps                    * 3.1        2020-01-16 [1] CRAN (R 4.3.0)
#>  lifecycle                  1.0.4      2023-11-07 [1] CRAN (R 4.3.0)
#>  limma                      3.58.1     2023-10-31 [1] Bioconductor
#>  listenv                    0.9.1      2024-01-29 [1] CRAN (R 4.3.2)
#>  lme4                       1.1-35.1   2023-11-05 [1] CRAN (R 4.3.0)
#>  lmerTest                   3.1-3      2020-10-23 [1] CRAN (R 4.3.0)
#>  lmom                       3.0        2023-08-29 [1] CRAN (R 4.3.0)
#>  locfit                     1.5-9.8    2023-06-11 [1] CRAN (R 4.3.0)
#>  lubridate                * 1.9.3      2023-09-27 [1] CRAN (R 4.3.0)
#>  magrittr                   2.0.3      2022-03-30 [1] CRAN (R 4.3.0)
#>  MASS                     * 7.3-60     2023-05-04 [1] CRAN (R 4.3.1)
#>  Matrix                   * 1.6-5      2024-01-11 [1] CRAN (R 4.3.0)
#>  MatrixGenerics             1.14.0     2023-10-24 [1] Bioconductor
#>  matrixStats                1.2.0      2023-12-11 [1] CRAN (R 4.3.0)
#>  memoise                    2.0.1      2021-11-26 [1] CRAN (R 4.3.0)
#>  metagenomeSeq              1.43.0     2023-04-25 [1] Bioconductor
#>  mgcv                       1.8-42     2023-03-02 [1] CRAN (R 4.3.1)
#>  mia                        1.10.0     2023-10-24 [1] Bioconductor
#>  MicrobiomeAnalysis         1.0.3      2024-02-06 [1] Github (HuaZou/MicrobiomeAnalysis@fd2a6a2)
#>  mime                       0.12       2021-09-28 [1] CRAN (R 4.3.0)
#>  miniUI                     0.1.1.1    2018-05-18 [1] CRAN (R 4.3.0)
#>  minqa                      1.2.6      2023-09-11 [1] CRAN (R 4.3.0)
#>  mitools                  * 2.4        2019-04-26 [1] CRAN (R 4.3.0)
#>  ModelMetrics               1.2.2.2    2020-03-17 [1] CRAN (R 4.3.0)
#>  multcomp                   1.4-25     2023-06-20 [1] CRAN (R 4.3.0)
#>  MultiAssayExperiment       1.28.0     2023-10-24 [1] Bioconductor
#>  multtest                   2.58.0     2023-10-24 [1] Bioconductor
#>  munsell                    0.5.0      2018-06-12 [1] CRAN (R 4.3.0)
#>  mvtnorm                    1.2-4      2023-11-27 [1] CRAN (R 4.3.0)
#>  nlme                       3.1-162    2023-01-31 [1] CRAN (R 4.3.1)
#>  nloptr                     2.0.3      2022-05-26 [1] CRAN (R 4.3.0)
#>  nnet                       7.3-19     2023-05-03 [1] CRAN (R 4.3.1)
#>  numDeriv                   2016.8-1.1 2019-06-06 [1] CRAN (R 4.3.0)
#>  parallelly                 1.36.0     2023-05-26 [1] CRAN (R 4.3.0)
#>  permute                    0.9-7      2022-01-27 [1] CRAN (R 4.3.0)
#>  phyloseq                 * 1.46.0     2023-10-24 [1] Bioconductor
#>  pillar                     1.9.0      2023-03-22 [1] CRAN (R 4.3.0)
#>  pkgbuild                   1.4.3      2023-12-10 [1] CRAN (R 4.3.0)
#>  pkgconfig                  2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
#>  pkgload                    1.3.4      2024-01-16 [1] CRAN (R 4.3.0)
#>  plyr                       1.8.9      2023-10-02 [1] CRAN (R 4.3.0)
#>  pROC                       1.18.5     2023-11-01 [1] CRAN (R 4.3.0)
#>  prodlim                    2023.08.28 2023-08-28 [1] CRAN (R 4.3.0)
#>  profvis                    0.3.8      2023-05-02 [1] CRAN (R 4.3.0)
#>  promises                   1.2.1      2023-08-10 [1] CRAN (R 4.3.0)
#>  proxy                      0.4-27     2022-06-09 [1] CRAN (R 4.3.0)
#>  purrr                    * 1.0.2      2023-08-10 [1] CRAN (R 4.3.0)
#>  R6                         2.5.1      2021-08-19 [1] CRAN (R 4.3.0)
#>  rbibutils                  2.2.16     2023-10-25 [1] CRAN (R 4.3.0)
#>  RColorBrewer               1.1-3      2022-04-03 [1] CRAN (R 4.3.0)
#>  Rcpp                       1.0.12     2024-01-09 [1] CRAN (R 4.3.0)
#>  RCurl                      1.98-1.14  2024-01-09 [1] CRAN (R 4.3.0)
#>  Rdpack                     2.6        2023-11-08 [1] CRAN (R 4.3.0)
#>  readr                    * 2.1.5      2024-01-10 [1] CRAN (R 4.3.0)
#>  readxl                     1.4.3      2023-07-06 [1] CRAN (R 4.3.0)
#>  recipes                    1.0.9      2023-12-13 [1] CRAN (R 4.3.0)
#>  relaimpo                 * 2.2-7      2023-10-04 [1] CRAN (R 4.3.0)
#>  remotes                    2.4.2.1    2023-07-18 [1] CRAN (R 4.3.0)
#>  reshape2                   1.4.4      2020-04-09 [1] CRAN (R 4.3.0)
#>  rhdf5                      2.46.1     2023-11-29 [1] Bioconductor
#>  rhdf5filters               1.14.1     2023-11-06 [1] Bioconductor
#>  Rhdf5lib                   1.24.1     2023-12-12 [1] Bioconductor 3.18 (R 4.3.2)
#>  rlang                      1.1.3      2024-01-10 [1] CRAN (R 4.3.0)
#>  rmarkdown                  2.25       2023-09-18 [1] CRAN (R 4.3.0)
#>  Rmpfr                      0.9-5      2024-01-21 [1] CRAN (R 4.3.0)
#>  rngtools                   1.5.2      2021-09-20 [1] CRAN (R 4.3.0)
#>  rootSolve                  1.8.2.4    2023-09-21 [1] CRAN (R 4.3.0)
#>  rpart                      4.1.19     2022-10-21 [1] CRAN (R 4.3.1)
#>  RSQLite                    2.3.5      2024-01-21 [1] CRAN (R 4.3.0)
#>  rstudioapi                 0.15.0     2023-07-07 [1] CRAN (R 4.3.0)
#>  rsvd                       1.0.5      2021-04-16 [1] CRAN (R 4.3.0)
#>  S4Arrays                   1.2.0      2023-10-24 [1] Bioconductor
#>  S4Vectors                  0.40.2     2023-11-23 [1] Bioconductor
#>  sandwich                   3.1-0      2023-12-11 [1] CRAN (R 4.3.0)
#>  sass                       0.4.8      2023-12-06 [1] CRAN (R 4.3.0)
#>  ScaledMatrix               1.10.0     2023-10-24 [1] Bioconductor
#>  scales                     1.3.0      2023-11-28 [1] CRAN (R 4.3.0)
#>  scater                     1.30.1     2023-12-06 [1] Bioconductor
#>  scuttle                    1.12.0     2023-10-24 [1] Bioconductor
#>  sessioninfo                1.2.2      2021-12-06 [1] CRAN (R 4.3.0)
#>  shape                      1.4.6      2021-05-19 [1] CRAN (R 4.3.0)
#>  shiny                      1.8.0      2023-11-17 [1] CRAN (R 4.3.0)
#>  SingleCellExperiment       1.24.0     2023-10-24 [1] Bioconductor
#>  SparseArray                1.2.3      2023-12-25 [1] Bioconductor 3.18 (R 4.3.2)
#>  sparseMatrixStats          1.14.0     2023-10-24 [1] Bioconductor
#>  statmod                    1.5.0      2023-01-06 [1] CRAN (R 4.3.0)
#>  stringi                    1.8.3      2023-12-11 [1] CRAN (R 4.3.0)
#>  stringr                  * 1.5.1      2023-11-14 [1] CRAN (R 4.3.0)
#>  SummarizedExperiment       1.32.0     2023-10-24 [1] Bioconductor
#>  survey                   * 4.2-1      2023-05-03 [1] CRAN (R 4.3.0)
#>  survival                 * 3.5-5      2023-03-12 [1] CRAN (R 4.3.1)
#>  TH.data                    1.1-2      2023-04-17 [1] CRAN (R 4.3.0)
#>  tibble                   * 3.2.1      2023-03-20 [1] CRAN (R 4.3.0)
#>  tidyr                    * 1.3.1      2024-01-24 [1] CRAN (R 4.3.2)
#>  tidyselect                 1.2.0      2022-10-10 [1] CRAN (R 4.3.0)
#>  tidytree                   0.4.6      2023-12-12 [1] CRAN (R 4.3.0)
#>  tidyverse                * 2.0.0      2023-02-22 [1] CRAN (R 4.3.0)
#>  timechange                 0.3.0      2024-01-18 [1] CRAN (R 4.3.0)
#>  timeDate                   4032.109   2023-12-14 [1] CRAN (R 4.3.0)
#>  treeio                     1.26.0     2023-10-24 [1] Bioconductor
#>  TreeSummarizedExperiment   2.10.0     2023-10-24 [1] Bioconductor
#>  tzdb                       0.4.0      2023-05-12 [1] CRAN (R 4.3.0)
#>  urlchecker                 1.0.1      2021-11-30 [1] CRAN (R 4.3.0)
#>  usethis                    2.2.2      2023-07-06 [1] CRAN (R 4.3.0)
#>  utf8                       1.2.4      2023-10-22 [1] CRAN (R 4.3.0)
#>  vctrs                      0.6.5      2023-12-01 [1] CRAN (R 4.3.0)
#>  vegan                      2.6-4      2022-10-11 [1] CRAN (R 4.3.0)
#>  vipor                      0.4.7      2023-12-18 [1] CRAN (R 4.3.0)
#>  viridis                    0.6.5      2024-01-29 [1] CRAN (R 4.3.2)
#>  viridisLite                0.4.2      2023-05-02 [1] CRAN (R 4.3.0)
#>  withr                      3.0.0      2024-01-16 [1] CRAN (R 4.3.0)
#>  Wrench                     1.20.0     2023-10-24 [1] Bioconductor
#>  xfun                       0.41       2023-11-01 [1] CRAN (R 4.3.0)
#>  xml2                       1.3.6      2023-12-04 [1] CRAN (R 4.3.0)
#>  xtable                     1.8-4      2019-04-21 [1] CRAN (R 4.3.0)
#>  XVector                    0.42.0     2023-10-24 [1] Bioconductor
#>  yaml                       2.3.8      2023-12-11 [1] CRAN (R 4.3.0)
#>  yulab.utils                0.1.4      2024-01-28 [1] CRAN (R 4.3.2)
#>  zlibbioc                   1.48.0     2023-10-24 [1] Bioconductor
#>  zoo                        1.8-12     2023-04-13 [1] CRAN (R 4.3.0)
#> 
#>  [1] /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library
#> 
#> ──────────────────────────────────────────────────────────
```


## Reference

+ Balance between community assembly processes mediates species coexistence in agricultural soil microbiomes across eastern China

+ [仿一篇文献的相关性分析和线性模型评估影响群落组成的重要环境变量](https://mp.weixin.qq.com/s/EP63-mi9GJSTRS35dQSgMw)

+ [Stepwise Regression Essentials in R](http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/154-stepwise-regression-essentials-in-r/)
