


# (PART) Metabolomics Data Analysis {.unnumbered}


# Data Processing {#dataprocessing}


代谢组数据一般是搜库后的质谱峰度谱数据，用峰强intensity表示。Raw intensity通常不直接用于假设检验或线性回归等统计方法，需要对其做数据预处理。

本次应用到的数据是Zeybel 2022年发布的文章_Multiomics Analysis Reveals the Impact of Microbiota on Host Metabolism in Hepatic Steatosis_的粪便代谢组学质谱数据。

> 55份粪便代谢组，1032个代谢物

## 处理流程

+ 数据检查：1.核查所有代谢物intensity value是数值型；2.缺失值的比例

+ 补充缺失值

+ 数据过滤（针对代谢物或样本）

+ 数据标准化（针对代谢物或样本）


## 加载R包

```r
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
library(tidyverse)
library(SummarizedExperiment)

# rm(list = ls())
options(stringsAsFactors = F)
options(future.globals.maxSize = 1000 * 1024^2)
```


## 导入数据

对数据[OmicsDataSet-Zeybel et al. - 2022.xlsx](https://github.com/HuaZou/DraftNotes/blob/main/InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx)处理后生成的输入文件，详细情况可参考**Data Set**具体章节。

> ```R
> saveRDS(se_metabolite, "./InputData/result/Zeybel_2022_fecal_metabolite_se.RDS", compress = TRUE)
> ```


```r
data_meta <- readRDS("./InputData/result/Zeybel_2022_fecal_metabolite_se.RDS")

data_meta
#> class: SummarizedExperiment 
#> dim: 1032 55 
#> metadata(0):
#> assays(1): ''
#> rownames(1032): Chem_100002945 Chem_100002356 ...
#>   Chem_100015836 Chem_826
#> rowData names(13): metabolitesID BIOCHEMICAL ... KEGG
#>   SampleIDHMDBID
#> colnames(55): P101001 P101003 ... P101095 P101096
#> colData names(47): PatientID Gender ...
#>   Right_leg_fat_free_mass Right_leg_total_body_water
```


```r
# colData(data_meta)
# assay(data_meta)
# rowData(data_meta)
```


## 构建QC样本

该数据集不存在QC样本，这导致不能做QC的变化范围的数据过滤，因此构建新的QC样本。以上述随机抽取6个样本作为QC样本。

QC样本一般是送测样本混合后再分成N份样本（迈维非靶或广靶均是这样做）再测，它可以评估每次质谱的效果或做代谢物过滤。


```r
meta_tab <- colData(data_meta) |>
  as.data.frame()
feature_tab <- rowData(data_meta) 
assay_tab <- assay(data_meta) |>
  as.data.frame()

rand_sample <- sample(colnames(assay_tab), 6)
QC_assay <- assay_tab[, rand_sample]
colnames(QC_assay) <- paste0("QC", 1:6)
assay_tab_new <- cbind(assay_tab, QC_assay)

QC_meta <- data.frame(matrix(NA, nrow = 6, ncol = ncol(meta_tab))) 
colnames(QC_meta) <- colnames(meta_tab)
rownames(QC_meta) <- colnames(QC_assay)
QC_meta$LiverFatClass <- "QC"
QC_meta$PatientID <- rownames(QC_meta)
meta_tab_new <- rbind(meta_tab, QC_meta)

data_meta_new <- SummarizedExperiment(
  assays = assay_tab_new,
  rowData = feature_tab,
  colData = meta_tab_new,
  checkDimnames = TRUE)

data_meta_new
#> class: SummarizedExperiment 
#> dim: 1032 61 
#> metadata(0):
#> assays(1): ''
#> rownames(1032): Chem_100002945 Chem_100002356 ...
#>   Chem_100015836 Chem_826
#> rowData names(13): metabolitesID BIOCHEMICAL ... KEGG
#>   SampleIDHMDBID
#> colnames(61): P101001 P101003 ... QC5 QC6
#> colData names(47): PatientID Gender ...
#>   Right_leg_fat_free_mass Right_leg_total_body_water
```


## 数据过滤


```r
CheckData <- function(object) {
  
  # object = data_meta_new
  
  # features are in rows and Samples in columns
  DataAssay <- SummarizedExperiment::assay(object)
  
  # numeric & missing values
  int_mat <- DataAssay
  rowNms <- rownames(int_mat)
  colNms <- colnames(int_mat)
  naNms <- sum(is.na(int_mat))
  for (i in 1:ncol(int_mat)) {
    if (class(int_mat[, i]) == "integer64") {
      int_mat[, i] <- as.double(int_mat[, i])
    }
  }
  
  num_mat <- apply(int_mat, 2, as.numeric)
  if (sum(is.na(num_mat)) > naNms) {
    num_mat <- apply(int_mat, 2, function(x) as.numeric(gsub(",",  "", x)))
    if (sum(is.na(num_mat)) > naNms) {
      message("<font color=\"red\">Non-numeric values were found and replaced by NA.</font>")
    } else {
      message("All data values are numeric.")
    }
  } else {
    message("All data values are numeric.")
  }
  
  int_mat <- num_mat
  rownames(int_mat) <- rowNms
  colnames(int_mat) <- colNms
  varCol <- apply(int_mat, 2, var, na.rm = T)
  constCol <- (varCol == 0 | is.na(varCol))
  constNum <- sum(constCol, na.rm = T)
  if (constNum > 0) {
    print(paste("<font color=\"red\">", constNum, 
      "features with a constant or single value across samples were found and deleted.</font>"))
    int_mat <- int_mat[, !constCol, drop = FALSE]
  }
  
  totalCount <- nrow(int_mat) * ncol(int_mat)
  naCount <- sum(is.na(int_mat))
  naPercent <- round(100 * naCount/totalCount, 1)

  print(paste("A total of ", naCount, " (", naPercent, 
    "%) missing values were detected.", sep = ""))  

  DataMeta <- colData(object) |>
    as.data.frame()
  DataFeature <- rowData(object)
  
  res <- SummarizedExperiment(
    assays = int_mat,
    rowData = DataFeature,
    colData = DataMeta,
    checkDimnames = TRUE)
  
  return(object)
}

se_check <- CheckData(object = data_meta_new)
#> [1] "A total of 7846 (12.5%) missing values were detected."
se_check
#> class: SummarizedExperiment 
#> dim: 1032 61 
#> metadata(0):
#> assays(1): ''
#> rownames(1032): Chem_100002945 Chem_100002356 ...
#>   Chem_100015836 Chem_826
#> rowData names(13): metabolitesID BIOCHEMICAL ... KEGG
#>   SampleIDHMDBID
#> colnames(61): P101001 P101003 ... QC5 QC6
#> colData names(47): PatientID Gender ...
#>   Right_leg_fat_free_mass Right_leg_total_body_water
```

结果：12.5%的缺失值存在，下面进行缺失值补充。


## 补缺失值

Missing Value 的产生原因主要有两个：1）一个代谢峰在某些生物样品中存在而在另外一些生物样品中不存在； 2）某些代谢物在生物样品中的浓度低于质谱的检测限。

Missing Value 在数据中的表现形式为 NA值。对于大规模代谢组学来说，因为其长时间的数据采集，质谱灵敏度的漂移使MV的问题更加严重。一般来说，对一个代谢组学数据来说， Missing Value 会占到所有数据点的 20%左右 。首先需要对 Missing Value 进行过滤，对于 Missing Value超过一定比例的代谢峰来说，该代谢峰很有可能是一个偶然出现的噪声信号，因此将其从数据中删除。 比如，在代谢组学数据中，一般采用的标准为代谢峰需要在 80%的质量控制（quality control， QC） 样品中出现，否则删除。对于 Missing Value 超过一定比例的生物样品来说，该样品很有可能是在样品制备或者数据采集过程中出现了误差，如该样品稀释比例异常或者进样体积异常，这些样品需要被删除掉。

删除掉异常的生物样品和代谢峰之后，剩余的 Missing Value 需要统计学方法进行补齐（MVimputation），不同的 Missing Value补齐方法对数据的结构影响非常大，最好的 Missing Value 补齐方法是那些可以最好的重构出数据原本结构的方法。 Gromski (The influence of scaling metabolomics data on model classification accuracy) 通过使用完整的代谢组学数据构建 Missing Value 数据，然后使用不同的 Missing Value补齐方法对数据进行补齐，然后对不同 Missing Value 补齐方法补齐的数据进行多元统计学分析， 最终发现 K 值临近方法（K-nearest neighbor，KNN）对代谢组学数据的补齐效果最好。

缺失值补充方法有很多种，如下

+ “none”: all missing values will be replaced by zero.

+ “LOD”: specific Limit Of Detection which provides by user.

+ “half_min”: half minimal values across samples except zero.

+ “median”: median values across samples except zero.

+ “mean”: mean values across samples except zero.

+ “min”: minimal values across samples except zero.

+ “knn”: k-nearest neighbors samples.

+ “rf”: nonparametric missing value imputation using Random Forest.

+ “QRILC”: missing values imputation based quantile regression. (default: “none”).


一般采用k近邻的方法，它的原理是离该缺失值样本最近的K个样本具有类似的属性，使用它们的平均值填补缺失值是相对可靠的方法，但是也需要注意该方法的阈值适用范围。

这里使用**[MicrobiomeAnalysis](https://zouhua.top/MicrobiomeAnalysis/)**提供的[impute_abundance](https://zouhua.top/MicrobiomeAnalysis/reference/impute_abundance.html)函数，先安装此包。
```R
if (!requireNamespace(c("remotes", "devtools"), quietly=TRUE)) {
  install.packages(c("devtools", "remotes"))
}

remotes::install_github("HuaZou/MicrobiomeAnalysis")

# library(MicrobiomeAnalysis)
```


```r
library(MicrobiomeAnalysis)

se_impute <- impute_abundance(
  object = se_check,
  group = "LiverFatClass",
  method = "knn",
  cutoff = 50,
  knum = 10)

se_impute
#> class: SummarizedExperiment 
#> dim: 957 61 
#> metadata(0):
#> assays(1): ''
#> rownames(957): Chem_100002945 Chem_100002356 ...
#>   Chem_100015836 Chem_826
#> rowData names(13): metabolitesID BIOCHEMICAL ... KEGG
#>   SampleIDHMDBID
#> colnames(61): P101001 P101003 ... QC5 QC6
#> colData names(47): PatientID Gender ...
#>   Right_leg_fat_free_mass Right_leg_total_body_water
```

结果：代谢物的缺失值在任何组大于50%会被移除，最后移除的代谢物表达矩阵用于缺失值补充。


## 数据过滤

在非靶向代谢组或蛋白质组经常会使用该方法，目的是过滤掉不太可能用于分析的代谢物或蛋白质。

过滤会基于QC样本的相对丰度标准方差relative standard deviation (RSD = SD / mean)，可以理解为特征的数据波动范围。 LC-MS或GC-MS对不同样本可能存在不同偏好行，采用QC样本可以得到波动范围，那些具有高RSD的特征可能受到质谱操作的影响较大，因此它们的可靠性相对较低不需要用于后续下游分析。一般情况下，RSD阈值在LC-MS和GC-MS分别是20%和30%（*保留波动小于该阈值的代谢物*）。

在通过QC的RSD过滤完后，还可以通过以下方法过滤噪声（过滤低丰度或高变化的特征）

+ 过滤方法

  - Interquantile range (IQR)（过滤常数特征，即波动较小的特征，它们可能是常态表达）;
  
  - Standard deviation (SD) （过滤常数特征，即波动较小的特征，它们可能是常态表达）;
  
  - Median absolute deviation (MAD) （过滤常数特征，即波动较小的特征，它们可能是常态表达）;
  
  - Relative standard deviation (RSD = SD/mean) （过滤低重复性特征，它们可能受到质谱操作影响）;
  
  - Non-parametric relative standard deviation (MAD/median) （过滤低重复性特征，它们可能受到质谱操作影响）;
  
  - Mean intensity value （过滤低丰度特征，它们可能是噪声或仪器测量极限值）;
  
  - Median intensity value （过滤低丰度特征，它们可能是噪声或仪器测量极限值）;


+ 一般过滤的阈值设置

  - 少于 250 个特征: 5%
  
  - 介于 250 到 500 个特征: 10%
  
  - 介于 500 到 1000 个特征: 25%
  
  - 超过 1000 个特征: 40%
  


```r
FilterFeature <- function(
    object,
    group,    
    qc_label,
    method = c("none", "iqr", "rsd", 
               "nrsd", "mean", "sd",
               "mad", "median"),
    rsd_cutoff = 25) {
  
  # object = se_impute
  # group = "LiverFatClass"  
  # qc_label = "QC"
  # method = "iqr"
  # rsd_cutoff = 25  
  
  # row->features; col->samples  
  features_tab <- SummarizedExperiment::assay(object) 
  metadata_tab <- SummarizedExperiment::colData(object) 
  
  # QC samples
  colnames(metadata_tab)[which(colnames(metadata_tab) == group)] <- "TempGroup"
  qc_samples <- metadata_tab %>% 
    as.data.frame() %>%
    dplyr::filter(TempGroup == qc_label)
  if (dim(qc_samples)[1] == 0) {
    stop("No qc samples have been chosen, please check your input")
  }
  
  # QC samples' feature table
  qc_feature <- features_tab[, colnames(features_tab) %in% 
                               rownames(qc_samples)] %>%
    t()
  
  # filter features by QC RSD
  rsd <- rsd_cutoff / 100
  sds <- apply(qc_feature, 2, sd, na.rm = T)
  mns <- apply(qc_feature, 2, mean, na.rm = T)
  rsd_vals <- abs(sds/mns) %>% na.omit()
  gd_inx <- rsd_vals < rsd
  int_mat <- features_tab[gd_inx, ]
  print(paste("Removed ", (dim(qc_feature)[2] - dim(int_mat)[1]), 
  " features based on QC RSD values. QC samples are excluded from downstream functional analysis."))
  
  # whether to filter features by percentage according to the number
  PerformFeatureFilter <- function(datMatrix, 
                                   qc_method = method,
                                   remain_num = NULL) {
    
    dat <- datMatrix
    feat_num <- ncol(dat)
    feat_nms <- colnames(dat)
    nm <- NULL
    if (qc_method == "none" && feat_num < 5000) { # only allow for less than 4000
      remain <- rep(TRUE, feat_num)
      nm <- "No filtering was applied"
    } else {
      if (qc_method == "rsd"){
        sds <- apply(dat, 2, sd, na.rm = T)
        mns <- apply(dat, 2, mean, na.rm = T)
        filter_val <- abs(sds/mns)
        nm <- "Relative standard deviation"
      } else if (qc_method == "nrsd" ) {
        mads <- apply(dat, 2, mad, na.rm = T)
        meds <- apply(dat, 2, median, na.rm = T)
        filter_val <- abs(mads/meds)
        nm <- "Non-paramatric relative standard deviation"
      } else if (qc_method == "mean") {
        filter_val <- apply(dat, 2, mean, na.rm = T)
        nm <- "mean"
      } else if (qc_method == "sd") {
        filter_val <- apply(dat, 2, sd, na.rm = T)
        nm <- "standard deviation"
      } else if (qc_method == "mad") {
        filter_val <- apply(dat, 2, mad, na.rm = T)
        nm <- "Median absolute deviation"
      } else if (qc_method == "median") {
        filter_val <- apply(dat, 2, median, na.rm = T)
        nm <- "median"
      } else if (qc_method == "iqr") { # iqr
        filter_val <- apply(dat, 2, IQR, na.rm = T)
        nm <- "Interquantile Range"
      }
      
      # get the rank of the filtered variables
      rk <- rank(-filter_val, ties.method = "random")
      
      if (is.null(remain_num)) { # apply empirical filtering based on data size
          if (feat_num < 250) { # reduce 5%
            remain <- rk < feat_num * 0.95
            message("Further feature filtering based on ", nm)
          } else if (feat_num < 500) { # reduce 10%
            remain <- rk < feat_num * 0.9
            message("Further feature filtering based on ", nm)
          } else if (feat_num < 1000) { # reduce 25%
            remain <- rk < feat_num * 0.75
            message("Further feature filtering based on ", nm)
          } else { # reduce 40%, if still over 5000, then only use top 5000
            remain <- rk < feat_num * 0.6
            message("Further feature filtering based on ", nm)
          }
      } else {
        remain <- rk < remain_num
      }
    }
    
    res <- datMatrix[, remain]
    
    return(res)
  }  
  
  feature_res <- PerformFeatureFilter(t(int_mat))
  
  # remove QC samples 
  feature_final <- feature_res[!rownames(feature_res) %in% rownames(qc_samples), ]
  
  # save int_mat into se object 
  datarow <- object@elementMetadata %>% 
    as.data.frame() 
  rownames(datarow) <- datarow$metabolitesID
  res <- import_SE(
    object = t(feature_final),
    rowdata = datarow,
    coldata = object@colData)
  
  return(res) 
}

se_filter <- FilterFeature(
  object = se_impute,
  group = "LiverFatClass",
  qc_label = "QC",
  method = "iqr",
  rsd_cutoff = 90)
#> [1] "Removed  69  features based on QC RSD values. QC samples are excluded from downstream functional analysis."

se_filter
#> class: SummarizedExperiment 
#> dim: 665 55 
#> metadata(0):
#> assays(1): ''
#> rownames(665): Chem_100002945 Chem_100002356 ...
#>   Chem_1004 Chem_100015836
#> rowData names(13): metabolitesID BIOCHEMICAL ... KEGG
#>   SampleIDHMDBID
#> colnames(55): P101001 P101003 ... P101095 P101096
#> colData names(47): PatientID Gender ...
#>   Right_leg_fat_free_mass Right_leg_total_body_water
```

**结果**：过滤掉不符合要求的代谢物以及也过滤掉了QC样本 (25%过滤太多了，这里选择90%)

+ 根据QC样本的代谢物RSD过滤代谢物

+ 再根据代谢物波动范围过滤不符合的代谢物


## 数据标准化

数据标准化为了三部分：

+ 基于单个数值本身的数据转换，目的是将数据进行各种形式的转换从而提高数据的正态分布性，校正奇异值，达到减少分析误差的效果。代谢组学数据大都为偏倚分布，因此数据转换是非常常见的数据处理方式之一，常用的方法为对数变换（log）。经过数据转换之后的数据仍然没有处在同一个标准量度上，对于代谢组学数据来说，不同代谢峰的强度可以相差几个甚至十几个数量级，如此大的差别导致在多元统计分析时，强度大的代谢峰有可能掩盖强度小的代谢峰的贡献。因此，在统计分析前，尤其是多元统计分析，为了将所有代谢峰统一到同一个量度，需要对数据进行中心化和标度化。

  - Log transformation (base 10)

  - Square root transformation (square root of data values)

  - Cube root transformation (cube root of data values)

使用**MicrobiomeAnalysis**提供的`transform_abundances`，选择*log10p*，改变数据的偏态分布。

```r
se_tran <- MicrobiomeAnalysis::transform_abundances(
  object = se_filter,
  transform = "log10p")

se_tran
#> class: SummarizedExperiment 
#> dim: 665 55 
#> metadata(0):
#> assays(1): ''
#> rownames(665): Chem_100002945 Chem_100002356 ...
#>   Chem_1004 Chem_100015836
#> rowData names(13): metabolitesID BIOCHEMICAL ... KEGG
#>   SampleIDHMDBID
#> colnames(55): P101001 P101003 ... P101095 P101096
#> colData names(47): PatientID Gender ...
#>   Right_leg_fat_free_mass Right_leg_total_body_water
```


+ 基于样本自身的标准化，目的是去除样本间的系统差异（比如不同测序深度），例如通过均值中心化处理， 代谢峰转变为与自己平均值之间的差值，且所有的变量都以零为中心变化，因此中心化的数据就能直接反应变量的变化情况，有利于观察组间差异和聚类分析。

  - Sample-specific normalization (i.e. weight, volume)
  
  - Normalization by sum (relative abundance)
  
  - Normalization by median
  
  - Normalization by a reference sample (PQN)
  
  - Normalization by a pooled sample from group (group PQN)
  
  - Normalization by reference feature
  
  - Quantile normalization (suggested only for > 1000 features)


```r
NormalizeData <- function(
    object,
    rowNorm = c("Quantile", "GroupPQN", "SamplePQN",
                "CompNorm", "SumNorm", "MedianNorm",
                "SpecNorm", "None"),
    ref = NULL,
    SpeWeight = 1) {
  
  # object = se_tran
  # rowNorm = "SumNorm"
  # ref = NULL
  # SpeWeight = 1
  
  # row->features; col->samples 
  features_tab <- SummarizedExperiment::assay(object) 
  metadata_tab <- SummarizedExperiment::colData(object)   
  
  # row->samples; col->features 
  feaTab <- t(features_tab)
  colNames <- colnames(feaTab)
  rowNames <- rownames(feaTab)
  
  #############################################
  # Sample normalization
  # perform quantile normalization on the raw data (can be log transformed later by user)
  QuantileNormalize <- function(data) {
    return(t(preprocessCore::normalize.quantiles(t(data), copy=FALSE)));
  }
  # normalize by a reference sample (probability quotient normalization)
  # ref should be the name of the reference sample
  ProbNorm <- function(x, ref_smpl) {
    return(x/median(as.numeric(x/ref_smpl), na.rm = T))
  }
  
  # normalize by a reference reference (i.e. creatinine)
  # ref should be the name of the cmpd
  CompNorm <- function(x, ref) {
    return(1000 * x/x[ref])
  }
  
  # normalize by sum (relative abundance)
  SumNorm <- function(x) {
    #return(1000 * x/sum(x, na.rm = T))
    return(x/sum(x, na.rm = T))
  }
  
  # normalize by median
  MedianNorm <- function(x) {
    return(x/median(x, na.rm = T))
  }  
  
  # row-wise normalization (samples)
  if (rowNorm == "Quantile") {
    datrowNorm <- QuantileNormalize(feaTab)
    # this can introduce constant variables if a variable is 
    # at the same rank across all samples (replaced by its average across all)
    varCol <- apply(datrowNorm, 2, var, na.rm = T)
    constCol <- (varCol == 0 | is.na(varCol))
    constNum <- sum(constCol, na.rm = T)
    if (constNum > 0) {
      message(paste("After quantile normalization", constNum, 
                    "features with a constant value were found and deleted."))
      datrowNorm <- datrowNorm[, !constCol, drop = FALSE]
      colNames <- colnames(datrowNorm)
      rowNames <- rownames(datrowNorm)
    }
    rownm <- "Quantile Normalization"
  } else if (rowNorm == "GroupPQN") {
    grp_inx <- metadata_tab$group == ref
    ref.smpl <- apply(feaTab[grp_inx, , drop = FALSE], 2, mean)
    datrowNorm <- t(apply(feaTab, 1, ProbNorm, ref.smpl))
    rownm <- "Probabilistic Quotient Normalization by a reference group"
  } else if (rowNorm == "SamplePQN") {
    ref.smpl <- feaTab[ref, , drop = FALSE]
    datrowNorm <- t(apply(feaTab, 1, ProbNorm, ref.smpl))
    rownm <- "Probabilistic Quotient Normalization by a reference sample"
  } else if (rowNorm == "CompNorm") {
    datrowNorm <- t(apply(t(feaTab), 1, CompNorm, ref))
    rownm <- "Normalization by a reference feature";
  } else if (rowNorm == "SumNorm") {
    datrowNorm <- t(apply(feaTab, 1, SumNorm))
    rownm <- "Normalization to constant sum"
  } else if (rowNorm == "MedianNorm") {
    datrowNorm <- t(apply(feaTab, 1, MedianNorm))
    rownm <- "Normalization to sample median"
  } else if(rowNorm == "SpecNorm") {
    norm.vec <- rep(SpeWeight, nrow(feaTab)) # default all same weight vec to prevent error
    datrowNorm <- feaTab / norm.vec
    message("No sample specific information were given, all set to 1.0")
    rownm <- "Normalization by sample-specific factor"
  } else {
    # nothing to do
    rownm <- "N/A"
    datrowNorm <- feaTab
  }
  ################################################ 
  
  # use apply will lose dimension info (i.e. row names and colnames)
  # row->samples; col->features 
  rownames(datrowNorm) <- rowNames
  colnames(datrowNorm) <- colNames
  
  # if the reference by feature, the feature column should be removed, since it is all 1
  if(rowNorm == "CompNorm" && !is.null(ref)){
    inx <- match(ref, colnames(datrowNorm))
    datrowNorm <- datrowNorm[, -inx, drop=FALSE]
    colNames <- colNames[-inx]
  }

  DataMeta <- colData(object) |>
    as.data.frame()
  DataFeature <- rowData(object)
  
  datrowNorm[is.na(datrowNorm)] <- 0
  
  se <- SummarizedExperiment(
    assays = t(datrowNorm),
    rowData = DataFeature,
    colData = DataMeta,
    checkDimnames = TRUE)  
  
  # need to do some sanity check, for log there may be Inf values introduced
  res <- CheckData(se)
  
  return(res)
}

se_norm <- NormalizeData(
  object = se_tran,
  rowNorm = "SumNorm")
#> [1] "A total of 0 (0%) missing values were detected."

se_norm
#> class: SummarizedExperiment 
#> dim: 665 55 
#> metadata(0):
#> assays(1): ''
#> rownames(665): Chem_100002945 Chem_100002356 ...
#>   Chem_1004 Chem_100015836
#> rowData names(13): metabolitesID BIOCHEMICAL ... KEGG
#>   SampleIDHMDBID
#> colnames(55): P101001 P101003 ... P101095 P101096
#> colData names(47): PatientID Gender ...
#>   Right_leg_fat_free_mass Right_leg_total_body_water
```


+ 基于特征的数据标准化，目的是增加各个变量在不同样品中的可比性

  - Mean centering (mean-centered only)

  - Auto scaling (mean-centered and divided by the standard deviation of each variable)

  - Pareto scaling (mean-centered and divided by the square root of the standard deviation of each variable)

  - Range scaling (mean-centered and divided by the range of each variable)


```r
se_scale <- scale_variables(
  object = se_norm,
  method = "zscore")

se_scale
#> class: SummarizedExperiment 
#> dim: 665 55 
#> metadata(0):
#> assays(1): ''
#> rownames(665): Chem_100002945 Chem_100002356 ...
#>   Chem_1004 Chem_100015836
#> rowData names(13): metabolitesID BIOCHEMICAL ... KEGG
#>   SampleIDHMDBID
#> colnames(55): P101001 P101003 ... P101095 P101096
#> colData names(47): PatientID Gender ...
#>   Right_leg_fat_free_mass Right_leg_total_body_water
```

查看数据状态

```r
SummarizedExperiment::assay(se_filter)  %>% data.frame() %>% head()
#>                   P101001     P101003     P101004
#> Chem_100002945 51127588.0  42040432.0  34940596.0
#> Chem_100002356  5105020.5   4006120.2   3885477.0
#> Chem_100021502   756686.2    983889.2    851026.5
#> Chem_100020519 38451736.2  44188252.3  39133949.7
#> Chem_100008903 94392176.0 117463144.0 115155104.0
#> Chem_100009217   144048.2    246049.9    217577.4
#>                   P101007     P101009    P101010    P101011
#> Chem_100002945 58518636.0  51118832.0 83783688.0 29017984.0
#> Chem_100002356  4285129.5   6665653.5  9057441.0  2802655.2
#> Chem_100021502   726593.9    232959.5   650261.1   541954.8
#> Chem_100020519 45607350.7    421584.8 45395682.9 37810261.1
#> Chem_100008903 79582632.0 118408760.0 92508664.0 94076424.0
#> Chem_100009217   211262.5    431295.3   135226.5   230266.2
#>                   P101012    P101013    P101016     P101017
#> Chem_100002945 51222064.0 77550128.0 30949554.0  26923596.0
#> Chem_100002356  5996555.0 11367511.0  3874736.8   2817151.0
#> Chem_100021502   598491.0   438885.6  1625844.8    566466.9
#> Chem_100020519 41630497.4 38062663.5 37906237.4    134190.9
#> Chem_100008903 69473744.0 80567352.0 98766592.0 129547656.0
#> Chem_100009217   163974.4   162897.1   276264.6    375196.6
#>                    P101018     P101019     P101021
#> Chem_100002945  56720032.0 27956064.00  48723600.0
#> Chem_100002356   8029728.0  3766663.75   5174967.0
#> Chem_100021502    427850.6   519559.00   1301591.2
#> Chem_100020519  42499104.2    76295.75  42311221.7
#> Chem_100008903 118271584.0 37880820.00 163868720.0
#> Chem_100009217    538405.4   197953.11    242844.1
#>                    P101022    P101024     P101025
#> Chem_100002945  16282054.0 77028824.0  32022342.0
#> Chem_100002356   1746182.9  5519105.5   2557365.0
#> Chem_100021502   1474247.4   970475.8    628680.1
#> Chem_100020519    227035.0 46791052.6    105729.3
#> Chem_100008903 106189520.0 71475920.0 147776592.0
#> Chem_100009217    447587.9   145823.4   1936915.5
#>                    P101027     P101030     P101031
#> Chem_100002945  22589448.0  38449788.0 59134052.00
#> Chem_100002356   1882902.6   2860324.2  4721201.00
#> Chem_100021502    635516.6    367246.8   512037.88
#> Chem_100020519  34986824.8    563855.2    73751.73
#> Chem_100008903 127571792.0 128319128.0 90447848.00
#> Chem_100009217    322639.6    381134.2   180923.18
#>                   P101038      P101041     P101042
#> Chem_100002945 32038030.0  20833830.00 33809080.00
#> Chem_100002356  4011627.5   2938779.00  3017260.50
#> Chem_100021502   852000.1    634488.56  1680135.75
#> Chem_100020519 39731812.0     69179.21    54615.88
#> Chem_100008903 46622592.0 111919096.00 89762056.00
#> Chem_100009217   165989.5    288088.38   143925.16
#>                   P101047     P101050    P101051    P101052
#> Chem_100002945 18637508.0  21978476.0 24265162.0 52203780.0
#> Chem_100002356  1935144.1   2897211.0  2476279.5  5928454.0
#> Chem_100021502   326005.6    316650.2   737202.9   459385.9
#> Chem_100020519 32650828.9    114410.1 37838221.5    46874.5
#> Chem_100008903 97617984.0 112900000.0 71779912.0 64008512.0
#> Chem_100009217   112179.6    222175.4   171262.0   247185.3
#>                     P101054    P101056     P101057
#> Chem_100002945  12836384.00 18546636.0 32301820.00
#> Chem_100002356   1685760.62  1650011.0  3419157.75
#> Chem_100021502    346176.59   585470.2   417958.47
#> Chem_100020519     30365.21 38779954.0    71663.52
#> Chem_100008903 111279888.00 68771592.0 77140264.00
#> Chem_100009217    371306.81   208351.8   254912.39
#>                    P101059    P101061    P101062
#> Chem_100002945  22645984.0 23683254.0 29027646.0
#> Chem_100002356   2196044.2  3217499.2  4060367.8
#> Chem_100021502    734586.3   337035.7   982299.4
#> Chem_100020519  36577523.8    76462.0   381492.1
#> Chem_100008903 113564872.0 96143304.0 98940424.0
#> Chem_100009217    416259.6   210019.8   325386.6
#>                    P101064    P101065    P101067    P101068
#> Chem_100002945  32629048.0 22950806.0 33555116.0 44283972.0
#> Chem_100002356   3031529.5  2467147.2  3567913.2  6525382.0
#> Chem_100021502   1255148.0   637699.1   284516.1   664800.1
#> Chem_100020519  45337079.7   284125.7 36151920.0    62530.2
#> Chem_100008903 108473368.0 86418592.0 95236288.0 71289048.0
#> Chem_100009217    380257.6   217008.5   149385.6   121522.4
#>                    P101069     P101071     P101072
#> Chem_100002945 52685972.00  32415040.0  34170948.0
#> Chem_100002356  3984333.50   3001414.5   4679519.0
#> Chem_100021502   684813.56    596846.1    316855.0
#> Chem_100020519 43167938.82  39401928.6    679952.5
#> Chem_100008903 74526792.00 115519872.0 127401592.0
#> Chem_100009217    99392.66    301355.7    590274.3
#>                    P101074    P101075    P101076
#> Chem_100002945  22550616.0 22058076.0 24455466.0
#> Chem_100002356   2529255.5  2583265.5  3515218.2
#> Chem_100021502    646136.8   198381.7   255897.7
#> Chem_100020519    794265.5 34178910.4 34230397.0
#> Chem_100008903 108255040.0 83989120.0 77315256.0
#> Chem_100009217    199429.7   179174.6   204990.0
#>                    P101077    P101079      P101080
#> Chem_100002945  25225170.0 15718590.0  29120336.00
#> Chem_100002356   3272875.0  2449462.5   2695001.50
#> Chem_100021502    547243.4   508791.6   1256550.25
#> Chem_100020519  36913595.2 34413781.2     21754.81
#> Chem_100008903 158257952.0 78587928.0 163246832.00
#> Chem_100009217    767403.1   375777.7    514687.62
#>                    P101081     P101082    P101084
#> Chem_100002945  65904836.0  22908578.0 29140440.0
#> Chem_100002356   6474709.5   2110243.8  3648091.2
#> Chem_100021502    339909.3    596292.2   497300.8
#> Chem_100020519  43515920.2  37389441.4 41948877.2
#> Chem_100008903 124678488.0 100435064.0 86139200.0
#> Chem_100009217    360975.7    197618.8   293438.9
#>                     P101085     P101088     P101090
#> Chem_100002945  20427124.00  29199012.0  24042020.0
#> Chem_100002356   3253531.75   4154170.8   2396959.8
#> Chem_100021502    309859.31    601515.1    794206.0
#> Chem_100020519     33125.69  33291433.8  44136149.0
#> Chem_100008903 103513520.00 101921248.0 107571936.0
#> Chem_100009217    238898.77    347081.2    352560.8
#>                    P101094    P101095    P101096
#> Chem_100002945 36910084.00 35662068.0 66402192.0
#> Chem_100002356  4759584.50  3452283.2  6374383.0
#> Chem_100021502   414972.84  3606340.5  1077637.5
#> Chem_100020519    62620.58 53547228.3   118264.1
#> Chem_100008903 85426888.00 53107852.0 80095704.0
#> Chem_100009217   168247.59   202723.2   121388.2
```


```r
SummarizedExperiment::assay(se_scale) %>% as.data.frame() %>% head()
#>                    P101001    P101003    P101004
#> Chem_100002945  1.07306917 0.86152338  0.5082323
#> Chem_100002356  0.82064850 0.46673458  0.4891795
#> Chem_100021502  0.38348742 1.03109633  0.8277830
#> Chem_100020519  0.80168130 0.88582795  0.8603053
#> Chem_100008903 -0.09730152 1.00326048  1.0979398
#> Chem_100009217 -1.12521659 0.07573911 -0.1028403
#>                    P101007    P101009     P101010
#> Chem_100002945  1.04933454  1.1822222  2.48445942
#> Chem_100002356  0.09580821  1.5621814  2.35582015
#> Chem_100021502  0.10697148 -1.7382319  0.19388958
#> Chem_100020519  0.81114382 -0.7005463  0.87897155
#> Chem_100008903 -1.13495047  0.7922921  0.05627634
#> Chem_100009217 -0.57017590  1.0807971 -1.16419173
#>                    P101011    P101012    P101013   P101016
#> Chem_100002945 -0.28405225  0.9148960  1.9319939 0.2751454
#> Chem_100002356 -0.57872637  1.0734888  2.6014330 0.5460211
#> Chem_100021502 -0.19478619 -0.1396695 -0.7230786 2.0808419
#> Chem_100020519  0.80531889  0.8068771  0.7739970 0.8596250
#> Chem_100008903 -0.01310616 -1.3136952 -0.8614738 0.6956945
#> Chem_100009217 -0.16992464 -0.9574710 -0.9821566 0.4104843
#>                   P101017    P101018    P101019    P101021
#> Chem_100002945 -0.7433067  0.9513809 -0.1129977  0.9228859
#> Chem_100002356 -0.7926733  1.5917109  0.3663488  0.8293007
#> Chem_100021502 -0.2666381 -0.8747766 -0.1235194  1.3692197
#> Chem_100020519 -1.1135092  0.7850518 -1.2550937  0.8299187
#> Chem_100008903  0.6506740  0.1219791 -2.6397483  1.6647504
#> Chem_100009217  0.6284458  1.2323728 -0.3195695 -0.1199946
#>                   P101022    P101024    P101025     P101027
#> Chem_100002945 -1.5728003  2.1789558 -0.6821484 -1.06861704
#> Chem_100002356 -1.5894014  1.0713201 -1.3341558 -1.67254046
#> Chem_100021502  1.7636231  0.8848582 -0.2867352  0.01129551
#> Chem_100020519 -0.8954834  0.8770938 -1.2263261  0.75833363
#> Chem_100008903  0.6260545 -0.9099553  0.5513188  0.76439961
#> Chem_100009217  1.2337086 -1.0635292  3.5792427  0.40137905
#>                   P101030    P101031     P101038    P101041
#> Chem_100002945  0.2926796  1.1997769 -0.29628072 -1.0340575
#> Chem_100002356 -0.6382651  0.4323314  0.07147409 -0.3905631
#> Chem_100021502 -0.9868480 -0.4676701  0.49089865  0.1467633
#> Chem_100020519 -0.6200918 -1.3155593  0.78665802 -1.3035970
#> Chem_100008903  0.8184063 -0.5549997 -2.65763306  0.6770235
#> Chem_100009217  0.7408365 -0.8054754 -0.95222048  0.3181630
#>                   P101042    P101047    P101050    P101051
#> Chem_100002945 -0.0958376 -0.7714505 -0.9031042 -0.8393171
#> Chem_100002356 -0.5628705 -0.9551309 -0.4288246 -0.9670793
#> Chem_100021502  1.7804216 -0.7856592 -1.1449344  0.3140216
#> Chem_100020519 -1.4069868  0.8426062 -1.1347284  0.7913141
#> Chem_100008903 -0.4368991  1.0047319  0.6997124 -1.0423101
#> Chem_100009217 -1.1949032 -1.2387747 -0.1938945 -0.8069717
#>                   P101052    P101054     P101056    P101057
#> Chem_100002945  1.0224999 -2.3147133 -1.47963462  0.3075335
#> Chem_100002356  1.0966509 -1.7937897 -1.92297496  0.1747746
#> Chem_100021502 -0.5945462 -1.0100058 -0.09086844 -0.4991725
#> Chem_100020519 -1.4558152 -1.5854144  0.80449715 -1.2716007
#> Chem_100008903 -1.5016496  0.5741487 -1.13225320 -0.2245752
#> Chem_100009217 -0.1255028  0.7811907 -0.40384813  0.2088243
#>                   P101059      P101061     P101062
#> Chem_100002945 -0.9041378 -0.830961764 -0.04269006
#> Chem_100002356 -1.1679764 -0.272491428  0.52868460
#> Chem_100021502  0.3712173 -1.093359526  1.04814884
#> Chem_100020519  0.7948727 -1.280571027 -0.71390394
#> Chem_100008903  0.6100619  0.009185413  0.48632983
#> Chem_100009217  0.9922981 -0.369602391  0.64507766
#>                    P101064     P101065      P101067
#> Chem_100002945 -0.29816470 -0.96715774  0.005732514
#> Chem_100002356 -0.64564813 -0.96580922 -0.055647831
#> Chem_100021502  1.17546656  0.05328984 -1.428569874
#> Chem_100020519  0.82422940 -0.84584830  0.780033076
#> Chem_100008903  0.02282356 -0.41934317 -0.078237408
#> Chem_100009217  0.63353253 -0.33768720 -1.057981752
#>                    P101068    P101069    P101071    P101072
#> Chem_100002945  0.60689290  1.0368525 -0.2442709 -0.1072135
#> Chem_100002356  1.32706991  0.1263946 -0.6112047  0.4674189
#> Chem_100021502  0.08575028  0.1369518 -0.1513205 -1.3147063
#> Chem_100020519 -1.35963789  0.8257427  0.7869318 -0.5680769
#> Chem_100008903 -1.15371004 -1.0172760  0.3253629  0.6506648
#> Chem_100009217 -1.51159030 -1.9068013  0.2212904  1.5332536
#>                   P101074    P101075    P101076    P101077
#> Chem_100002945 -0.5653190 -0.6045877 -0.4761183 -1.0070099
#> Chem_100002356 -0.5306034 -0.4645580  0.1769096 -0.5154926
#> Chem_100021502  0.3364421 -1.8590448 -1.4555677 -0.3873513
#> Chem_100020519 -0.4537732  0.8220845  0.8044084  0.7469109
#> Chem_100008903  0.9489675  0.1351589 -0.3237120  1.1571875
#> Chem_100009217 -0.2530297 -0.4551003 -0.2640308  1.9586874
#>                    P101079    P101080    P101081    P101082
#> Chem_100002945 -1.42377523 -0.3557742  1.2846662 -0.2028952
#> Chem_100002356 -0.55867549 -0.7410245  1.0386309 -0.7075808
#> Chem_100021502 -0.07531674  1.3138423 -1.3179406  0.3723564
#> Chem_100020519  0.83032979 -1.7057653  0.7875397  0.8943481
#> Chem_100008903 -0.02373964  1.6753779  0.2373957  1.1565220
#> Chem_100009217  1.03322697  1.3570251  0.4346865 -0.0908827
#>                   P101084    P101085   P101088    P101090
#> Chem_100002945 -0.5571438 -0.6432060 0.2680558 -0.8631570
#> Chem_100002356 -0.1780117  0.2362866 0.8376409 -1.0465957
#> Chem_100021502 -0.5129976 -0.9397795 0.3032375  0.4512293
#> Chem_100020519  0.8016100 -1.5150492 0.8345681  0.8427964
#> Chem_100008903 -0.6930137  1.0481281 0.9958996  0.2774736
#> Chem_100009217  0.1433649  0.2026186 0.9426616  0.6029496
#>                    P101094    P101095    P101096
#> Chem_100002945  0.60534398 -0.1098016  1.4375694
#> Chem_100002356  0.95059563 -0.3602152  1.1136725
#> Chem_100021502 -0.53483695  3.0939556  0.8708722
#> Chem_100020519 -1.32067449  0.8750618 -1.1626365
#> Chem_100008903  0.05615684 -2.3429673 -1.0188086
#> Chem_100009217 -0.63260805 -0.6079294 -1.6079520
```


## 数据分布变化

+ 函数来自于**POMA**包

```r
POMABoxplots <- function(
    data,
    group,
    feature_type = "samples",    
    jitter = FALSE,
    feature_name = NULL,
    show_number = NULL,
    label_size = 10,
    legend_position = "bottom") {
  
  # data = se_impute
  # group = "LiverFatClass"
  # feature_type = "samples"
  # jitter = FALSE
  # feature_name = NULL
  # show_number = 10
  # label_size = 10
  # legend_position = "bottom"
  
  if (missing(data)) {
    stop("data argument is empty!")
  }
  
  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee SummarizedExperiment::SummarizedExperiment")
  }
  
  if (!(feature_type %in% c("samples", "features"))) {
    stop("Incorrect value for group argument!")
  }
  
  if (!is.null(feature_name)) {
    if(!any(feature_name %in% rownames(SummarizedExperiment::assay(data)))) {
      stop("None of the specified features found")
    }
    if(!all(feature_name %in% rownames(SummarizedExperiment::assay(data)))){
      warning(paste0("Feature/s ",
                     paste0(feature_name[!feature_name %in% rownames(SummarizedExperiment::assay(data))], collapse = ", "),
                     " not found"))
    }
  }
  
  
  if(!(legend_position %in% c("none", "top", "bottom", "left", "right"))) {
    stop("Incorrect value for legend_position argument!")
  }
  
  e <- t(SummarizedExperiment::assay(data))
  target <- SummarizedExperiment::colData(data) %>%
    as.data.frame() %>% 
    tibble::rownames_to_column("ID") %>%
    dplyr::select(c("ID", group))
  colnames(target)[which(colnames(target) == group)] <- "TempGroup"
  data <- cbind(target, e) %>%
    dplyr::arrange(desc(ID))

  if (feature_type == "samples") {
    if (is.null(show_number)) {
      plot_data <- data %>%
        tidyr::pivot_longer(cols = -c(ID, TempGroup)) %>%
        ggplot2::ggplot(ggplot2::aes(ID, value, color = TempGroup))      
    } else {
      selected_samples <- target$ID[1:show_number]
      plot_data <- data %>%
        tidyr::pivot_longer(cols = -c(ID, TempGroup)) %>%
        dplyr::filter(ID %in% selected_samples) %>%
        ggplot2::ggplot(ggplot2::aes(ID, value, color = TempGroup))       
    }
  } else {
    if(all(is.null(feature_name), is.null(show_number))) {
      plot_data <- data %>%
        dplyr::select(-ID) %>%
        tidyr::pivot_longer(cols = -TempGroup) %>%
        ggplot2::ggplot(ggplot2::aes(name, value, color = TempGroup))
      
    } else if (all(!is.null(feature_name), is.null(show_number))) {
      plot_data <- data %>%
        dplyr::select(-ID) %>%
        tidyr::pivot_longer(cols = -TempGroup) %>%
        dplyr::filter(name %in% feature_name) %>%
        ggplot2::ggplot(ggplot2::aes(name, value, color = TempGroup))
    } else if (all(is.null(feature_name), !is.null(show_number))) {
      selected_features <- colnames(e)[1:show_number]
      plot_data <- data %>%
        dplyr::select(-ID) %>%
        tidyr::pivot_longer(cols = -TempGroup) %>%
        dplyr::filter(name %in% selected_features) %>%
        ggplot2::ggplot(ggplot2::aes(name, value, color = TempGroup))
    } else if (all(!is.null(feature_name), !is.null(show_number))) {
      plot_data <- data %>%
        dplyr::select(-ID) %>%
        tidyr::pivot_longer(cols = -TempGroup) %>%
        dplyr::filter(name %in% feature_name) %>%
        ggplot2::ggplot(ggplot2::aes(name, value, color = TempGroup))
    }
  }
  
  plot_complete <- plot_data +
    ggplot2::geom_boxplot() +
    {if(jitter)ggplot2::geom_jitter(alpha = 0.5, position = ggplot2::position_jitterdodge())} +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "", 
                  y = "Value") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = label_size),
                   legend.title = ggplot2::element_blank(),
                   legend.position = legend_position) +
    ggplot2::scale_colour_viridis_d(begin = 0, end = 0.8)
  
  return(plot_complete)
}


POMADensity <- function(
    data,
    group,
    feature_type = "features",    
    feature_name = NULL,
    show_number = NULL,
    legend_position = "bottom") {
  
  # data = se_impute
  # group = "LiverFatClass"
  # feature_type = "features"
  # feature_name = NULL
  # show_number = 10
  # legend_position = "bottom"

  if (missing(data)) {
    stop("data argument is empty!")
  }
  
  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee SummarizedExperiment::SummarizedExperiment")
  }
  
  if (!(feature_type %in% c("samples", "features"))) {
    stop("Incorrect value for group argument!")
  }
  
  if (!is.null(feature_name)) {
    if(!any(feature_name %in% rownames(SummarizedExperiment::assay(data)))) {
      stop("None of the specified features found")
    }
    if(!all(feature_name %in% rownames(SummarizedExperiment::assay(data)))){
      warning(paste0("Feature/s ",
                     paste0(feature_name[!feature_name %in% rownames(SummarizedExperiment::assay(data))], collapse = ", "),
                     " not found"))
    }
  }
  
  if(!(legend_position %in% c("none", "top", "bottom", "left", "right"))) {
    stop("Incorrect value for legend_position argument!")
  }
  
  e <- t(SummarizedExperiment::assay(data))
  target <- SummarizedExperiment::colData(data) %>%
    as.data.frame() %>% 
    tibble::rownames_to_column("ID") %>%
    dplyr::select(c("ID", group))
  colnames(target)[which(colnames(target) == group)] <- "TempGroup"
  data <- cbind(target, e)
  
  if (feature_type == "samples") {
    if (is.null(show_number)) {
      plot_data <- data %>%
        tidyr::pivot_longer(cols = -c(ID, TempGroup)) %>%
        ggplot2::ggplot(ggplot2::aes(value, fill = TempGroup))      
    } else {
      selected_samples <- target$ID[1:show_number]      
      plot_data <- data %>%
        tidyr::pivot_longer(cols = -c(ID, TempGroup)) %>%
        dplyr::filter(ID %in% selected_samples) %>%
        ggplot2::ggplot(ggplot2::aes(value, fill = TempGroup))       
    }
  } else {
    if(all(is.null(feature_name), is.null(show_number))) {
      plot_data <- data %>%
        dplyr::select(-ID) %>%
        tidyr::pivot_longer(cols = -TempGroup) %>%
        ggplot2::ggplot(ggplot2::aes(value, fill = TempGroup))
      
    } else if (all(!is.null(feature_name), is.null(show_number))) {
      plot_data <- data %>%
        dplyr::select(-ID) %>%
        tidyr::pivot_longer(cols = -TempGroup) %>%
        dplyr::filter(name %in% feature_name) %>%
        ggplot2::ggplot(ggplot2::aes(value, fill = TempGroup))
    } else if (all(is.null(feature_name), !is.null(show_number))) {
      selected_features <- colnames(e)[1:show_number]      
      plot_data <- data %>%
        dplyr::select(-ID) %>%
        tidyr::pivot_longer(cols = -TempGroup) %>%
        dplyr::filter(name %in% selected_features) %>%
        ggplot2::ggplot(ggplot2::aes(value, fill = TempGroup))
    } else if (all(!is.null(feature_name), !is.null(show_number))) {
      plot_data <- data %>%
        dplyr::select(-ID) %>%
        tidyr::pivot_longer(cols = -TempGroup) %>%
        dplyr::filter(name %in% feature_name) %>%
        ggplot2::ggplot(ggplot2::aes(value, fill = TempGroup))
    }
  }  
  
  plot_complete <- plot_data +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Value",
                  y = "Density") +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = legend_position) +
    ggplot2::scale_fill_viridis_d(begin = 0, end = 0.8)
  
  return(plot_complete)
  
}


get_distribution <- function(
    datset,
    Type = c("raw", "check", "filter", 
             "impute", "norm_relative", "norm_log10",
             "norm_scale")) {
  
  # datset = se
  # Type = "raw"
  
  dat <- SummarizedExperiment::assay(datset) %>% 
      data.frame() %>%
      rownames_to_column("name") %>%
      tidyr::gather(key = "sample", value = "value", -name)    
  
  if (Type == "raw") {
    pl <- ggplot(dat, aes(x = value)) + 
            geom_histogram(color = "black", fill = "white") +
            labs(title = "Distribution of Raw \n Metabolites Intensity", 
                 x = "Raw Intensity", y = "Frequency") +
            theme_bw()    
  } else if (Type == "check") {
    pl <- ggplot(dat, aes(x = value)) + 
            geom_histogram(color = "black", fill = "white") +
            labs(title = "Distribution of check \n Metabolites Intensity", 
                 x = "checked Intensity", y = "Frequency") +
            theme_bw()    
  } else if (Type == "filter") {
    pl <- ggplot(dat, aes(x = value)) + 
            geom_histogram(color = "black", fill = "white")+
            labs(title = "Distribution of filter\n Metabolites Intensity", 
                 x = "filtered Intensity", y = "Frequency")+
            theme_bw()    
  } else if (Type == "impute") {
    pl <- ggplot(dat, aes(x = value)) + 
            geom_histogram(color = "black", fill = "white")+
            labs(title = "Distribution of impute\n Metabolites Intensity", 
                 x = "imputed Intensity", y = "Frequency")+
            theme_bw()    
  } else if (Type == "norm_relative") {
    pl <- ggplot(dat, aes(x = value)) + 
            geom_histogram(color = "black", fill = "white")+
            labs(title = "Distribution of norm\n Metabolites Intensity", 
                 x = "norm \n relative abundance", y="Frequency")+
            theme_bw()    
  } else if (Type == "norm_log10") {
    pl <- ggplot(dat, aes(x = value)) + 
            geom_histogram(color = "black", fill = "white")+
            labs(title = "Distribution of norm\n Metabolites Intensity", 
                 x = "norm \n log10(Intensity)", y = "Frequency")+
            theme_bw()    
  } else if (Type == "norm_scale") {
    pl <- ggplot(dat, aes(x = value)) + 
            geom_histogram(color = "black", fill = "white")+
            labs(title = "Distribution of norm\n Metabolites Intensity", 
                 x = "norm \n Scale(Intensity)", y = "Frequency")+
            theme_bw()    
  }
  
  return(pl)
}
```

+ 以样本为基点的代谢物数据分布情况（组间箱线图）

```r
pl_unnor <- POMABoxplots(data = se_impute, group = "LiverFatClass", 
                         feature_type = "samples", jitter = FALSE, show_number = 10) +
  ggtitle("Not Normalized (imputation)") +
  theme(legend.position = "none") 

pl_nor_log <- POMABoxplots(data = se_tran, group = "LiverFatClass", 
                          feature_type = "samples", jitter = FALSE, show_number = 10) +
  ggtitle("Normalized (log10)")

pl_nor_rb <- POMABoxplots(data = se_norm, group = "LiverFatClass", 
                          feature_type = "samples", jitter = FALSE, show_number = 10) +
  ggtitle("Normalized (relative abundance)")

pl_nor_zscore <- POMABoxplots(data = se_scale, group = "LiverFatClass", 
                          feature_type = "samples", jitter = FALSE, show_number = 10) +
  ggtitle("Normalized (Zscore)")

cowplot::plot_grid(pl_unnor, pl_nor_log, 
                   pl_nor_rb, pl_nor_zscore, 
                   ncol = 1, align = "v",
                   labels = LETTERS[1:4])
```

<img src="20-Metabolome_DataProcess_files/figure-html/unnamed-chunk-15-1.png" width="100%" />

+ 以代谢物为基点的代谢物数据分布情况（组间箱线图）

```r
pl_unnor <- POMADensity(data = se_impute, group = "LiverFatClass", feature_type = "features") +
  ggtitle("Not Normalized") +
  theme(legend.position = "none") # data before normalization

pl_nor_log <- POMADensity(data = se_tran, group = "LiverFatClass", feature_type = "features") +
  ggtitle("Normalized (log10)") # data after normalization

pl_nor_rb <- POMADensity(data = se_norm, group = "LiverFatClass", feature_type = "features") +
  ggtitle("Normalized (relative abundance)") # data after normalization

pl_nor_zscore <- POMADensity(data = se_scale, group = "LiverFatClass", feature_type = "features") +
  ggtitle("Normalized (Zscore)") # data after normalization

cowplot::plot_grid(pl_unnor, pl_nor_log, 
                   pl_nor_rb, pl_nor_zscore, 
                   ncol = 1, align = "v",
                   labels = LETTERS[1:4])
```

<img src="20-Metabolome_DataProcess_files/figure-html/unnamed-chunk-16-1.png" width="100%" />

+ 所有以代谢物为基点的数据分布


```r
raw_pl <- get_distribution(datset = data_meta_new, Type = "raw")
check_pl <- get_distribution(datset = se_check, Type = "check")
filter_pl <- get_distribution(datset = se_filter, Type = "filter")
impute_pl <- get_distribution(datset = se_impute, Type = "impute")
norm_log10_pl <- get_distribution(datset = se_tran, Type = "norm_log10")
norm_relative_pl <- get_distribution(datset = se_norm, Type = "norm_relative")
norm_scale_pl <- get_distribution(datset = se_scale, Type = "norm_scale")

cowplot::plot_grid(raw_pl, check_pl, 
                   filter_pl, impute_pl, 
                   norm_log10_pl, norm_relative_pl,
                   norm_scale_pl,
                   align = "hv", nrow = 2,
                   labels = LETTERS[1:7])
```

<img src="20-Metabolome_DataProcess_files/figure-html/unnamed-chunk-17-1.png" width="100%" />


结果：

+ 补缺失值没有改变数据分布状态，仍然是偏态分布

+ log10单个数据转换将偏态分布转成偏向正态分布，而在该基础上的相对丰度则由偏向了正态分布

+ relative abundance在样本内部归一化其所有的特征，数据分布没有发生任何变化

+ Zscore是跨样本针对特征归一化其，数据分布呈现标准正态分布


## 保存数据

```r
if (!dir.exists("./InputData/result/QC")) {
  dir.create("./InputData/result/QC", recursive = TRUE)
}

saveRDS(data_meta_new, "./InputData/result/QC/se_raw.RDS", compress = TRUE)
saveRDS(se_check, "./InputData/result/QC/se_check.RDS", compress = TRUE)
saveRDS(se_impute, "./InputData/result/QC/se_impute.RDS", compress = TRUE)
saveRDS(se_filter, "./InputData/result/QC/se_filter.RDS", compress = TRUE)

saveRDS(se_tran, "./InputData/result/QC/se_tran.RDS", compress = TRUE)
saveRDS(se_norm, "./InputData/result/QC/se_norm.RDS", compress = TRUE)
saveRDS(se_scale, "./InputData/result/QC/se_scale.RDS", compress = TRUE)
```



## 总结

数据预处理是代谢组分析较为重要的步骤，通常建议使用log+zscore的方法对数据归一化处理，用于后续的统计分析，但看到有关文章在计算Log2FoldChange的时候使用原始intensity值计算。


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
#>  date     2023-12-08
#>  pandoc   3.1.3 @ /Users/zouhua/opt/anaconda3/bin/ (via rmarkdown)
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
#>  cowplot                1.1.1     2020-12-30 [2] CRAN (R 4.1.0)
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
#>  farver                 2.1.1     2022-07-06 [2] CRAN (R 4.1.2)
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
#>  impute                 1.68.0    2021-10-26 [2] Bioconductor
#>  IRanges              * 2.28.0    2021-10-26 [2] Bioconductor
#>  iterators              1.0.14    2022-02-05 [2] CRAN (R 4.1.2)
#>  jquerylib              0.1.4     2021-04-26 [2] CRAN (R 4.1.0)
#>  jsonlite               1.8.7     2023-06-29 [2] CRAN (R 4.1.3)
#>  KEGGREST               1.34.0    2021-10-26 [2] Bioconductor
#>  KernSmooth             2.23-22   2023-07-10 [2] CRAN (R 4.1.3)
#>  knitr                  1.43      2023-05-25 [2] CRAN (R 4.1.3)
#>  labeling               0.4.2     2020-10-20 [2] CRAN (R 4.1.0)
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
#>  viridisLite            0.4.2     2023-05-02 [2] CRAN (R 4.1.2)
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

