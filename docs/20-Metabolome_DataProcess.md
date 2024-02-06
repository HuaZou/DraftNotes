


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
#> [1] "A total of 7873 (12.5%) missing values were detected."
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
#> dim: 955 61 
#> metadata(0):
#> assays(1): ''
#> rownames(955): Chem_100002945 Chem_100002356 ...
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
#> [1] "Removed  100  features based on QC RSD values. QC samples are excluded from downstream functional analysis."

se_filter
#> class: SummarizedExperiment 
#> dim: 641 55 
#> metadata(0):
#> assays(1): ''
#> rownames(641): Chem_100002945 Chem_100002356 ...
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
#> dim: 641 55 
#> metadata(0):
#> assays(1): ''
#> rownames(641): Chem_100002945 Chem_100002356 ...
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
#> dim: 641 55 
#> metadata(0):
#> assays(1): ''
#> rownames(641): Chem_100002945 Chem_100002356 ...
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
#> dim: 641 55 
#> metadata(0):
#> assays(1): ''
#> rownames(641): Chem_100002945 Chem_100002356 ...
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
#>                   P101001      P101003      P101004
#> Chem_100002945 51127588.0  42040432.00  34940596.00
#> Chem_100002356  5105020.5   4006120.25   3885477.00
#> Chem_100008903 94392176.0 117463144.00 115155104.00
#> Chem_100009217   144048.2    246049.92    217577.44
#> Chem_100000657 25632184.0  26952236.00  25106562.00
#> Chem_100001397   123026.4     86646.23     22810.19
#>                   P101007     P101009    P101010    P101011
#> Chem_100002945 58518636.0  51118832.0 83783688.0 29017984.0
#> Chem_100002356  4285129.5   6665653.5  9057441.0  2802655.2
#> Chem_100008903 79582632.0 118408760.0 92508664.0 94076424.0
#> Chem_100009217   211262.5    431295.3   135226.5   208435.8
#> Chem_100000657 31371314.0  27787270.0 26685844.0 27780988.0
#> Chem_100001397   375686.5    118662.2   130040.2   258790.7
#>                   P101012    P101013     P101016
#> Chem_100002945 51222064.0 77550128.0 30949554.00
#> Chem_100002356  5996555.0 11367511.0  3874736.75
#> Chem_100008903 69473744.0 80567352.0 98766592.00
#> Chem_100009217   163974.4   162897.1   276264.59
#> Chem_100000657 30158644.0 27591838.0 21537830.00
#> Chem_100001397   121337.8   526472.9    51847.78
#>                    P101017     P101018    P101019
#> Chem_100002945  26923596.0  56720032.0 27956064.0
#> Chem_100002356   2817151.0   8029728.0  3766663.8
#> Chem_100008903 129547656.0 118271584.0 37880820.0
#> Chem_100009217    375196.6    538405.4   724898.5
#> Chem_100000657  34510160.0  24383964.0 25435894.0
#> Chem_100001397    309598.5    142944.8   122802.0
#>                     P101021     P101022     P101024
#> Chem_100002945  48723600.00  16282054.0 77028824.00
#> Chem_100002356   5174967.00   1746182.9  5519105.50
#> Chem_100008903 163868720.00 106189520.0 71475920.00
#> Chem_100009217    242844.06    447587.9   116601.31
#> Chem_100000657  33230896.00  28310448.0 24481250.00
#> Chem_100001397     36066.18    143211.6    41351.85
#>                    P101025     P101027     P101030
#> Chem_100002945  32022342.0  22589448.0  38449788.0
#> Chem_100002356   2557365.0   1882902.6   2860324.2
#> Chem_100008903 147776592.0 127571792.0 128319128.0
#> Chem_100009217   1936915.5    322639.6    381134.2
#> Chem_100000657  36634552.0  26803284.0  38250240.0
#> Chem_100001397    299205.7    338827.8    240744.7
#>                   P101031    P101038      P101041
#> Chem_100002945 59134052.0 32038030.0  20833830.00
#> Chem_100002356  4721201.0  4011627.5   2938779.00
#> Chem_100008903 90447848.0 46622592.0 111919096.00
#> Chem_100009217   171138.8   203968.3    288088.38
#> Chem_100000657 31415658.0 36232668.0  27869314.00
#> Chem_100001397   748813.9   497035.7     45508.17
#>                   P101042    P101047     P101050
#> Chem_100002945 33809080.0 18637508.0  21978476.0
#> Chem_100002356  3017260.5  1935144.1   2897211.0
#> Chem_100008903 89762056.0 97617984.0 112900000.0
#> Chem_100009217   143925.2   112179.6    222175.4
#> Chem_100000657 29975144.0 30274784.0  29706444.0
#> Chem_100001397   141089.5   149316.7    167735.3
#>                    P101051    P101052     P101054
#> Chem_100002945 24265162.00 52203780.0  12836384.0
#> Chem_100002356  2476279.50  5928454.0   1685760.6
#> Chem_100008903 71779912.00 64008512.0 111279888.0
#> Chem_100009217   171872.37   252141.2    371306.8
#> Chem_100000657 28526700.00 19797506.0  22277808.0
#> Chem_100001397    83878.49   153736.2    131393.6
#>                   P101056     P101057      P101059
#> Chem_100002945 18546636.0 32301820.00  22645984.00
#> Chem_100002356  1650011.0  3419157.75   2196044.25
#> Chem_100008903 68771592.0 77140264.00 113564872.00
#> Chem_100009217   211104.1   254912.39    416259.56
#> Chem_100000657 20182492.0 22884156.00  28125854.00
#> Chem_100001397   125445.7    36627.55     67012.24
#>                   P101061    P101062      P101064
#> Chem_100002945 23683254.0 29027646.0  32629048.00
#> Chem_100002356  3217499.2  4060367.8   3031529.50
#> Chem_100008903 96143304.0 98940424.0 108473368.00
#> Chem_100009217   197874.2   325386.6    380257.62
#> Chem_100000657 27426428.0 27010240.0  29414698.00
#> Chem_100001397   248423.8    35064.2     32855.53
#>                   P101065    P101067    P101068     P101069
#> Chem_100002945 22950806.0 33555116.0 44283972.0 52685972.00
#> Chem_100002356  2467147.2  3567913.2  6525382.0  3984333.50
#> Chem_100008903 86418592.0 95236288.0 71289048.0 74526792.00
#> Chem_100009217   217008.5   149385.6   121522.4    99392.66
#> Chem_100000657 25252692.0 33943364.0 28292364.0 31375248.00
#> Chem_100001397   137605.2   155844.5   109021.4   178198.00
#>                     P101071     P101072     P101074
#> Chem_100002945  32415040.00  34170948.0  22550616.0
#> Chem_100002356   3001414.50   4679519.0   2529255.5
#> Chem_100008903 115519872.00 127401592.0 108255040.0
#> Chem_100009217    301355.72    590274.3    199429.7
#> Chem_100000657  33193432.00  27121940.0  24320710.0
#> Chem_100001397     74016.23    404225.1    124975.2
#>                   P101075    P101076     P101077    P101079
#> Chem_100002945 22058076.0 24455466.0  25225170.0 15718590.0
#> Chem_100002356  2583265.5  3515218.2   3272875.0  2449462.5
#> Chem_100008903 83989120.0 77315256.0 158257952.0 78587928.0
#> Chem_100009217   179174.6   204990.0    767403.1   375777.7
#> Chem_100000657 25176120.0 25405452.0  25782852.0 26353500.0
#> Chem_100001397   167129.7   137572.7    410891.4   292623.9
#>                    P101080     P101081     P101082
#> Chem_100002945  29120336.0  65904836.0  22908578.0
#> Chem_100002356   2695001.5   6474709.5   2110243.8
#> Chem_100008903 163246832.0 124678488.0 100435064.0
#> Chem_100009217    514687.6    360975.7    197618.8
#> Chem_100000657  25819136.0  34524436.0  27735140.0
#> Chem_100001397    162495.1    237136.2    154358.6
#>                   P101084     P101085      P101088
#> Chem_100002945 29140440.0  20427124.0  29199012.00
#> Chem_100002356  3648091.2   3253531.8   4154170.75
#> Chem_100008903 86139200.0 103513520.0 101921248.00
#> Chem_100009217   293438.9    238898.8    347081.22
#> Chem_100000657 30228504.0  29317998.0  28259446.00
#> Chem_100001397   181003.1    126785.3     36267.95
#>                    P101090     P101094     P101095
#> Chem_100002945  24042020.0 36910084.00 35662068.00
#> Chem_100002356   2396959.8  4759584.50  3452283.25
#> Chem_100008903 107571936.0 85426888.00 53107852.00
#> Chem_100009217    352560.8   168247.59   551062.83
#> Chem_100000657  32471818.0 25804370.00 25684144.00
#> Chem_100001397    217137.2    75320.37    59115.54
#>                    P101096
#> Chem_100002945 66402192.00
#> Chem_100002356  6374383.00
#> Chem_100008903 80095704.00
#> Chem_100009217   121388.24
#> Chem_100000657 28136850.00
#> Chem_100001397    81256.43
```


```r
SummarizedExperiment::assay(se_scale) %>% as.data.frame() %>% head()
#>                   P101001      P101003    P101004
#> Chem_100002945  1.0449213  0.831722997  0.5381731
#> Chem_100002356  0.8063931  0.449454631  0.5192783
#> Chem_100008903 -0.1198802  0.983934326  1.1581141
#> Chem_100009217 -1.1485476 -0.009281135 -0.1495838
#> Chem_100000657 -0.7460131  0.309588904  0.2793535
#> Chem_100001397 -0.1200006 -0.489903593 -2.2163154
#>                   P101007    P101009     P101010
#> Chem_100002945  1.0779000  1.2189442  2.48823486
#> Chem_100002356  0.1293967  1.6041200  2.38092871
#> Chem_100008903 -1.0899843  0.8686275  0.09896180
#> Chem_100009217 -0.5921313  0.9841933 -1.16233297
#> Chem_100000657 -0.1793197  0.2491561  0.06532847
#> Chem_100001397  1.2426876 -0.1067152  0.02705011
#>                    P101011    P101012    P101013    P101016
#> Chem_100002945 -0.20197614  0.9122919  1.9294958  0.2720216
#> Chem_100002356 -0.51157735  1.0792773  2.6160224  0.5457083
#> Chem_100008903  0.09889906 -1.3142293 -0.8431960  0.7020107
#> Chem_100009217 -0.37558543 -0.9771321 -0.9950725  0.3196353
#> Chem_100000657  0.23045073 -0.0143958 -0.6384568 -0.6569034
#> Chem_100001397  0.92753156 -0.1864659  1.7466355 -1.1048399
#>                   P101017     P101018     P101019
#> Chem_100002945 -0.7827304  0.98592911 -0.21249088
#> Chem_100002356 -0.8321839  1.62977030  0.27949585
#> Chem_100008903  0.5878083  0.18567326 -2.79746243
#> Chem_100009217  0.5002016  1.12599829  1.99752238
#> Chem_100000657  0.6141785 -1.89522369 -0.14173146
#> Chem_100001397  1.0110060 -0.03354461 -0.03339455
#>                   P101021    P101022    P101024    P101025
#> Chem_100002945  1.0284134 -1.5214384  2.2499589 -0.7782133
#> Chem_100002356  0.9273960 -1.5592975  1.1499730 -1.4191452
#> Chem_100008903  1.8431647  0.6818344 -0.7917393  0.4071834
#> Chem_100009217 -0.1241391  1.1222149 -1.4524581  3.2676296
#> Chem_100000657  1.2671270  0.7023805 -0.5865167 -0.1406785
#> Chem_100001397 -1.7124181  0.1880810 -1.5062022  0.8011251
#>                  P101027    P101030     P101031    P101038
#> Chem_100002945 -1.029400  0.2351081  1.13748744 -0.1712154
#> Chem_100002356 -1.648225 -0.6826398  0.39116366  0.1760698
#> Chem_100008903  0.811321  0.7489955 -0.62560628 -2.5150752
#> Chem_100009217  0.326320  0.6039177 -0.96384269 -0.5303699
#> Chem_100000657 -0.560201  1.6665523 -0.08223301  1.4304311
#> Chem_100001397  1.204952  0.7310413  2.16348025  1.7076679
#>                   P101041      P101042    P101047
#> Chem_100002945 -1.0289888 -0.140250512 -0.7406186
#> Chem_100002356 -0.3955722 -0.600380116 -0.9353836
#> Chem_100008903  0.6739563 -0.503717132  1.0461782
#> Chem_100009217  0.2283293 -1.229952162 -1.2368496
#> Chem_100000657  0.2534160 -0.114099301  2.3816195
#> Chem_100001397 -1.3828830  0.004630317  0.4112869
#>                   P101050     P101051    P101052
#> Chem_100002945 -0.9151959 -0.76594474  0.9820107
#> Chem_100002356 -0.4472592 -0.91202351  1.0710136
#> Chem_100008903  0.6740586 -0.96048979 -1.5529674
#> Chem_100009217 -0.2669785 -0.79733879 -0.1683491
#> Chem_100000657  0.6287515  0.08996049 -2.7700821
#> Chem_100001397  0.3437965 -0.61178487  0.1340711
#>                    P101054    P101056    P101057    P101059
#> Chem_100002945 -2.36847759 -1.4082840  0.3006656 -0.8890482
#> Chem_100002356 -1.85999355 -1.8754002  0.1718931 -1.1625818
#> Chem_100008903  0.46048497 -1.0629946 -0.2308007  0.6230809
#> Chem_100009217  0.62565059 -0.4001338  0.1260828  0.8757329
#> Chem_100000657 -1.61532062 -2.1590584 -0.4543482  0.1284644
#> Chem_100001397 -0.02308227 -0.0680733 -1.5957044 -0.8923146
#>                   P101061     P101062     P101064
#> Chem_100002945 -0.9028297 -0.09052683 -0.26875362
#> Chem_100002356 -0.3418457  0.48673474 -0.62297503
#> Chem_100008903 -0.1064402  0.42244966  0.06046142
#> Chem_100009217 -0.5763039  0.51607422  0.54621720
#> Chem_100000657 -0.3817293  0.33790109 -0.34845091
#> Chem_100001397  0.7978595 -1.69380781 -1.93124220
#>                    P101065     P101067    P101068
#> Chem_100002945 -0.95950517  0.01604351  0.5850408
#> Chem_100002356 -0.96743778 -0.04660294  1.3130450
#> Chem_100008903 -0.42681751 -0.06385235 -1.1835973
#> Chem_100009217 -0.39333584 -1.07118105 -1.5152935
#> Chem_100000657 -0.87685244  1.18825654 -0.3526311
#> Chem_100001397  0.02408398  0.20013929 -0.3159273
#>                   P101069    P101071    P101072     P101074
#> Chem_100002945  1.1824673 -0.2007193 -0.1980506 -0.57508523
#> Chem_100002356  0.2574325 -0.5761728  0.3877867 -0.54359676
#> Chem_100008903 -0.8134078  0.3866606  0.5256379  0.93434237
#> Chem_100009217 -1.8076170  0.1620392  1.3345367 -0.31981063
#> Chem_100000657  0.7738432  0.6857132 -1.0061525  0.02996046
#> Chem_100001397  0.3934027 -0.8314649  1.3602053  0.05458039
#>                   P101075    P101076    P101077     P101079
#> Chem_100002945 -0.6183263 -0.4228775 -0.9825753 -1.41183410
#> Chem_100002356 -0.4813953  0.2183863 -0.5035003 -0.56136716
#> Chem_100008903  0.1075950 -0.2581992  1.1869791 -0.02905182
#> Chem_100009217 -0.5142242 -0.2958941  1.7995391  0.90953082
#> Chem_100000657  0.2972914  0.1861575 -1.4339898  0.76855287
#> Chem_100001397  0.4479513  0.1639357  1.3682446  1.22081846
#>                   P101080   P101081    P101082     P101084
#> Chem_100002945 -0.3190578 1.2446150 -0.1765249 -0.62423078
#> Chem_100002356 -0.7127031 1.0160673 -0.6866051 -0.24047674
#> Chem_100008903  1.7378650 0.2045521  1.2015743 -0.80152840
#> Chem_100009217  1.2383436 0.3289937 -0.1440551  0.02563072
#> Chem_100000657 -0.5992344 0.1284076  1.8976686 -0.36988796
#> Chem_100001397  0.2621112 0.5910000  0.4716894  0.29003762
#>                    P101085    P101088    P101090    P101094
#> Chem_100002945 -0.71844665  0.1882081 -0.8018269  0.4592286
#> Chem_100002356  0.16473295  0.7707010 -1.0019370  0.8300756
#> Chem_100008903  0.93798723  0.8957131  0.3548741 -0.1362861
#> Chem_100009217  0.07678649  0.7830551  0.5317072 -0.7439557
#> Chem_100000657  1.58471748  1.3529428  0.9240312 -0.1070827
#> Chem_100001397  0.11288325 -1.5630066  0.6442162 -0.6925749
#>                    P101095    P101096
#> Chem_100002945 -0.04499693  1.3615580
#> Chem_100002356 -0.30605646  1.0620177
#> Chem_100008903 -2.27633093 -1.1073245
#> Chem_100009217  1.23442053 -1.6285536
#> Chem_100000657 -1.24409768 -0.9775123
#> Chem_100001397 -1.15569756 -0.7808386
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
#>  Biobase                  * 2.62.0     2023-10-24 [1] Bioconductor
#>  BiocGenerics             * 0.48.1     2023-11-01 [1] Bioconductor
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
#>  boot                       1.3-28.1   2022-11-22 [1] CRAN (R 4.3.1)
#>  bslib                      0.6.1      2023-11-28 [1] CRAN (R 4.3.0)
#>  cachem                     1.0.8      2023-05-01 [1] CRAN (R 4.3.0)
#>  caTools                    1.18.2     2021-03-28 [1] CRAN (R 4.3.0)
#>  cellranger                 1.1.0      2016-07-27 [1] CRAN (R 4.3.0)
#>  checkmate                  2.3.1      2023-12-04 [1] CRAN (R 4.3.0)
#>  class                      7.3-22     2023-05-03 [1] CRAN (R 4.3.1)
#>  cli                        3.6.2      2023-12-11 [1] CRAN (R 4.3.0)
#>  cluster                    2.1.4      2022-08-22 [1] CRAN (R 4.3.1)
#>  codetools                  0.2-19     2023-02-01 [1] CRAN (R 4.3.1)
#>  colorspace                 2.1-0      2023-01-23 [1] CRAN (R 4.3.0)
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
#>  generics                   0.1.3      2022-07-05 [1] CRAN (R 4.3.0)
#>  GenomeInfoDb             * 1.38.5     2023-12-28 [1] Bioconductor 3.18 (R 4.3.2)
#>  GenomeInfoDbData           1.2.11     2024-01-24 [1] Bioconductor
#>  GenomicRanges            * 1.54.1     2023-10-29 [1] Bioconductor
#>  ggbeeswarm                 0.7.2      2023-04-29 [1] CRAN (R 4.3.0)
#>  ggplot2                  * 3.4.4      2023-10-12 [1] CRAN (R 4.3.0)
#>  ggrepel                    0.9.5      2024-01-10 [1] CRAN (R 4.3.0)
#>  gld                        2.6.6      2022-10-23 [1] CRAN (R 4.3.0)
#>  glmnet                     4.1-8      2023-08-22 [1] CRAN (R 4.3.0)
#>  glue                       1.7.0      2024-01-09 [1] CRAN (R 4.3.0)
#>  gmp                        0.7-4      2024-01-15 [1] CRAN (R 4.3.0)
#>  gplots                     3.1.3.1    2024-02-02 [1] CRAN (R 4.3.2)
#>  gridExtra                  2.3        2017-09-09 [1] CRAN (R 4.3.0)
#>  gsl                        2.1-8      2023-01-24 [1] CRAN (R 4.3.0)
#>  gtable                     0.3.4      2023-08-21 [1] CRAN (R 4.3.0)
#>  gtools                     3.9.5      2023-11-20 [1] CRAN (R 4.3.0)
#>  highr                      0.10       2022-12-22 [1] CRAN (R 4.3.0)
#>  Hmisc                      5.1-1      2023-09-12 [1] CRAN (R 4.3.0)
#>  hms                        1.1.3      2023-03-21 [1] CRAN (R 4.3.0)
#>  htmlTable                  2.4.2      2023-10-29 [1] CRAN (R 4.3.0)
#>  htmltools                  0.5.7      2023-11-03 [1] CRAN (R 4.3.0)
#>  htmlwidgets                1.6.4      2023-12-06 [1] CRAN (R 4.3.0)
#>  httpuv                     1.6.14     2024-01-26 [1] CRAN (R 4.3.2)
#>  httr                       1.4.7      2023-08-15 [1] CRAN (R 4.3.0)
#>  igraph                     2.0.1.1    2024-01-30 [1] CRAN (R 4.3.2)
#>  impute                     1.76.0     2023-10-24 [1] Bioconductor
#>  IRanges                  * 2.36.0     2023-10-24 [1] Bioconductor
#>  irlba                      2.3.5.1    2022-10-03 [1] CRAN (R 4.3.0)
#>  iterators                  1.0.14     2022-02-05 [1] CRAN (R 4.3.0)
#>  jquerylib                  0.1.4      2021-04-26 [1] CRAN (R 4.3.0)
#>  jsonlite                   1.8.8      2023-12-04 [1] CRAN (R 4.3.0)
#>  KernSmooth                 2.23-21    2023-05-03 [1] CRAN (R 4.3.1)
#>  knitr                      1.45       2023-10-30 [1] CRAN (R 4.3.0)
#>  labeling                   0.4.3      2023-08-29 [1] CRAN (R 4.3.0)
#>  later                      1.3.2      2023-12-06 [1] CRAN (R 4.3.0)
#>  lattice                    0.21-8     2023-04-05 [1] CRAN (R 4.3.1)
#>  lazyeval                   0.2.2      2019-03-15 [1] CRAN (R 4.3.0)
#>  lifecycle                  1.0.4      2023-11-07 [1] CRAN (R 4.3.0)
#>  limma                      3.58.1     2023-10-31 [1] Bioconductor
#>  lme4                       1.1-35.1   2023-11-05 [1] CRAN (R 4.3.0)
#>  lmerTest                   3.1-3      2020-10-23 [1] CRAN (R 4.3.0)
#>  lmom                       3.0        2023-08-29 [1] CRAN (R 4.3.0)
#>  locfit                     1.5-9.8    2023-06-11 [1] CRAN (R 4.3.0)
#>  lubridate                * 1.9.3      2023-09-27 [1] CRAN (R 4.3.0)
#>  magrittr                   2.0.3      2022-03-30 [1] CRAN (R 4.3.0)
#>  MASS                       7.3-60     2023-05-04 [1] CRAN (R 4.3.1)
#>  Matrix                     1.6-5      2024-01-11 [1] CRAN (R 4.3.0)
#>  MatrixGenerics           * 1.14.0     2023-10-24 [1] Bioconductor
#>  matrixStats              * 1.2.0      2023-12-11 [1] CRAN (R 4.3.0)
#>  memoise                    2.0.1      2021-11-26 [1] CRAN (R 4.3.0)
#>  metagenomeSeq              1.43.0     2023-04-25 [1] Bioconductor
#>  mgcv                       1.8-42     2023-03-02 [1] CRAN (R 4.3.1)
#>  mia                        1.10.0     2023-10-24 [1] Bioconductor
#>  MicrobiomeAnalysis       * 1.0.3      2024-02-06 [1] Github (HuaZou/MicrobiomeAnalysis@fd2a6a2)
#>  mime                       0.12       2021-09-28 [1] CRAN (R 4.3.0)
#>  miniUI                     0.1.1.1    2018-05-18 [1] CRAN (R 4.3.0)
#>  minqa                      1.2.6      2023-09-11 [1] CRAN (R 4.3.0)
#>  multcomp                   1.4-25     2023-06-20 [1] CRAN (R 4.3.0)
#>  MultiAssayExperiment       1.28.0     2023-10-24 [1] Bioconductor
#>  multtest                   2.58.0     2023-10-24 [1] Bioconductor
#>  munsell                    0.5.0      2018-06-12 [1] CRAN (R 4.3.0)
#>  mvtnorm                    1.2-4      2023-11-27 [1] CRAN (R 4.3.0)
#>  nlme                       3.1-162    2023-01-31 [1] CRAN (R 4.3.1)
#>  nloptr                     2.0.3      2022-05-26 [1] CRAN (R 4.3.0)
#>  nnet                       7.3-19     2023-05-03 [1] CRAN (R 4.3.1)
#>  numDeriv                   2016.8-1.1 2019-06-06 [1] CRAN (R 4.3.0)
#>  permute                    0.9-7      2022-01-27 [1] CRAN (R 4.3.0)
#>  phyloseq                   1.46.0     2023-10-24 [1] Bioconductor
#>  pillar                     1.9.0      2023-03-22 [1] CRAN (R 4.3.0)
#>  pkgbuild                   1.4.3      2023-12-10 [1] CRAN (R 4.3.0)
#>  pkgconfig                  2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
#>  pkgload                    1.3.4      2024-01-16 [1] CRAN (R 4.3.0)
#>  plyr                       1.8.9      2023-10-02 [1] CRAN (R 4.3.0)
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
#>  S4Vectors                * 0.40.2     2023-11-23 [1] Bioconductor
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
#>  SummarizedExperiment     * 1.32.0     2023-10-24 [1] Bioconductor
#>  survival                   3.5-5      2023-03-12 [1] CRAN (R 4.3.1)
#>  TH.data                    1.1-2      2023-04-17 [1] CRAN (R 4.3.0)
#>  tibble                   * 3.2.1      2023-03-20 [1] CRAN (R 4.3.0)
#>  tidyr                    * 1.3.1      2024-01-24 [1] CRAN (R 4.3.2)
#>  tidyselect                 1.2.0      2022-10-10 [1] CRAN (R 4.3.0)
#>  tidytree                   0.4.6      2023-12-12 [1] CRAN (R 4.3.0)
#>  tidyverse                * 2.0.0      2023-02-22 [1] CRAN (R 4.3.0)
#>  timechange                 0.3.0      2024-01-18 [1] CRAN (R 4.3.0)
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

