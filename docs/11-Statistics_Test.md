


# Hypothesis Testing Methods {#HypothesisTestingMethods}


显著性检验方法也即是假设检验方法。一般假设检验方法使用需要根据三个条件判断：

1. 数据分组数目也即是处理组数目，以2为阈值；

2. 数据的分布是否符合正态分布，符合则选择参数检验方法，否则选择非参数检验方法；

3. 数据是否是配对数据。


<div class="figure" style="text-align: center">
<img src="./InputData/figures/HypothesisTestingMethods/HypothesisTestingMethods.png" alt="Hypothesis Testing Methods" width="100%" height="100%" />
<p class="caption">(\#fig:unnamed-chunk-2)Hypothesis Testing Methods</p>
</div>


## 输入数据

对数据[OmicsDataSet-Zeybel et al. - 2022.xlsx](https://github.com/HuaZou/DraftNotes/blob/main/InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx)处理后生成的输入文件，详细情况可参考**Data Set**具体章节。

> ```R
> saveRDS(long_se_protein, "./InputData/result/Zeybel_2022_plasma_protein_se_Paired.RDS", compress = TRUE)
> ```

本次使用配对的血清蛋白质组学数据作为案例数据。



```r
knitr::opts_chunk$set(message = FALSE, warning = FALSE)

library(tidyverse)
library(SummarizedExperiment)
library(ggpubr)


long_se_protein <- readRDS("./InputData/result/Zeybel_2022_plasma_protein_se_Paired.RDS")

long_se_protein
#> class: SummarizedExperiment 
#> dim: 72 42 
#> metadata(0):
#> assays(1): ''
#> rownames(72): IL8 VEGFA ... TNFB CSF_1
#> rowData names(3): Protein ID LOD prop
#> colnames(42): P101001_After P101001_Before ...
#>   P101077_After P101077_Before
#> colData names(49): SampleID PatientID ...
#>   Right_leg_fat_free_mass Right_leg_total_body_water
```


+ 每组样本的数目 (配对分组：Mild， Moderate， Severe)


```r
with(colData(long_se_protein) %>% data.frame(), table(Stage, LiverFatClass))
#>         LiverFatClass
#> Stage    Mild Moderate None Severe
#>   After     5        7    2      7
#>   Before    6        7    0      8
```

+ 准备数据


```r
profile <- assay(long_se_protein) %>%
  as.data.frame()
metadata <- colData(long_se_protein) %>%
  as.data.frame()

# 两组不配对数据
metadata_2_unpaired <- metadata %>%
  dplyr::filter(Stage == "Before") %>%
  dplyr::filter(LiverFatClass %in% c("Mild", "Moderate"))
profile_2_unpaired <- profile[, pmatch(rownames(metadata_2_unpaired), colnames(profile)), F]
merge_2_unpaired <- metadata_2_unpaired %>%
  dplyr::select(SampleID, LiverFatClass) %>%
  dplyr::inner_join(profile_2_unpaired %>% 
                      t() %>%
                      as.data.frame() %>%
                      tibble::rownames_to_column("SampleID"),
                    by = "SampleID")


# 两组配对数据
metadata_2_paired <- metadata %>%
  dplyr::filter(LiverFatClass == "Moderate")
profile_2_paired <- profile[, pmatch(rownames(metadata_2_paired), colnames(profile)), F]
merge_2_paired <- metadata_2_paired %>%
  dplyr::select(SampleID, PatientID, Stage) %>%
  dplyr::inner_join(profile_2_paired %>% 
                      t() %>%
                      as.data.frame() %>%
                      tibble::rownames_to_column("SampleID"),
                    by = "SampleID") %>%
  dplyr::filter(!SampleID %in% c("P101052_After", "P101052_Before"))

# 三组数据
metadata_3_unpaired <- metadata %>%
  dplyr::filter(Stage == "Before") %>%
  dplyr::filter(LiverFatClass %in% c("Mild", "Moderate", "Severe"))
profile_3_unpaired <- profile[, pmatch(rownames(metadata_3_unpaired), colnames(profile)), F]
merge_3_unpaired <- metadata_3_unpaired %>%
  dplyr::select(SampleID, LiverFatClass) %>%
  dplyr::inner_join(profile_3_unpaired %>% 
                      t() %>%
                      as.data.frame() %>%
                      tibble::rownames_to_column("SampleID"),
                    by = "SampleID")
```


## 正态性评估

+ 可视化探索: density plot 密度图提供了一个关于分布是否呈钟形(正态分布)的直观判断

```r
ggdensity(merge_2_unpaired$IL8, 
          main = "Density plot of IL8",
          xlab = "IL8")
```

<img src="11-Statistics_Test_files/figure-html/unnamed-chunk-6-1.png" width="100%" />


+ 可视化探索: histogram 如果直方图大致呈“钟形”，则假定数据为正态分布

```r
gghistogram(merge_2_unpaired$IL8, 
          main = "Density plot of IL8",
          xlab = "IL8")
```

<img src="11-Statistics_Test_files/figure-html/unnamed-chunk-7-1.png" width="100%" />


+ 可视化探索: Q-Q plot Q-Q图描绘了给定样本与正态分布之间的相关性

```r
ggqqplot(merge_2_unpaired$IL8, 
          main = "Density plot of IL8",
          xlab = "IL8")
```

<img src="11-Statistics_Test_files/figure-html/unnamed-chunk-8-1.png" width="100%" />


+ 正态检验: `shapiro.test`提供单变量的正态分性检验方法（Shapiro-Wilk test）

> If the p-value of the test is greater than α = 0.05, then the data is assumed to be normally distributed.


```r
shapiro.test(merge_2_unpaired$IL8)
#> 
#> 	Shapiro-Wilk normality test
#> 
#> data:  merge_2_unpaired$IL8
#> W = 0.93269, p-value = 0.3694
```


+ 正态检验: `ks.test`提供单变量的正态分性检验方法（Kolmogorov-Smirnov test）

> If the p-value of the test is greater than α = 0.05, then the data is assumed to be normally distributed.


```r
ks.test(merge_2_unpaired$IL8,
        "pnorm")
#> 
#> 	One-sample Kolmogorov-Smirnov test
#> 
#> data:  merge_2_unpaired$IL8
#> D = 0.99996, p-value = 2.22e-16
#> alternative hypothesis: two-sided
```

结果：

+ 第二次检验的p值小于0.05，说明*IL8*数据不是正态分布。


## 非正态数据转换方法

如果给定的数据集不是正态分布，通常可以执行以下转换之一，使其更符合正态分布:

1. **Log Transformation**: 将$x$做$log(x)$转换

2. **Square Root Transformation**: 将$x$做$\sqrt{x}$开平方根转换

3. **Cube Root Transformation**: 将$x$做$x^{1/3}$开立方根转换


## Paired student's t-test

配对T检验适合两组数据且它们是正态分布和配对，计算t统计量。


![](./InputData/figures/HypothesisTestingMethods/t_statistics.jpg)


+ R基础函数`t.test`

```r
t.test(IL8 ~ Stage, data = merge_2_paired, paired = TRUE, alternative = "two.sided")
#> 
#> 	Paired t-test
#> 
#> data:  IL8 by Stage
#> t = 0.13325, df = 5, p-value = 0.8992
#> alternative hypothesis: true difference in means is not equal to 0
#> 95 percent confidence interval:
#>  -0.566393  0.628323
#> sample estimates:
#> mean of the differences 
#>                0.030965
```

+ `rstatix`提供的`t_test()`

```r
library(rstatix)

stat.test <- merge_2_paired  %>% 
  t_test(IL8 ~ Stage, paired = TRUE, detailed = TRUE) %>%
  add_significance()

stat.test
#> # A tibble: 1 × 14
#>   estimate .y.   group1 group2    n1    n2 statistic     p
#>      <dbl> <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl>
#> 1   0.0310 IL8   After  Before     6     6     0.133 0.899
#> # ℹ 6 more variables: df <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, method <chr>, alternative <chr>,
#> #   p.signif <chr>
```


p-value = 0.899比alpha = 0.05要大，接受零假设。IL8在前后没有显著差异。另外也可以计算Effect size效应值 (平均值和方差的比值)。

$d = \frac{mean_{D}}{SD_{D}}$


```r
merge_2_paired %>% cohens_d(IL8 ~ Stage, paired = TRUE)
#> # A tibble: 1 × 7
#>   .y.   group1 group2 effsize    n1    n2 magnitude 
#> * <chr> <chr>  <chr>    <dbl> <int> <int> <ord>     
#> 1 IL8   After  Before  0.0544     6     6 negligible
```


+ 可视化

```r
stat_label <- stat.test %>% add_xy_position(x = "Stage")

ggpaired(merge_2_paired, 
         x = "Stage", 
         y = "IL8", 
         order = c("Before", "After"),
         ylab = "IL8", 
         xlab = "Stage",
         fill = "Stage") + 
  stat_pvalue_manual(stat_label, tip.length = 0) +
  labs(subtitle = get_test_label(stat_label, detailed = TRUE))
```

<img src="11-Statistics_Test_files/figure-html/unnamed-chunk-14-1.png" width="100%" />


## Student's t-test

T检验适合两组数据且它们是正态分布，计算t统计量。


```r
t.test(IL8 ~ LiverFatClass, data = merge_2_unpaired, paired = FALSE, alternative = "two.sided")
#> 
#> 	Welch Two Sample t-test
#> 
#> data:  IL8 by LiverFatClass
#> t = -0.32904, df = 6.193, p-value = 0.753
#> alternative hypothesis: true difference in means between group Mild and group Moderate is not equal to 0
#> 95 percent confidence interval:
#>  -1.0116344  0.7702116
#> sample estimates:
#>     mean in group Mild mean in group Moderate 
#>               4.716680               4.837391

stat.test <- merge_2_unpaired  %>% 
  t_test(IL8 ~ LiverFatClass, paired = FALSE, detailed = TRUE) %>%
  add_significance()
# stat.test

merge_2_unpaired %>% cohens_d(IL8 ~ LiverFatClass, paired = FALSE)
#> # A tibble: 1 × 7
#>   .y.   group1 group2   effsize    n1    n2 magnitude 
#> * <chr> <chr>  <chr>      <dbl> <int> <int> <ord>     
#> 1 IL8   Mild   Moderate  -0.188     6     7 negligible
```

+ 可视化

```r
stat_label <- stat.test %>% add_xy_position(x = "LiverFatClass")

ggboxplot(merge_2_unpaired, 
         x = "LiverFatClass", 
         y = "IL8", 
         order = c("Mild", "Moderate"),
         ylab = "IL8", 
         xlab = "LiverFatClass",
         fill = "LiverFatClass") + 
  stat_pvalue_manual(stat_label, tip.length = 0) +
  labs(subtitle = get_test_label(stat_label, detailed = TRUE))
```

<img src="11-Statistics_Test_files/figure-html/unnamed-chunk-16-1.png" width="100%" />


## Wilcoxon signed rank t-test

配对Wilcoxon检验适合两组数据且它们是配对，对数据分布没有正态分布要求。

+ R基础函数`wilcox.test`

```r
wilcox.test(IL8 ~ Stage, data = merge_2_paired, paired = TRUE, alternative = "two.sided")
#> 
#> 	Wilcoxon signed rank exact test
#> 
#> data:  IL8 by Stage
#> V = 12, p-value = 0.8438
#> alternative hypothesis: true location shift is not equal to 0
```

+ `rstatix`提供的`wilcox_test()`

```r
library(rstatix)

stat.test <- merge_2_paired  %>% 
  wilcox_test(IL8 ~ Stage, paired = TRUE, detailed = TRUE) %>%
  add_significance()

stat.test
#> # A tibble: 1 × 13
#>   estimate .y.   group1 group2    n1    n2 statistic     p
#>      <dbl> <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl>
#> 1   0.0401 IL8   After  Before     6     6        12 0.844
#> # ℹ 5 more variables: conf.low <dbl>, conf.high <dbl>,
#> #   method <chr>, alternative <chr>, p.signif <chr>
```


## Mann-Whitney U test

> A Mann-Whitney U test (sometimes called the Wilcoxon rank-sum test) is used to compare the differences between two independent samples when the sample distributions are not normally distributed and the sample sizes are small (n <30).


```r
wilcox.test(IL8 ~ LiverFatClass, data = merge_2_unpaired, paired = FALSE, alternative = "two.sided")
#> 
#> 	Wilcoxon rank sum exact test
#> 
#> data:  IL8 by LiverFatClass
#> W = 13, p-value = 0.2949
#> alternative hypothesis: true location shift is not equal to 0
```

+ 可视化


```r
ggplot(merge_2_unpaired, aes(x = LiverFatClass, y = IL8)) + 
  geom_boxplot(width=0.3) +
  stat_summary(fun = mean, geom = "point", col = "black") +  
  stat_summary(fun = mean, geom = "text", col = "black", size = 3, 
               vjust = 3, aes(label = paste("Mean:", round(after_stat(y), digits = 2)))) +
  xlab("LiverFatClass") +
  ylab("IL8") +
  theme_bw()
```

<img src="11-Statistics_Test_files/figure-html/unnamed-chunk-20-1.png" width="100%" />


## Repeated measures One-way ANOVA

重复测量单因素方差分析需要满足以下假设：

1.在设计的任何block中都没有显著的异常值 (`rstatix::identify_outliers()`可查看离群点)

2.数据服从正态分布 (`rstatix::shapiro_test()`可查看检验结果)

3.方差齐性: 组间差异的方差应该相等 (`rstatix::anova_test()`可查看检验结果)

4.处理水平大于2


+ 输入数据


```r
data("selfesteem", package = "datarium")

selfesteem <- selfesteem %>%
  gather(key = "time", value = "score", t1, t2, t3) %>%
  convert_as_factor(id, time)

head(selfesteem, 3)
#> # A tibble: 3 × 3
#>   id    time  score
#>   <fct> <fct> <dbl>
#> 1 1     t1     4.01
#> 2 2     t1     2.56
#> 3 3     t1     3.24
```

+ 基本统计特征

```r
selfesteem %>%
  group_by(time) %>%
  get_summary_stats(score, type = "mean_sd")
#> # A tibble: 3 × 5
#>   time  variable     n  mean    sd
#>   <fct> <fct>    <dbl> <dbl> <dbl>
#> 1 t1    score       10  3.14 0.552
#> 2 t2    score       10  4.93 0.863
#> 3 t3    score       10  7.64 1.14
```

+ 可视化

```r
ggboxplot(selfesteem, x = "time", y = "score", add = "point")
```

<img src="11-Statistics_Test_files/figure-html/unnamed-chunk-23-1.png" width="100%" />

+ 假设评估：离群点

```r
selfesteem %>%
  group_by(time) %>%
  identify_outliers(score)
#> # A tibble: 2 × 5
#>   time  id    score is.outlier is.extreme
#>   <fct> <fct> <dbl> <lgl>      <lgl>     
#> 1 t1    6      2.05 TRUE       FALSE     
#> 2 t2    2      6.91 TRUE       FALSE
```

+ 假设评估：正态分布

```r
selfesteem %>%
  group_by(time) %>%
  shapiro_test(score)
#> # A tibble: 3 × 4
#>   time  variable statistic     p
#>   <fct> <chr>        <dbl> <dbl>
#> 1 t1    score        0.967 0.859
#> 2 t2    score        0.876 0.117
#> 3 t3    score        0.923 0.380

# ggqqplot(selfesteem, "score", facet.by = "time")
```

+ 统计检验：ANOVA整体评估变量在所有处理水平的显著性


```r
res.aov <- anova_test(data = selfesteem, dv = score, wid = id, within = time)
get_anova_table(res.aov)
#> ANOVA Table (type III tests)
#> 
#>   Effect DFn DFd      F        p p<.05   ges
#> 1   time   2  18 55.469 2.01e-08     * 0.829
```

结果：

1. p-value=2.01e-08，表明个人的score在不同时间点是统计学显著差异的。

2. ges: 广义效应大小(由受试者内部因素引起的可变性量)。


+ 统计检验：后置检验评估具体组间两两差异结果并做了检验结果校正


```r
pwc <- selfesteem %>%
  pairwise_t_test(
    score ~ time, paired = TRUE,
    p.adjust.method = "bonferroni"
    )
pwc
#> # A tibble: 3 × 10
#>   .y.   group1 group2    n1    n2 statistic    df          p
#> * <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl>      <dbl>
#> 1 score t1     t2        10    10     -4.97     9    7.72e-4
#> 2 score t1     t3        10    10    -13.2      9    3.34e-7
#> 3 score t2     t3        10    10     -4.87     9    8.86e-4
#> # ℹ 2 more variables: p.adj <dbl>, p.adj.signif <chr>
```

结果：

1. 任意两组间的pvalue均小于显著性水平alpha = 0.05。


+ 可视化：加上显著性标记

```r
pwc_label <- pwc %>% add_xy_position(x = "time")

ggboxplot(selfesteem, x = "time", y = "score", add = "point") + 
  stat_pvalue_manual(pwc_label) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_label))
```

<img src="11-Statistics_Test_files/figure-html/unnamed-chunk-28-1.png" width="100%" />


## One-way ANOVA

单因素方差分析需要满足以下假设：

1.数据服从正态分布 (`rstatix::shapiro_test()`可查看检验结果)

2.方差齐性: 组间差异的方差应该相等 (`rstatix::anova_test()`可查看检验结果)

3.处理水平大于2

+ 数据探索

```r
ggboxplot(PlantGrowth, x = "group", y = "weight",
          color = "group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("ctrl", "trt1", "trt2"),
          ylab = "Weight", xlab = "Treatment")
```

<img src="11-Statistics_Test_files/figure-html/unnamed-chunk-29-1.png" width="100%" />


+ 检验：评估植物的平均weight是否在三组处理间是显著差异的


```r
res.aov <- aov(weight ~ group, data = PlantGrowth)

summary(res.aov)
#>             Df Sum Sq Mean Sq F value Pr(>F)  
#> group        2  3.766  1.8832   4.846 0.0159 *
#> Residuals   27 10.492  0.3886                 
#> ---
#> Signif. codes:  
#> 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

结果：由于p值小于0.05的显著性水平，可以得出模型中标注“*”的组之间存在显著性差异。


+ 后置检验：组均值之间的多重两两比较 by Tukey HSD (Tukey Honest Significant Differences)


```r
TukeyHSD(res.aov)
#>   Tukey multiple comparisons of means
#>     95% family-wise confidence level
#> 
#> Fit: aov(formula = weight ~ group, data = PlantGrowth)
#> 
#> $group
#>             diff        lwr       upr     p adj
#> trt1-ctrl -0.371 -1.0622161 0.3202161 0.3908711
#> trt2-ctrl  0.494 -0.1972161 1.1852161 0.1979960
#> trt2-trt1  0.865  0.1737839 1.5562161 0.0120064
```


+ 后置检验2: 采用`multcomp::glht`方法

> The function glht() [in the multcomp package] can be used to do multiple comparison processes for an ANOVA. General linear hypothesis tests are abbreviated as glht.



```r
library(multcomp)

summary(glht(res.aov, linfct = mcp(group = "Tukey")))
#> 
#> 	 Simultaneous Tests for General Linear Hypotheses
#> 
#> Multiple Comparisons of Means: Tukey Contrasts
#> 
#> 
#> Fit: aov(formula = weight ~ group, data = PlantGrowth)
#> 
#> Linear Hypotheses:
#>                  Estimate Std. Error t value Pr(>|t|)  
#> trt1 - ctrl == 0  -0.3710     0.2788  -1.331   0.3909  
#> trt2 - ctrl == 0   0.4940     0.2788   1.772   0.1980  
#> trt2 - trt1 == 0   0.8650     0.2788   3.103   0.0122 *
#> ---
#> Signif. codes:  
#> 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> (Adjusted p values reported -- single-step method)
```


+ 后置检验3: T-test with pairs


```r
pairwise.t.test(PlantGrowth$weight, 
                PlantGrowth$group,
                p.adjust.method = "BH")
#> 
#> 	Pairwise comparisons using t tests with pooled SD 
#> 
#> data:  PlantGrowth$weight and PlantGrowth$group 
#> 
#>      ctrl  trt1 
#> trt1 0.194 -    
#> trt2 0.132 0.013
#> 
#> P value adjustment method: BH
```


+ 可视化：加上假设检验的结果


```r
pwc_label2 <- PlantGrowth %>%
  pairwise_t_test(
    weight ~ group,
    p.adjust.method = "BH"
    ) %>% add_xy_position(x = "group")

ggboxplot(PlantGrowth, x = "group", y = "weight", add = "point") + 
  stat_pvalue_manual(pwc_label2) +
  labs(
    #subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc_label2))
```

<img src="11-Statistics_Test_files/figure-html/unnamed-chunk-34-1.png" width="100%" />


## Friedman test

Friedman test，是一种非参数检验的方法，用于评估三个或更多成对组的分布之间是否存在统计学上的显著差异。当不满足单向重复测量ANOVA检验的正态性假设或因变量在有序量表上测量时，建议使用该方法。


+ 输入数据


```r
data("selfesteem", package = "datarium")

selfesteem <- selfesteem %>%
  gather(key = "time", value = "score", t1, t2, t3) %>%
  convert_as_factor(id, time)

head(selfesteem, 3)
#> # A tibble: 3 × 3
#>   id    time  score
#>   <fct> <fct> <dbl>
#> 1 1     t1     4.01
#> 2 2     t1     2.56
#> 3 3     t1     3.24
```

+ 统计检验：friedman test整体评估变量在所有处理水平的显著性


```r
res.fried <- selfesteem %>% friedman_test(score ~ time | id)

res.fried
#> # A tibble: 1 × 6
#>   .y.       n statistic    df        p method       
#> * <chr> <int>     <dbl> <dbl>    <dbl> <chr>        
#> 1 score    10      18.2     2 0.000112 Friedman test

# selfesteem %>% friedman_effsize(score ~ time | id)
```

+ 统计检验：后置检验评估具体组间两两差异结果并做了检验结果校正


```r
pwc <- selfesteem %>%
  wilcox_test(
    score ~ time, paired = TRUE,
    p.adjust.method = "bonferroni"
    )
pwc
#> # A tibble: 3 × 9
#>   .y.   group1 group2    n1    n2 statistic     p p.adj
#> * <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl>
#> 1 score t1     t2        10    10         0 0.002 0.006
#> 2 score t1     t3        10    10         0 0.002 0.006
#> 3 score t2     t3        10    10         1 0.004 0.012
#> # ℹ 1 more variable: p.adj.signif <chr>
```

+ 可视化：加上假设检验的结果


```r
pwc_label <- pwc %>% add_xy_position(x = "time")

ggboxplot(selfesteem, x = "time", y = "score", add = "point") + 
  stat_pvalue_manual(pwc_label, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.fried, detailed = TRUE),
    caption = get_pwc_label(pwc_label))
```

<img src="11-Statistics_Test_files/figure-html/unnamed-chunk-38-1.png" width="100%" />

## Kruskal-Wallis test

Kruskal-Wallis检验是单向方差分析检验的非参数替代检验。在有两个以上的组进行比较的情况下，它扩展了两样本Wilcoxon检验。当不满足单因素方差分析的假设时，建议使用。

+ 数据探索

```r
PlantGrowth %>% 
  group_by(group) %>%
  get_summary_stats(weight, type = "common")
#> # A tibble: 3 × 11
#>   group variable     n   min   max median   iqr  mean    sd
#>   <fct> <fct>    <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl>
#> 1 ctrl  weight      10  4.17  6.11   5.16 0.743  5.03 0.583
#> 2 trt1  weight      10  3.59  6.03   4.55 0.662  4.66 0.794
#> 3 trt2  weight      10  4.92  6.31   5.44 0.467  5.53 0.443
#> # ℹ 2 more variables: se <dbl>, ci <dbl>
```

+ 检验：评估植物的平均weight是否在三组处理间是显著差异的


```r
res.kruskal <- PlantGrowth %>% kruskal_test(weight ~ group)

res.kruskal
#> # A tibble: 1 × 6
#>   .y.        n statistic    df      p method        
#> * <chr>  <int>     <dbl> <int>  <dbl> <chr>         
#> 1 weight    30      7.99     2 0.0184 Kruskal-Wallis

# Effect size
# PlantGrowth %>% kruskal_effsize(weight ~ group)
```

结果：由于p值小于0.05的显著性水平，可以得出模型中标注“*”的组之间存在显著性差异。


+ 后置检验：组均值之间的多重两两比较 by Dunn’s test


```r
pwc <- PlantGrowth %>% 
  dunn_test(weight ~ group, p.adjust.method = "bonferroni") 
pwc
#> # A tibble: 3 × 9
#>   .y.    group1 group2    n1    n2 statistic       p  p.adj
#> * <chr>  <chr>  <chr>  <int> <int>     <dbl>   <dbl>  <dbl>
#> 1 weight ctrl   trt1      10    10     -1.12 0.264   0.791 
#> 2 weight ctrl   trt2      10    10      1.69 0.0912  0.273 
#> 3 weight trt1   trt2      10    10      2.81 0.00500 0.0150
#> # ℹ 1 more variable: p.adj.signif <chr>
```


+ 后置检验2: 采用`wilcox_test`方法



```r
pwc2 <- PlantGrowth %>% 
  wilcox_test(weight ~ group, p.adjust.method = "bonferroni")
pwc2
#> # A tibble: 3 × 9
#>   .y.    group1 group2    n1    n2 statistic     p p.adj
#> * <chr>  <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl>
#> 1 weight ctrl   trt1      10    10      67.5 0.199 0.597
#> 2 weight ctrl   trt2      10    10      25   0.063 0.189
#> 3 weight trt1   trt2      10    10      16   0.009 0.027
#> # ℹ 1 more variable: p.adj.signif <chr>
```

+ 可视化：加上假设检验的结果


```r
pwc_label2 <- pwc %>% add_xy_position(x = "group")

ggboxplot(PlantGrowth, x = "group", y = "weight", add = "point") + 
  stat_pvalue_manual(pwc_label2, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc_label2))
```

<img src="11-Statistics_Test_files/figure-html/unnamed-chunk-43-1.png" width="100%" />


## Blocked Wilcoxon rank-sum test

**Two-sided Wilcoxon tests blocked for ‘study’**是Wilcoxon检验是在考虑不同研究来源（study）的影响下进行的差异检验。这不同于单纯地在每个研究内部分别进行Wilcoxon检验，因为它试图控制或调整来自不同研究的潜在影响，从而提供更准确和可靠的整体分析结果。


+ **formula**: a formula of the form y ~ x | block where y is a numeric variable, x is a factor and block is an optional factor for stratification.


```r
# 安装并加载coin包
library(coin)

# 示例数据，确保group和study列均为因子类型
data <- data.frame(
  variable = rnorm(100),
  group = factor(rep(c("A", "B"), 50)),
  study = factor(rep(c("Study1", "Study2"), each = 50))
)

# 使用wilcox_test进行分组Wilcoxon检验
result <- coin::wilcox_test(variable ~ group | study, data = data)
print(result)
#> 
#> 	Asymptotic Wilcoxon-Mann-Whitney Test
#> 
#> data:  variable by
#> 	 group (A, B) 
#> 	 stratified by study
#> Z = 0.60445, p-value = 0.5455
#> alternative hypothesis: true mu is not equal to 0
```



## 总结

1.两组使用t-test或wilcox-test，前者适合正态分布数据后者为非参数检验方法。

2.三组及以上使用ANOVA, friedman或KW检验，前者适合正态分布数据后者为非参数检验方法。

3.三组数据的初步检验结果需要做后置检验才能解析出具体组间差异。

4.在假设检验前，可以对数据进行探索，如正态性评估以及箱线图组间分布评估。


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
#>  package              * version   date (UTC) lib source
#>  abind                  1.4-5     2016-07-21 [2] CRAN (R 4.1.0)
#>  backports              1.4.1     2021-12-13 [2] CRAN (R 4.1.0)
#>  Biobase              * 2.54.0    2021-10-26 [2] Bioconductor
#>  BiocGenerics         * 0.40.0    2021-10-26 [2] Bioconductor
#>  bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.1.0)
#>  bookdown               0.34      2023-05-09 [2] CRAN (R 4.1.2)
#>  broom                  1.0.5     2023-06-09 [2] CRAN (R 4.1.3)
#>  bslib                  0.6.0     2023-11-21 [1] CRAN (R 4.1.3)
#>  cachem                 1.0.8     2023-05-01 [2] CRAN (R 4.1.2)
#>  callr                  3.7.3     2022-11-02 [2] CRAN (R 4.1.2)
#>  car                    3.1-2     2023-03-30 [2] CRAN (R 4.1.2)
#>  carData                3.0-5     2022-01-06 [2] CRAN (R 4.1.2)
#>  cli                    3.6.1     2023-03-23 [2] CRAN (R 4.1.2)
#>  codetools              0.2-19    2023-02-01 [2] CRAN (R 4.1.2)
#>  coin                 * 1.4-2     2021-10-08 [2] CRAN (R 4.1.0)
#>  colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.1.2)
#>  crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.1.2)
#>  DelayedArray           0.20.0    2021-10-26 [2] Bioconductor
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
#>  fs                     1.6.2     2023-04-25 [2] CRAN (R 4.1.2)
#>  generics               0.1.3     2022-07-05 [2] CRAN (R 4.1.2)
#>  GenomeInfoDb         * 1.30.1    2022-01-30 [2] Bioconductor
#>  GenomeInfoDbData       1.2.7     2022-03-09 [2] Bioconductor
#>  GenomicRanges        * 1.46.1    2021-11-18 [2] Bioconductor
#>  ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.1.2)
#>  ggpubr               * 0.6.0     2023-02-10 [2] CRAN (R 4.1.2)
#>  ggsignif               0.6.4     2022-10-13 [2] CRAN (R 4.1.2)
#>  glue                   1.6.2     2022-02-24 [2] CRAN (R 4.1.2)
#>  gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.1.2)
#>  highr                  0.10      2022-12-22 [2] CRAN (R 4.1.2)
#>  hms                    1.1.3     2023-03-21 [2] CRAN (R 4.1.2)
#>  htmltools              0.5.7     2023-11-03 [1] CRAN (R 4.1.3)
#>  htmlwidgets            1.6.2     2023-03-17 [2] CRAN (R 4.1.2)
#>  httpuv                 1.6.11    2023-05-11 [2] CRAN (R 4.1.3)
#>  IRanges              * 2.28.0    2021-10-26 [2] Bioconductor
#>  jquerylib              0.1.4     2021-04-26 [2] CRAN (R 4.1.0)
#>  jsonlite               1.8.7     2023-06-29 [2] CRAN (R 4.1.3)
#>  knitr                  1.43      2023-05-25 [2] CRAN (R 4.1.3)
#>  labeling               0.4.2     2020-10-20 [2] CRAN (R 4.1.0)
#>  later                  1.3.1     2023-05-02 [2] CRAN (R 4.1.2)
#>  lattice                0.21-8    2023-04-05 [2] CRAN (R 4.1.2)
#>  libcoin                1.0-9     2021-09-27 [2] CRAN (R 4.1.0)
#>  lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.1.2)
#>  lubridate            * 1.9.2     2023-02-10 [2] CRAN (R 4.1.2)
#>  magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.1.2)
#>  MASS                 * 7.3-60    2023-05-04 [2] CRAN (R 4.1.2)
#>  Matrix                 1.6-0     2023-07-08 [2] CRAN (R 4.1.3)
#>  MatrixGenerics       * 1.6.0     2021-10-26 [2] Bioconductor
#>  matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.1.3)
#>  memoise                2.0.1     2021-11-26 [2] CRAN (R 4.1.0)
#>  mime                   0.12      2021-09-28 [2] CRAN (R 4.1.0)
#>  miniUI                 0.1.1.1   2018-05-18 [2] CRAN (R 4.1.0)
#>  modeltools             0.2-23    2020-03-05 [2] CRAN (R 4.1.0)
#>  multcomp             * 1.4-25    2023-06-20 [2] CRAN (R 4.1.3)
#>  munsell                0.5.0     2018-06-12 [2] CRAN (R 4.1.0)
#>  mvtnorm              * 1.2-2     2023-06-08 [2] CRAN (R 4.1.3)
#>  pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.1.2)
#>  pkgbuild               1.4.2     2023-06-26 [2] CRAN (R 4.1.3)
#>  pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.1.0)
#>  pkgload                1.3.2.1   2023-07-08 [2] CRAN (R 4.1.3)
#>  prettyunits            1.1.1     2020-01-24 [2] CRAN (R 4.1.0)
#>  processx               3.8.2     2023-06-30 [2] CRAN (R 4.1.3)
#>  profvis                0.3.8     2023-05-02 [2] CRAN (R 4.1.2)
#>  promises               1.2.0.1   2021-02-11 [2] CRAN (R 4.1.0)
#>  ps                     1.7.5     2023-04-18 [2] CRAN (R 4.1.2)
#>  purrr                * 1.0.1     2023-01-10 [2] CRAN (R 4.1.2)
#>  R6                     2.5.1     2021-08-19 [2] CRAN (R 4.1.0)
#>  Rcpp                   1.0.11    2023-07-06 [1] CRAN (R 4.1.3)
#>  RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.1.2)
#>  readr                * 2.1.4     2023-02-10 [2] CRAN (R 4.1.2)
#>  remotes                2.4.2     2021-11-30 [2] CRAN (R 4.1.0)
#>  rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.1.2)
#>  rmarkdown              2.23      2023-07-01 [2] CRAN (R 4.1.3)
#>  rstatix              * 0.7.2     2023-02-01 [2] CRAN (R 4.1.2)
#>  rstudioapi             0.15.0    2023-07-07 [2] CRAN (R 4.1.3)
#>  S4Vectors            * 0.32.4    2022-03-29 [2] Bioconductor
#>  sandwich               3.0-2     2022-06-15 [2] CRAN (R 4.1.2)
#>  sass                   0.4.6     2023-05-03 [2] CRAN (R 4.1.2)
#>  scales                 1.2.1     2022-08-20 [2] CRAN (R 4.1.2)
#>  sessioninfo            1.2.2     2021-12-06 [2] CRAN (R 4.1.0)
#>  shiny                  1.7.4.1   2023-07-06 [2] CRAN (R 4.1.3)
#>  stringi                1.7.12    2023-01-11 [2] CRAN (R 4.1.2)
#>  stringr              * 1.5.0     2022-12-02 [2] CRAN (R 4.1.2)
#>  SummarizedExperiment * 1.24.0    2021-10-26 [2] Bioconductor
#>  survival             * 3.5-5     2023-03-12 [2] CRAN (R 4.1.2)
#>  TH.data              * 1.1-2     2023-04-17 [2] CRAN (R 4.1.2)
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
#>  withr                  2.5.0     2022-03-03 [2] CRAN (R 4.1.2)
#>  xfun                   0.40      2023-08-09 [1] CRAN (R 4.1.3)
#>  xml2                   1.3.5     2023-07-06 [2] CRAN (R 4.1.3)
#>  xtable                 1.8-4     2019-04-21 [2] CRAN (R 4.1.0)
#>  XVector                0.34.0    2021-10-26 [2] Bioconductor
#>  yaml                   2.3.7     2023-01-23 [2] CRAN (R 4.1.2)
#>  zlibbioc               1.40.0    2021-10-26 [2] Bioconductor
#>  zoo                    1.8-12    2023-04-13 [2] CRAN (R 4.1.2)
#> 
#>  [1] /Users/zouhua/Library/R/x86_64/4.1/library
#>  [2] /Library/Frameworks/R.framework/Versions/4.1/Resources/library
#> 
#> ──────────────────────────────────────────────────────────
```


## Reference

+ [data-science-for-beginners](https://bookdown.org/BaktiSiregar/data-science-for-beginners-part-2/)

