--- 
title: "生物信息随记"
author: "Hua Zou"
date: "2023-08-24 and updated on 2024-02-06"
site: bookdown::bookdown_site
documentclass: book
bibliography: [assets/book.bib, assets/packages.bib]
url: https://zouhua.top/DraftNotes/
#cover-image: images/cover.png
description: |
  记录平时做生信分析涉及到的一些知识点.
biblio-style: apalike
csl: chicago-fullnote-bibliography.csl
---


# 介绍


学习了很多东西，但是不成体系，又不知道记录在哪里，就随手创建了这个bookdown。暂时的想法是记录日常学习到的一些东西，后期有时间再重新梳理成博客，并将它们分门别类发布出来。


希望这里是自己学习的天堂～

加油～


## 输入数据

这些笔记的输入数据统一路径都在**[InputData](https://github.com/HuaZou/DraftNotes/blob/main/InputData/)**对应的子目录下面。如果想测试某些代码或者复现结果，可到对应的子目录内下载相关的数据。

+ **InputData**: https://github.com/HuaZou/DraftNotes/blob/main/InputData/

## 需要额外安装的R包

**[MicrobiomeAnalysis](https://zouhua.top/MicrobiomeAnalysis/)**提供一些数据预处理的方法如[impute_abundance](https://zouhua.top/MicrobiomeAnalysis/reference/impute_abundance.html)等

```R
if (!requireNamespace(c("remotes", "devtools"), quietly=TRUE)) {
  install.packages(c("devtools", "remotes"))
}

remotes::install_github("HuaZou/MicrobiomeAnalysis")

# library(MicrobiomeAnalysis)
```


