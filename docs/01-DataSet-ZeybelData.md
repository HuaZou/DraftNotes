


# (PART) Data Set {.unnumbered}

该书涉及到的数据集都会在本章节展示。


# Zeybel Dataset {#ZeybelDataset}


该数据集是来自于Zeybel 2022年发布的文章_Multiomics Analysis Reveals the Impact of Microbiota on Host Metabolism in Hepatic Steatosis_，它包含了多种组学数据，如

+ 微生物组（粪便和口腔）

+ 宿主人体学指标

+ 宿主临床学指标

+ 宿主血浆代谢组

+ 宿主血浆靶向免疫因子

+ 22位患者纵向时间序列数据


本脚本目的是生成符合数据分析的下游数据对象，主要有如下两类：


+ **phyloseq**: phyloseq-class object

  - otu table

  - taxa table

  - sample table

  - tree file

  - Representative sequences (ASV OTU)

+ **SummarizedExperiment**: SummarizedExperiment-class object

  - otu table

  - taxa table

  - sample table


## 加载R包

```r
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
library(dplyr)
library(tibble)
library(phyloseq)
library(SummarizedExperiment)

# rm(list = ls())
options(stringsAsFactors = F)
options(future.globals.maxSize = 1000 * 1024^2)
```


## Importing Data

可以点击此处下载数据[OmicsDataSet-Zeybel et al. - 2022.xlsx](https://github.com/HuaZou/DraftNotes/blob/main/InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx)

### cross-section data

56 participants


+ metadata

  - Patients' Phenotypic information (sheet 2)

  - Clinical and Physical variables (sheet 3)

+ taxa of fecal samples

  - Metaphlan2 profile of Gut microbiota (sheet 4)

  - Metaphlan2 profile of Oral microbiota (sheet 5)

+ metabolites of fecal samples

  - Metabolomics of Patients' plasma samples (sheet 6)

+ proteomics of host plasma

  - Inflammatory Proteomics of Patients' plasma samples (sheet 7)



```r
meta_info <- readxl::read_xlsx("./InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx", sheet = 2)
meta_clin <- readxl::read_xlsx("./InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx", sheet = 3)

taxa_gut <- readxl::read_xlsx("./InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx", sheet = 4)
taxa_ora <- readxl::read_xlsx("./InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx", sheet = 5)

metabolites_fecal <- readxl::read_xlsx("./InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx", sheet = 6)

protein_plasma <- readxl::read_xlsx("./InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx", sheet = 7)
```


### Longitudinal Data

22位患者在治疗前后的时间序列样本

+ metadata

  - Patients' Phenotypic information (sheet 8)

  - Clinical and Physical variables (sheet 9)

+ taxa of fecal samples

  - Metaphlan2 profile of Gut microbiota (sheet 10)

  - Metaphlan2 profile of Oral microbiota (sheet 11)

+ metabolites of fecal samples

  - Metabolomics of Patients' plasma samples (sheet 12)

+ proteomics of host plasma

  - Inflammatory Proteomics of Patients' plasma samples (sheet 13)



```r
long_meta_info <- readxl::read_xlsx("./InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx", sheet = 8)
long_meta_clin <- readxl::read_xlsx("./InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx", sheet = 9)

long_taxa_gut <- readxl::read_xlsx("./InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx", sheet = 10)
long_taxa_ora <- readxl::read_xlsx("./InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx", sheet = 11)

long_metabolites_fecal <- readxl::read_xlsx("./InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx", sheet = 12)

long_protein_plasma <- readxl::read_xlsx("./InputData/Zeybel-2022/OmicsDataSet-Zeybel et al. - 2022.xlsx", sheet = 13)
```



## Cross-section Data

56位患者的横截面也即是基线时期的数据

+ metadata

```r
metadata <- meta_info %>%
  dplyr::select(-`Patient ID`) %>%
  # dplyr::rename(`Patient ID` = `Patient ID...2`) %>%
  dplyr::inner_join(meta_clin,
                    by = "PatientID") %>%
  dplyr::select(-c("Stage", "Metabolomics", "Proteomics",
                   "GutMetagenomics", "OralMetagenomics")) %>%
  dplyr::select(PatientID, Gender, Age, everything())

head(metadata)
#> # A tibble: 6 × 47
#>   PatientID Gender   Age LiverFatClass AlcoholConsumption
#>   <chr>     <chr>  <dbl> <chr>         <chr>             
#> 1 P101001   Male      52 Severe        No                
#> 2 P101003   Female    31 None          No                
#> 3 P101004   Male      43 Moderate      Yes               
#> 4 P101007   Female    61 Severe        No                
#> 5 P101009   Male      51 Moderate      No                
#> 6 P101010   Male      27 Mild          Yes               
#> # ℹ 42 more variables: Smoker <chr>, Liver_fat <dbl>,
#> #   Sodium <dbl>, Potassium <dbl>, Creatinine <dbl>,
#> #   Urea_BUN <dbl>, Uric_acid <dbl>, ALT <dbl>, AST <dbl>,
#> #   GGT <dbl>, Alkaline_phosphatase <dbl>,
#> #   Total_Bilirubin <dbl>, Albumin <dbl>,
#> #   Creatine_Kinase <dbl>, Total_cholesterol <dbl>,
#> #   High_Density_Lipoprotein <dbl>, …
```

+ phyloseq-class object

```r
import_metaphlan_taxa <- function(data_metaphlan2,
                                  taxa_level = NULL) {

  taxa_rank <- c("Kingdom", "Phylum", "Class",
                 "Order", "Family", "Genus",
                 "Species", "Strain")

  # rename 1st column into "ID"
  colnames(data_metaphlan2)[which(colnames(data_metaphlan2) == "ID" |
                                    colnames(data_metaphlan2) == "clade_name")] <- "ID"

  # remove the "NCBI_tax_id" column
  if (any(colnames(data_metaphlan2) %in% "NCBI_tax_id")) {
    data_metaphlan2 <- data_metaphlan2 %>%
      dplyr::select(-NCBI_tax_id)
  }

  # remove sample names with "metaphlan_bugs_list" suffix
  colnames(data_metaphlan2) <- gsub("_metaphlan_bugs_list$", "",
                                    colnames(data_metaphlan2))

  if (is.null(taxa_level)) {
    taxa_level <- "Species"
  }
  ind_number <- match(taxa_level, taxa_rank)

  # remove k__Archaea & k__Viruses & k__Eukaryota & Others which are not bacteria
  if (length(grep("k__Bacteria", data_metaphlan2$ID)) > 0) {
    taxa_bacteria <- data_metaphlan2 %>%
      dplyr::slice(grep("k__Bacteria", ID)) %>%
      tibble::column_to_rownames("ID")
  } else {
    stop("There are no bacteria identified in this profile")
  }

  # dividing by 100 to normalize the sum into 1
  taxa_bac_per <- (taxa_bacteria / 100)  %>%
    tibble::rownames_to_column("ID")

  taxa <- taxa_bac_per %>%
    dplyr::group_by(ID) %>%
    dplyr::mutate(Number=length(unlist(strsplit(ID, "|", fixed = TRUE)))) %>%
    dplyr::select(ID, Number) %>%
    dplyr::filter(Number == ind_number)

  abundance_table <- taxa_bac_per %>%
    dplyr::filter(ID%in%taxa$ID)

  taxa_table <- lapply(abundance_table$ID, strsplit, split = "|", fixed = TRUE)
  taxa_table <- lapply(taxa_table, unlist)
  taxa_table <- do.call(rbind, taxa_table)
  taxa_table <- data.frame(taxa_table)
  names(taxa_table) <- taxa_rank[1:ind_number]
  rownames(taxa_table) <- taxa_table[, ind_number]

  abundance_table_res <- abundance_table %>%
    tibble::column_to_rownames("ID")
  rownames(abundance_table_res) <- rownames(taxa_table)
  abundance_table_res <- abundance_table_res[match(rownames(taxa_table),
                                                   row.names(abundance_table_res)), ,
                                             drop = FALSE]

  res <- list(tax_tab=taxa_table,
              abu_tab=abundance_table_res)

  return(res)
}

get_metaphlan_phyloseq <- function(otu_tab,
                                   sam_tab,
                                   tax_tab = NULL) {

  sid <- intersect(sam_tab$PatientID, colnames(otu_tab))

  sam_tab <- sam_tab[sam_tab$PatientID %in% sid, ] %>%
    tibble::column_to_rownames("PatientID")
  otu_tab <- otu_tab %>%
    dplyr::select(all_of(rownames(sam_tab)))


  if (!is.null(tax_tab)) {
    tax_tab <- phyloseq::tax_table(as(tax_tab, "matrix"))
  }

  if (!is.null(sam_tab)) {
    sam_tab <- phyloseq::sample_data(sam_tab %>% data.frame())
  }

  asv_tab <- phyloseq::otu_table(as(otu_tab, "matrix"),
                                 taxa_are_rows = TRUE)

  ps <- phyloseq::phyloseq(asv_tab, tax_tab, sam_tab)

  return(ps)
}


taxa_gut_list <- import_metaphlan_taxa(data_metaphlan2 = taxa_gut)
taxa_ora_list <- import_metaphlan_taxa(data_metaphlan2 = taxa_ora)

ps_gut <- get_metaphlan_phyloseq(otu_tab = taxa_gut_list$abu_tab,
                                 sam_tab = metadata,
                                 tax_tab = taxa_ora_list$tax_tab)
ps_gut
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 547 taxa and 42 samples ]
#> sample_data() Sample Data:       [ 42 samples by 46 sample variables ]
#> tax_table()   Taxonomy Table:    [ 547 taxa by 7 taxonomic ranks ]

ps_ora <- get_metaphlan_phyloseq(otu_tab = taxa_ora_list$abu_tab,
                                 sam_tab = metadata,
                                 tax_tab = taxa_ora_list$tax_tab)
ps_ora
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 547 taxa and 41 samples ]
#> sample_data() Sample Data:       [ 41 samples by 46 sample variables ]
#> tax_table()   Taxonomy Table:    [ 547 taxa by 7 taxonomic ranks ]
```


+ SummarizedExperiment: metabolites

```r
metabolites_fecal$metabolitesID <- paste0("Chem_", metabolites_fecal$CHEMICALID)

# Row->metabolites; Col->Samples
metabolites_fecal_assay <- metabolites_fecal %>%
  dplyr::select(c("metabolitesID", c(starts_with("P10")))) %>%
  tibble::column_to_rownames("metabolitesID")

metabolites_fecal_rowData <- metabolites_fecal %>%
  dplyr::select(-c(starts_with("P10"))) %>%
  dplyr::select(metabolitesID, everything())

metabolites_fecal_colData <- metadata  %>%
  dplyr::filter(PatientID %in% colnames(metabolites_fecal_assay)) %>%
  dplyr::arrange(PatientID)

metabolites_fecal_assay <- metabolites_fecal_assay %>%
  dplyr::select(all_of(metabolites_fecal_colData$PatientID))


se_metabolite <- SummarizedExperiment(
                     assays = metabolites_fecal_assay,
                     rowData = metabolites_fecal_rowData,
                     colData = metabolites_fecal_colData,
                     checkDimnames = TRUE)

se_metabolite
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

+ SummarizedExperiment: proteomics

```r
protein_plasma_new <- protein_plasma %>%
  dplyr::rename(ProteinID=SampleName)
protein_plasma_new$ProteinID <- gsub("-", "_", protein_plasma_new$ProteinID)

# Row->proteins; Col->Samples
protein_plasma_assay <- protein_plasma_new %>%
  dplyr::select(c("ProteinID", c(starts_with("P10")))) %>%
  tibble::column_to_rownames("ProteinID")

protein_plasma_rowData <- protein_plasma_new %>%
  dplyr::select(-c(starts_with("P10"))) %>%
  dplyr::select(ProteinID, everything())

protein_plasma_colData <- metadata %>%
  dplyr::filter(PatientID %in% colnames(protein_plasma_assay)) %>%
  dplyr::arrange(PatientID)

protein_plasma_assay <- protein_plasma_assay %>%
  dplyr::select(all_of(protein_plasma_colData$PatientID))

se_protein <- SummarizedExperiment(
                     assays = protein_plasma_assay,
                     rowData = protein_plasma_rowData,
                     colData = protein_plasma_colData,
                     checkDimnames = TRUE)

se_protein
#> class: SummarizedExperiment 
#> dim: 72 54 
#> metadata(0):
#> assays(1): ''
#> rownames(72): IL8 VEGFA ... TNFB CSF_1
#> rowData names(3): ProteinID LOD prop
#> colnames(54): P101001 P101003 ... P101095 P101096
#> colData names(47): PatientID Gender ...
#>   Right_leg_fat_free_mass Right_leg_total_body_water
```


+ output

```r
if (!dir.exists("./InputData/result/")) {
  dir.create("./InputData/result/", recursive = TRUE)
}

saveRDS(ps_gut, "./InputData/result/Zeybel_2022_gut_MGS_ps.RDS", compress = TRUE)
saveRDS(ps_ora, "./InputData/result/Zeybel_2022_oral_MGS_ps.RDS", compress = TRUE)
saveRDS(se_metabolite, "./InputData/result/Zeybel_2022_fecal_metabolite_se.RDS", compress = TRUE)
saveRDS(se_protein, "./InputData/result/Zeybel_2022_plasma_protein_se.RDS", compress = TRUE)
```



## longitudinal data

+ metadata

```r
# after
metadata_after <- long_meta_info %>%
  dplyr::inner_join(long_meta_clin,
                    by = "PatientID") %>%
  dplyr::select(PatientID, Gender, Age, everything())


# before
metadata_before <- metadata %>%
  dplyr::filter(PatientID %in% metadata_after$PatientID) %>%
  dplyr::mutate(Stage = "Before")

# cbind
long_metadata <- rbind(metadata_before, metadata_after) %>%
  dplyr::mutate(SampleID = paste(PatientID, Stage, sep = "_")) %>%
  dplyr::select(SampleID, PatientID, Gender, Age, Stage,
                LiverFatClass,
                Weight, Body_mass_index,
                Hip_circumference,
                Waist_circumference, everything()) %>%
  dplyr::arrange(PatientID, desc(Stage))

head(long_metadata)
#> # A tibble: 6 × 49
#>   SampleID PatientID Gender   Age Stage LiverFatClass Weight
#>   <chr>    <chr>     <chr>  <dbl> <chr> <chr>          <dbl>
#> 1 P101001… P101001   Male      52 Befo… Severe          83  
#> 2 P101001… P101001   Male      52 After Severe          82.8
#> 3 P101004… P101004   Male      43 Befo… Moderate        88  
#> 4 P101004… P101004   Male      43 After Moderate        88.6
#> 5 P101007… P101007   Female    61 Befo… Severe          95  
#> 6 P101007… P101007   Female    61 After Moderate        96.7
#> # ℹ 42 more variables: Body_mass_index <dbl>,
#> #   Hip_circumference <dbl>, Waist_circumference <dbl>,
#> #   AlcoholConsumption <chr>, Smoker <chr>,
#> #   Liver_fat <dbl>, Sodium <dbl>, Potassium <dbl>,
#> #   Creatinine <dbl>, Urea_BUN <dbl>, Uric_acid <dbl>,
#> #   ALT <dbl>, AST <dbl>, GGT <dbl>,
#> #   Alkaline_phosphatase <dbl>, Total_Bilirubin <dbl>, …
```

+ phyloseq-class object

```r
get_long_taxa <- function(x, y) {

  # x = long_taxa_gut
  # y = taxa_gut_list

  dat <- x
  otu_tab <- y$abu_tab
  tax_tab <- y$tax_tab

  dat$TaxaID <- paste0("s__", dat$TaxaID)
  dat <- dat %>%
    tibble::column_to_rownames("TaxaID")
  dat_cln <- dat[, !is.na(colSums(dat))]
  sid <- intersect(colnames(dat_cln), colnames(otu_tab))
  # before
  dat_before <- otu_tab %>%
    dplyr::select(all_of(sid))
  dat_after <- dat_cln %>%
    dplyr::select(all_of(sid))

  colnames(dat_before) <- paste0(colnames(dat_before), "_Before")
  colnames(dat_after) <- paste0(colnames(dat_after), "_After")

  mdat <- dat_before %>%
    tibble::rownames_to_column("TaxaID") %>%
    dplyr::inner_join(dat_after %>%
                        tibble::rownames_to_column("TaxaID"),
                      by = "TaxaID") %>%
    tibble::column_to_rownames("TaxaID")

  tax_tab_mdat <- tax_tab[rownames(tax_tab)%in%rownames(mdat), , ]


  res <- list(abu_tab = mdat,
              tax_tab = tax_tab_mdat)

  return(res)
}


get_metaphlan_phyloseq2 <- function(
    otu_tab,
    sam_tab,
    tax_tab = NULL) {

  sid <- intersect(sam_tab$SampleID, colnames(otu_tab))

  sam_tab <- sam_tab[sam_tab$SampleID %in% sid, ] %>%
    tibble::column_to_rownames("SampleID")
  otu_tab <- otu_tab %>%
    dplyr::select(all_of(rownames(sam_tab)))


  if (!is.null(tax_tab)) {
    tax_tab <- phyloseq::tax_table(as(tax_tab, "matrix"))
  }

  if (!is.null(sam_tab)) {
    sam_tab <- phyloseq::sample_data(sam_tab %>% data.frame())
  }

  asv_tab <- phyloseq::otu_table(as(otu_tab, "matrix"),
                                 taxa_are_rows = TRUE)

  ps <- phyloseq::phyloseq(asv_tab, tax_tab, sam_tab)

  return(ps)
}

long_taxa_gut_list <- get_long_taxa(x = long_taxa_gut, y = taxa_gut_list)
long_taxa_ora_list <- get_long_taxa(x = long_taxa_ora, y = taxa_ora_list)

long_ps_gut <- get_metaphlan_phyloseq2(
                      otu_tab = long_taxa_gut_list$abu_tab,
                      sam_tab = long_metadata,
                      tax_tab = long_taxa_gut_list$tax_tab)
long_ps_gut
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 547 taxa and 34 samples ]
#> sample_data() Sample Data:       [ 34 samples by 48 sample variables ]
#> tax_table()   Taxonomy Table:    [ 547 taxa by 7 taxonomic ranks ]

long_ps_ora <- get_metaphlan_phyloseq2(
                      otu_tab = long_taxa_ora_list$abu_tab,
                      sam_tab = long_metadata,
                      tax_tab = long_taxa_ora_list$tax_tab)
long_ps_ora
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 547 taxa and 36 samples ]
#> sample_data() Sample Data:       [ 36 samples by 48 sample variables ]
#> tax_table()   Taxonomy Table:    [ 547 taxa by 7 taxonomic ranks ]
```

+ SummarizedExperiment: metabolites

```r
metabolites_ID <- intersect(metabolites_fecal$BIOCHEMICAL, long_metabolites_fecal$MetaboliteID)
long_metabolites_fecal_rowData <- metabolites_fecal_rowData %>%
  dplyr::filter(BIOCHEMICAL %in% metabolites_ID)

# Row->metabolites; Col->Samples
long_metabolites_fecal_assay_after <- long_metabolites_fecal %>%
  dplyr::inner_join(long_metabolites_fecal_rowData %>%
                      dplyr::select(metabolitesID, BIOCHEMICAL),
                    by = c(MetaboliteID = "BIOCHEMICAL")) %>%
  dplyr::select(-MetaboliteID) %>%
  dplyr::select(c("metabolitesID", c(starts_with("P10")))) %>%
  tibble::column_to_rownames("metabolitesID")

meta_sid <- intersect(colnames(long_metabolites_fecal_assay_after),
                 colnames(metabolites_fecal_assay))
long_metabolites_fecal_assay_after <- long_metabolites_fecal_assay_after %>%
  dplyr::select(all_of(meta_sid))

long_metabolites_fecal_assay_before <-
  metabolites_fecal_assay[rownames(metabolites_fecal_assay) %in% rownames(long_metabolites_fecal_assay_after),
                          colnames(metabolites_fecal_assay) %in% colnames(long_metabolites_fecal_assay_after)]

colnames(long_metabolites_fecal_assay_before) <- paste0(colnames(long_metabolites_fecal_assay_before), "_Before")
colnames(long_metabolites_fecal_assay_after) <- paste0(colnames(long_metabolites_fecal_assay_after), "_After")

long_metabolites_fecal_assay <- long_metabolites_fecal_assay_before %>%
  tibble::rownames_to_column("TempRowNames") %>%
  dplyr::inner_join(long_metabolites_fecal_assay_after %>%
                      tibble::rownames_to_column("TempRowNames"),
                    by = "TempRowNames") %>%
  tibble::column_to_rownames("TempRowNames")

long_metabolites_fecal_colData <- long_metadata  %>%
  dplyr::filter(SampleID %in% colnames(long_metabolites_fecal_assay)) %>%
  dplyr::arrange(SampleID)


long_metabolites_fecal_assay <- long_metabolites_fecal_assay %>%
  dplyr::select(all_of(long_metabolites_fecal_colData$SampleID))


long_se_metabolite <- SummarizedExperiment(
                     assays = long_metabolites_fecal_assay,
                     rowData = long_metabolites_fecal_rowData,
                     colData = long_metabolites_fecal_colData,
                     checkDimnames = TRUE)

long_se_metabolite
#> class: SummarizedExperiment 
#> dim: 929 42 
#> metadata(0):
#> assays(1): ''
#> rownames(929): Chem_100002945 Chem_100002356 ...
#>   Chem_100015836 Chem_826
#> rowData names(13): metabolitesID BIOCHEMICAL ... KEGG
#>   SampleIDHMDBID
#> colnames(42): P101001_After P101001_Before ...
#>   P101077_After P101077_Before
#> colData names(49): SampleID PatientID ...
#>   Right_leg_fat_free_mass Right_leg_total_body_water
```

+ SummarizedExperiment: proteomics

```r
long_protein_plasma$ProteinID <- gsub("-", "_", long_protein_plasma$ProteinID)
proteins_ID <- intersect(protein_plasma_new$ProteinID, long_protein_plasma$ProteinID)

# Row->proteins; Col->Samples
long_protein_plasma_assay_after <- long_protein_plasma %>%
  dplyr::filter(ProteinID %in% proteins_ID) %>%
  dplyr::select(c(ProteinID, c(starts_with("P10")))) %>%
  tibble::column_to_rownames("ProteinID")

pron_id <- intersect(colnames(long_protein_plasma_assay_after),
                     colnames(protein_plasma_assay))
long_protein_plasma_assay_before <- protein_plasma_assay %>%
  dplyr::select(all_of(pron_id))
long_protein_plasma_assay_after <- long_protein_plasma_assay_after %>%
  dplyr::select(all_of(pron_id))

colnames(long_protein_plasma_assay_before) <- paste0(colnames(long_protein_plasma_assay_before), "_Before")
colnames(long_protein_plasma_assay_after) <- paste0(colnames(long_protein_plasma_assay_after), "_After")

long_protein_plasma_assay <- long_protein_plasma_assay_before %>%
  tibble::rownames_to_column("TempRowNames") %>%
  dplyr::inner_join(long_protein_plasma_assay_after %>%
                      tibble::rownames_to_column("TempRowNames"),
                    by = "TempRowNames") %>%
  tibble::column_to_rownames("TempRowNames")

long_protein_plasma_rowData <- protein_plasma_rowData %>%
  dplyr::rename(`Protein ID` = ProteinID)

long_protein_plasma_colData <- long_metadata  %>%
  dplyr::filter(SampleID %in% colnames(long_protein_plasma_assay)) %>%
  dplyr::arrange(SampleID)

long_protein_plasma_assay <- long_protein_plasma_assay %>%
  dplyr::select(all_of(long_protein_plasma_colData$SampleID))

long_se_protein <- SummarizedExperiment(
                     assays = long_protein_plasma_assay,
                     rowData = long_protein_plasma_rowData,
                     colData = long_protein_plasma_colData,
                     checkDimnames = TRUE)

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


## 保存数据

```r
if (!dir.exists("./InputData/result/")) {
  dir.create("./InputData/result/", recursive = TRUE)
}

saveRDS(long_ps_gut, "./InputData/result/Zeybel_2022_gut_MGS_ps_Paired.RDS", compress = TRUE)
saveRDS(long_ps_ora, "./InputData/result/Zeybel_2022_oral_MGS_ps_Paired.RDS", compress = TRUE)
saveRDS(long_se_metabolite, "./InputData/result/Zeybel_2022_fecal_metabolite_se_Paired.RDS", compress = TRUE)
saveRDS(long_se_protein, "./InputData/result/Zeybel_2022_plasma_protein_se_Paired.RDS", compress = TRUE)
```


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
#>  package              * version   date (UTC) lib source
#>  ade4                   1.7-22    2023-02-06 [2] CRAN (R 4.1.2)
#>  ape                    5.7-1     2023-03-13 [2] CRAN (R 4.1.2)
#>  Biobase              * 2.54.0    2021-10-26 [2] Bioconductor
#>  BiocGenerics         * 0.40.0    2021-10-26 [2] Bioconductor
#>  biomformat             1.22.0    2021-10-26 [2] Bioconductor
#>  Biostrings             2.62.0    2021-10-26 [2] Bioconductor
#>  bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.1.0)
#>  bookdown               0.34      2023-05-09 [2] CRAN (R 4.1.2)
#>  bslib                  0.6.0     2023-11-21 [1] CRAN (R 4.1.3)
#>  cachem                 1.0.8     2023-05-01 [2] CRAN (R 4.1.2)
#>  callr                  3.7.3     2022-11-02 [2] CRAN (R 4.1.2)
#>  cellranger             1.1.0     2016-07-27 [2] CRAN (R 4.1.0)
#>  cli                    3.6.1     2023-03-23 [2] CRAN (R 4.1.2)
#>  cluster                2.1.4     2022-08-22 [2] CRAN (R 4.1.2)
#>  codetools              0.2-19    2023-02-01 [2] CRAN (R 4.1.2)
#>  colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.1.2)
#>  crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.1.2)
#>  data.table             1.14.8    2023-02-17 [2] CRAN (R 4.1.2)
#>  DBI                    1.1.3     2022-06-18 [2] CRAN (R 4.1.2)
#>  DelayedArray           0.20.0    2021-10-26 [2] Bioconductor
#>  devtools               2.4.5     2022-10-11 [2] CRAN (R 4.1.2)
#>  digest                 0.6.33    2023-07-07 [1] CRAN (R 4.1.3)
#>  downlit                0.4.3     2023-06-29 [2] CRAN (R 4.1.3)
#>  dplyr                * 1.1.2     2023-04-20 [2] CRAN (R 4.1.2)
#>  ellipsis               0.3.2     2021-04-29 [2] CRAN (R 4.1.0)
#>  evaluate               0.21      2023-05-05 [2] CRAN (R 4.1.2)
#>  fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.1.2)
#>  fastmap                1.1.1     2023-02-24 [2] CRAN (R 4.1.2)
#>  foreach                1.5.2     2022-02-02 [2] CRAN (R 4.1.2)
#>  fs                     1.6.2     2023-04-25 [2] CRAN (R 4.1.2)
#>  generics               0.1.3     2022-07-05 [2] CRAN (R 4.1.2)
#>  GenomeInfoDb         * 1.30.1    2022-01-30 [2] Bioconductor
#>  GenomeInfoDbData       1.2.7     2022-03-09 [2] Bioconductor
#>  GenomicRanges        * 1.46.1    2021-11-18 [2] Bioconductor
#>  ggplot2                3.4.2     2023-04-03 [2] CRAN (R 4.1.2)
#>  glue                   1.6.2     2022-02-24 [2] CRAN (R 4.1.2)
#>  gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.1.2)
#>  htmltools              0.5.7     2023-11-03 [1] CRAN (R 4.1.3)
#>  htmlwidgets            1.6.2     2023-03-17 [2] CRAN (R 4.1.2)
#>  httpuv                 1.6.11    2023-05-11 [2] CRAN (R 4.1.3)
#>  igraph                 1.5.0     2023-06-16 [1] CRAN (R 4.1.3)
#>  IRanges              * 2.28.0    2021-10-26 [2] Bioconductor
#>  iterators              1.0.14    2022-02-05 [2] CRAN (R 4.1.2)
#>  jquerylib              0.1.4     2021-04-26 [2] CRAN (R 4.1.0)
#>  jsonlite               1.8.7     2023-06-29 [2] CRAN (R 4.1.3)
#>  knitr                  1.43      2023-05-25 [2] CRAN (R 4.1.3)
#>  later                  1.3.1     2023-05-02 [2] CRAN (R 4.1.2)
#>  lattice                0.21-8    2023-04-05 [2] CRAN (R 4.1.2)
#>  lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.1.2)
#>  magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.1.2)
#>  MASS                   7.3-60    2023-05-04 [2] CRAN (R 4.1.2)
#>  Matrix                 1.6-0     2023-07-08 [2] CRAN (R 4.1.3)
#>  MatrixGenerics       * 1.6.0     2021-10-26 [2] Bioconductor
#>  matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.1.3)
#>  memoise                2.0.1     2021-11-26 [2] CRAN (R 4.1.0)
#>  mgcv                   1.8-42    2023-03-02 [2] CRAN (R 4.1.2)
#>  mime                   0.12      2021-09-28 [2] CRAN (R 4.1.0)
#>  miniUI                 0.1.1.1   2018-05-18 [2] CRAN (R 4.1.0)
#>  multtest               2.50.0    2021-10-26 [2] Bioconductor
#>  munsell                0.5.0     2018-06-12 [2] CRAN (R 4.1.0)
#>  nlme                   3.1-162   2023-01-31 [2] CRAN (R 4.1.2)
#>  permute                0.9-7     2022-01-27 [2] CRAN (R 4.1.2)
#>  phyloseq             * 1.38.0    2021-10-26 [2] Bioconductor
#>  pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.1.2)
#>  pkgbuild               1.4.2     2023-06-26 [2] CRAN (R 4.1.3)
#>  pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.1.0)
#>  pkgload                1.3.2.1   2023-07-08 [2] CRAN (R 4.1.3)
#>  plyr                   1.8.8     2022-11-11 [2] CRAN (R 4.1.2)
#>  prettyunits            1.1.1     2020-01-24 [2] CRAN (R 4.1.0)
#>  processx               3.8.2     2023-06-30 [2] CRAN (R 4.1.3)
#>  profvis                0.3.8     2023-05-02 [2] CRAN (R 4.1.2)
#>  promises               1.2.0.1   2021-02-11 [2] CRAN (R 4.1.0)
#>  ps                     1.7.5     2023-04-18 [2] CRAN (R 4.1.2)
#>  purrr                  1.0.1     2023-01-10 [2] CRAN (R 4.1.2)
#>  R6                     2.5.1     2021-08-19 [2] CRAN (R 4.1.0)
#>  Rcpp                   1.0.11    2023-07-06 [1] CRAN (R 4.1.3)
#>  RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.1.2)
#>  readxl                 1.4.3     2023-07-06 [2] CRAN (R 4.1.3)
#>  remotes                2.4.2     2021-11-30 [2] CRAN (R 4.1.0)
#>  reshape2               1.4.4     2020-04-09 [2] CRAN (R 4.1.0)
#>  rhdf5                  2.38.1    2022-03-10 [2] Bioconductor
#>  rhdf5filters           1.6.0     2021-10-26 [2] Bioconductor
#>  Rhdf5lib               1.16.0    2021-10-26 [2] Bioconductor
#>  rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.1.2)
#>  rmarkdown              2.23      2023-07-01 [2] CRAN (R 4.1.3)
#>  rstudioapi             0.15.0    2023-07-07 [2] CRAN (R 4.1.3)
#>  S4Vectors            * 0.32.4    2022-03-29 [2] Bioconductor
#>  sass                   0.4.6     2023-05-03 [2] CRAN (R 4.1.2)
#>  scales                 1.2.1     2022-08-20 [2] CRAN (R 4.1.2)
#>  sessioninfo            1.2.2     2021-12-06 [2] CRAN (R 4.1.0)
#>  shiny                  1.7.4.1   2023-07-06 [2] CRAN (R 4.1.3)
#>  stringi                1.7.12    2023-01-11 [2] CRAN (R 4.1.2)
#>  stringr                1.5.0     2022-12-02 [2] CRAN (R 4.1.2)
#>  SummarizedExperiment * 1.24.0    2021-10-26 [2] Bioconductor
#>  survival               3.5-5     2023-03-12 [2] CRAN (R 4.1.2)
#>  tibble               * 3.2.1     2023-03-20 [2] CRAN (R 4.1.2)
#>  tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.1.2)
#>  urlchecker             1.0.1     2021-11-30 [2] CRAN (R 4.1.0)
#>  usethis                2.2.2     2023-07-06 [2] CRAN (R 4.1.3)
#>  utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.1.2)
#>  vctrs                  0.6.3     2023-06-14 [1] CRAN (R 4.1.3)
#>  vegan                  2.6-4     2022-10-11 [2] CRAN (R 4.1.2)
#>  withr                  2.5.0     2022-03-03 [2] CRAN (R 4.1.2)
#>  xfun                   0.40      2023-08-09 [1] CRAN (R 4.1.3)
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

+ [Multiomics Analysis Reveals the Impact of Microbiota on Host Metabolism in Hepatic Steatosis](https://onlinelibrary.wiley.com/doi/full/10.1002/advs.202104373)

+ [Basic storage, access, and manipulation of phylogenetic sequencing data with phyloseq](https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-basics.html)

+ [Practical: Organizing data with SummarizedExperiment](https://bioconductor.org/help/course-materials/2019/BSS2019/04_Practical_CoreApproachesInBioconductor.html)
