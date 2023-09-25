---
title: "Rumen_Compar"
author: "Gerardo Diaz"
date: "`r Sys.Date()`"
output: html_document
editor_options: 
  chunk_output_type: console
---

# 0. Loading packages
Code to load all the needed packages:
```{r libraries, warning=FALSE, message=FALSE}
library("phyloseq")
library("forcats")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("microbiome")   # Extra data manipulation
library(metagenomeSeq)  # Data normalization (?) + Statistical testing
library(limma)          # Statistical testing
library("statmod")      # Will work with limma for "duplicateCorrelation" function
library(psych)          # Data exploration
library("tidyverse")    # ggplot2, dplyr, tidyr, forcats, tibble, string, purrr, readr
library(vegan)
library(lme4)           # Statistical testing
library(lmerTest)       # Statistical testing
library(emmeans)        # Statistical testing
library(car)            # Statistical testing
library(table1)         # for doing tables
library(lemon)          # for volcano plot
library("ggrepel")      # for volcano plot
library(decontam)       # decontamination using extraction blank samples reads
library(ade4)           # for procrustes 
library(grid)           # for plotting personalized procrustes plot
library(UpSetR)         # for upset plots of common taxa between methods
#library("MicrobiotaProcess") # only activate for get_upset
 library(ggExtra)       # for marginal histograms in scatterplot
```

# 1. Loading RAW data 
Change your working directory to where the files are located
Define you 3 tables to build your phyloseq object:
* OTU
* Taxonomy
* Samples
```{r}
  otu_mat<- read_excel("RumenDATA_MASTER.xlsx", sheet = "OTU_matrix_1")
  tax_mat<- read_excel("RumenDATA_MASTER.xlsx", sheet = "Taxonomy_table_1")
  samples_df<- read_excel("RumenDATA_MASTER.xlsx", sheet = "Samples")
```

Give some format to the sample table (re-leveling, etc)
```{r}
# formatting and re-leveling categorical variables 
  samples_df$collection_date= factor(samples_df$collection_date, levels=c("9.21.21", "10.18.21", "10.20.21"))
  samples_df$collection_day= factor(samples_df$collection_day, levels=c("Pre_weaning", "At_weaning", "Post_weaning"))
   samples_df$castration_group= factor(samples_df$castration_group, levels=c("Birth", "Turnout", "Pre_weaning","Weaning"))
   samples_df$castration= factor(samples_df$castration, levels=c("Birth", "Turnout", "Pre_weaning","Weaning", "Not_castrated"))
   samples_df$castration_status= factor(samples_df$castration_status, levels = c("Castrated", "Not_castrated"))
   samples_df$weaning= factor(samples_df$weaning, levels=c("Fence_line", "Truck","Not_weaned"))
   samples_df$weaning_group=factor(samples_df$weaning_group, levels=c("Fence_line", "Truck"))
   samples_df$weaning_status=factor(samples_df$weaning_status, levels=c("Weaned", "Not_weaned"))
  samples_df$cow_ID= as.factor(samples_df$cow_ID)
  samples_df$extraction_date= as.factor(samples_df$extraction_date)
  samples_df$extraction_run= as.factor(samples_df$extraction_run)
```

## 1.1. Phyloseq object
Create the RAW phyloseq object
```{r}
# Define the row names from the otu column in otu_mat and tax_mat
  otu_mat <- otu_mat %>%
    tibble::column_to_rownames("otu") 
  tax_mat <- tax_mat %>% 
    tibble::column_to_rownames("otu")
   samples_df <- samples_df %>% 
    tibble::column_to_rownames("samples")
# Transform into matrixes otu and tax tables (sample table can be left as data frame)
  otu_mat <- as.matrix(otu_mat)
  tax_mat <- as.matrix(tax_mat)
# Transform to phyloseq objects
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  samples = sample_data(samples_df)
  ps.object <- phyloseq(OTU, TAX, samples)
```

Give a name to the ps.object
```{r}
# Change the left part of the equivalence 
  R_K1 = ps.object
```

## 1.2. Decontamination of 16S and dropping NA

We will use *Decontam* package. For more information please see: 
https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
https://github.com/benjjneb/decontam
```{r}

Rumen = R_16S

# Method: FREQUENCY
# Not prevalence because I have just 3 blank samples and 2 of them did not register any classified read

# STEP 1: Inspect your reads for every sample (sample, mock and blank negative control)
df=as.data.frame(sample_data(Rumen))

str(df$seq_qPCR)
hist(sample_data(Rumen)$seq_qPCR)

df$LibrarySize=sample_sums(Rumen)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_type)) + geom_point()

# STEP 2: Run decontam using "Frequency" method and inspect your results
contamdf.freq <- isContaminant(Rumen, method="frequency", conc="seq_qPCR")
head(contamdf.freq)
table(contamdf.freq$contaminant)

# STEP 3: Trying to plot the model to confirm that my contaminants were REALLY contaminats (DID NOT WORK BECAUSE I HAVE 2 NEGATIVE CONTROLS WITH 0 READS)
set.seed(999)
plot_frequency(Rumen, taxa_names(Rumen)[sample(which(contamdf.freq$contaminant),2)],conc="seq_qPCR") 
sum(sample_sums(Rumen))

# STEP 4: Getting your DECONTAMINATED Phyloseq object 
Rumen_DC <- prune_taxa(!contamdf.freq$contaminant, Rumen) # you can get the phyloseq object of your contaminants just taking out the "!" symbol
Rumen_DC 

# STEP 5: Getting the number and percentage of total reads removed with the contaminants
sum(sample_sums(Rumen_DC))
#x=5053/5868161
#x*100

R_16S = Rumen_DC

```

Dropping NA from 16S
```{r}
  R_16S <- subset_taxa(R_16S, !Kingdom %in% c("NA"))
```

## 1.3. Dropping viruses, viroids, eukaryota reads from kraken classification
We will just leave bacteria and archaea classifications to have a fear comparison
```{r}
  R_K0 <- subset_taxa(R_K0, !Domain %in% c("Viruses","Viroids","Eukaryota")) # Raw 11588 OTU, clean 4706 OTU
  R_K0.1 <- subset_taxa(R_K0.1, !Domain %in% c("Viruses","Viroids","Eukaryota")) # Raw 5206 OTU, clean 4635 OTU
  R_K1 <- subset_taxa(R_K1, !Domain %in% c("Viruses","Viroids","Eukaryota")) # Raw 2796 OTU, clean 2772 OTU
```

*Clean phyloseq object
R_16S: 7653 taxa
R_K0: 4706 taxa
R_K0.1: 4635 taxa
R_K1: 2772 taxa

# 2. Summary statistics
Just to explore what is going on
```{r}
samples_2<-subset(samples_df, samples_df$sample_type=="Sample")
samples_2=data.frame(samples_2)
samples_D1<-subset(samples_2, samples_2$collection_day=="Pre_weaning")
samples_D1=data.frame(samples_D1)
samples_D2<-subset(samples_2, samples_2$collection_day=="At_weaning")
samples_D2=data.frame(samples_D2)
samples_D3<-subset(samples_2, samples_2$collection_day=="Post_weaning")
samples_D3=data.frame(samples_D3)

##################################################################
# kruskal.test(samples_D3$age_days~samples_D3$weaning)           #
# wilcox.test(samples_D3$age_days~samples_D3$weaning)            #
# t.test(samples_D3$age_days~samples_D3$weaning)                 #
# ANOVA<-aov(samples_D3$age_days~samples_D3$weaning)             #
# summary(ANOVA)                                                 #
##################################################################

# SEQUENCING READS
# Raw reads for sample type
describeBy(samples_df$seq_16SRAWreads, samples_df$sample_type)
describeBy(samples_df$seq_RAWreads, samples_df$sample_type)
# 16S
describeBy(samples_2$seq_qPCR, samples_2$collection_day)
ANOVA<-aov(samples_2$seq_qPCR~samples_2$collection_day)
summary(ANOVA)
kruskal.test(samples_2$seq_qPCR, samples_2$collection_day)
describeBy(samples_2$seq_16SRAWreads, samples_2$collection_day)
ANOVA<-aov(samples_2$seq_16SRAWreads~samples_2$collection_day)
summary(ANOVA)
kruskal.test(samples_2$seq_16SRAWreads, samples_2$collection_day)
describeBy(samples_2$seq_16SdenoisedANDmerged, samples_2$collection_day)
ANOVA<-aov(samples_2$seq_16SdenoisedANDmerged~samples_2$collection_day)
summary(ANOVA)
kruskal.test(samples_2$seq_16SdenoisedANDmerged, samples_2$collection_day)
describeBy(samples_2$seq_16Sclassified, samples_2$collection_day)
ANOVA<-aov(samples_2$seq_16Sclassified~samples_2$collection_day)
summary(ANOVA)
kruskal.test(samples_2$seq_16Sclassified, samples_2$collection_day)
# shotgun
describeBy(samples_2$seq_RAWreads, samples_2$collection_day)
ANOVA<-aov(samples_2$seq_RAWreads~samples_2$collection_day)
summary(ANOVA)
kruskal.test(samples_2$seq_RAWreads, samples_2$collection_day)
describeBy(samples_2$seq_nonHOSTreads, samples_2$collection_day)
ANOVA<-aov(samples_2$seq_nonHOSTreads~samples_2$collection_day)
summary(ANOVA)
kruskal.test(samples_2$seq_nonHOSTreads, samples_2$collection_day)
describeBy(samples_2$seq_K0classified, samples_2$collection_day)
ANOVA<-aov(samples_2$seq_K0classified~samples_2$collection_day)
summary(ANOVA)
kruskal.test(samples_2$seq_K0classified, samples_2$collection_day)
describeBy(samples_2$seq_K1classified, samples_2$collection_day)
ANOVA<-aov(samples_2$seq_K1classified~samples_2$collection_day)
summary(ANOVA)
kruskal.test(samples_2$seq_K1classified, samples_2$collection_day)

pairwise.t.test(samples_2$seq_K1classified,samples_2$collection_day, p.adj = "bonferroni")
```

# 3. Phyloseq object preparation
We will:
* Subset data: Samples (n=95) 
* Normalize (using the Cumulative Sum Scaling) the counts for all our samples to make the groups comparable between each other. We will normalize the negative and positive controls as well only to present that result.
* Our normalized samples will be subset only of: (1) Bacteria and (2) Archaea. 

## 3.1. Normalization (CSS)
Subset phyloseq object first so we only have the Samples (n=95) to CSS them
Obtain:Sample (95 samples from 3 collection days) or by each collection day
```{r}
  R_16S_Sample <- subset_samples(R_16S, sample_type =="Sample")
  R_K0_Sample <- subset_samples(R_K0, sample_type =="Sample")
  R_K0.1_Sample <- subset_samples(R_K0.1, sample_type =="Sample")
  R_K1_Sample <- subset_samples(R_K1, sample_type =="Sample")
```
Subset the outliers from all the phyloseq objects: R0015 and R0127
```{r}
  R_16S_Sample_clean <- prune_samples(!(sample_names(R_16S_Sample) %in% c("R0015", "R0127")), R_16S_Sample)
  R_K0_Sample_clean <- prune_samples(!(sample_names(R_K0_Sample) %in% c("R0015", "R0127")), R_K0_Sample)
  R_K0.1_Sample_clean <- prune_samples(!(sample_names(R_K0.1_Sample) %in% c("R0015", "R0127")), R_K0.1_Sample)
  R_K1_Sample_clean <- prune_samples(!(sample_names(R_K1_Sample) %in% c("R0015", "R0127")), R_K1_Sample)
```

Method Cumulative Sum Scaling (using metagenomeSeq package) 
Define the phyloseq object that you want to normalize: 
```{r}
# Change the right part of the equivalence 
  ps = R_K0_Sample_clean
```
Obtain the CSS normalized object:
```{r}
# Changing from phyloseq object to MRexperiment (metagenomeSeq format) for CSS normalization
rumen.metaseq <- phyloseq_to_metagenomeSeq(ps) # change here

# CSS normalization (cumNorm)
# cumNormStat: Calculates the percentile for which to sum counts up to and scale by. 
rumen.metaseq.norm<- cumNorm(rumen.metaseq, p=cumNormStat(rumen.metaseq))

# Extracting (MRcounts) the normalized table to be used as otu_table (norm=TRUE)
CSS_rumen.metaseq <- MRcounts(rumen.metaseq.norm, norm = TRUE)

# Merging the extracted normalized count table to use as otu_table in normalized phyloseq object "CSS_Rumen"
CSS_Rumen <- merge_phyloseq(otu_table(CSS_rumen.metaseq,
taxa_are_rows=T),sample_data(ps),tax_table(ps))
```
Recover your table
```{r}
# change the left part of the table
  CSS_K0_Sample_clean = CSS_Rumen
```
We obtained: CSS_16S_Sample; CSS_K0_Sample; CSS_K0.1_Sample; CSS_K1_Sample
Without R0015 and R0127: CSS_16S_Sample_clean; CSS_K0_Sample_clean; CSS_K0.1_Sample_clean; CSS_K1_Sample_clean

## 3.2. Agglomerate taxa to different levels
Define the ps object that you will agglomerate 
Also define the level you will agglomerate
```{r}
# Change the right part of the equivalence 
  ps = R_K1_Sample_clean
  lvl = "Phylum"
```
Agglomerate taxa:
```{r}
 # glom= tax_glom(ps, lvl)
#glom= tax_glom(ps, lvl,NArm = F)
# glom= tax_glom(ps, lvl, bad_empty = c("NA", NA)) # Remove NA from kraken
 glom= tax_glom(ps, lvl, bad_empty = c("NA", NA, "", " ", "\t")) # Remove NA from 16S
# NOTE: To obtain the number of classified reads at any given level, just agglomerate Rumen_Sample object to the level you want to investigate and then do: 
# get how many reads where classified to a given taxonomic rank
  sum(taxa_sums(glom))
```
Recover object -- Give a proper name to the agglomerated object according to the kingdom/domain you agglomerate and the level you agglomerated 
```{r}
# change the left part of the equivalence 
  RK1_div_phy = glom
```
You will get: 

For unique OTU's and number of reads classified FIGURE:
*R16S_King, R16S_Phy, R16S_Cla, R16S_Ord, R16S_Fam, R16S_Gen, R16S_Spe
*RK0_Dom, RK0_Phy, RK0_Cla, RK0_Ord, RK0_Fam, RK0_Gen, RK0_Spe
*RK0.1_Dom, RK0.1_Phy, RK0.1_Cla, RK0.1_Ord, RK0.1_Fam, RK0.1_Gen, RK0.1_Spe
*RK1_Dom, RK1_Phy, RK1_Cla, RK1_Ord, RK1_Fam, RK1_Gen, RK1_Spe

For relative abundance with clean phyloseq objects (without R0015 and R0127):
*R16S_Phylum_clean, R16S_Genus_clean
*RK0_Phylum_clean, RK0_Genus_clean
*RK1_Phylum_clean, RK1_Genus_clean

For alpha div analysis:
*R16S_div_gen, RK0_div_gen, RK1_div_gen

For beta div analysis with clean phyloseq objects (without R0015 and R0127):
*R16S_Genus_clean, RK0_Genus_clean, RK0.1_Genus_clean, RK1_Genus_clean

## FIGURE 1
```{r}
  reads_otu <- read_excel("RumenDATA_MASTER.xlsx", sheet = "reads_otu")
  reads_otu <- as.data.frame(reads_otu)
  reads_otu$x <- factor(reads_otu$x, levels=c("OTU/ASV","Domain/Kingdom","Phylum","Class","Order","Family","Genus","Species"))
  reads_otu$group <- factor(reads_otu$group, levels=c("16S","Kraken_0","Kraken_0.1","Kraken_1"))
  
  reads_otu2<-filter(reads_otu,!x=="OTU/ASV", !group=="Kraken_0.1")
levels(reads_otu2$group) <- c("16S-V4", "SMS cs0","SMS cs0.1", "SMS cs1") 

  F1A <- ggplot(reads_otu2, aes(x=x, y=y, group=group, color=group)) +
    geom_line()+
    geom_point()+
    theme_bw()+
    theme(legend.text = element_text(size = 13),axis.text.x = element_text(angle = 45, hjust = 1, size=11))+
    xlab("Taxonomic rank")+
    ylab("Classified reads (%)")+
    ggtitle("Proportion of classified reads resolved to each taxonomic rank")+
    guides(color=guide_legend(title="Method")) 

  reads_otu$x <- factor(reads_otu$x, levels=c("Domain/Kingdom","Phylum","Class","Order","Family","Genus","Species","OTU/ASV"))
  reads_otu3<-filter(reads_otu,!group=="Kraken_0.1")
levels(reads_otu3$group) <- c("16S-V4", "SMS cs0","SMS cs0.1", "SMS cs1") 

    F1B <- ggplot(reads_otu3, aes(x=x, y=unique, group=group, color=group)) +
    geom_line()+
    geom_point()+
    scale_y_continuous(trans='log10')+
    theme_bw()+
    theme(legend.text = element_text(size = 13),axis.text.x = element_text(angle = 45, hjust = 1,size = 11))+
    xlab("Taxonomic rank")+
      ylab("Nº of unique taxa (log10 scale)")+ggtitle("Nº of unique taxa identified at each taxonomic rank")+
    guides(color=guide_legend(title="Method"))  
```
## 3.3 Phyla in common
Agglomerating and formatting ps object before "get_upset" information
```{r}
# Agglomerate the merge phyloseq to phylum level
# FIRST tax_glom to genus level (dropping NAs)
ps_ALL_phy <- tax_glom(ps_ALL, "Phylum", bad_empty = c("NA", NA, "", " ", "\t"))
# changing replicated phylum (same bacteria but different phylum name)
tax_mat<-as.data.frame(tax_table(ps_ALL_phy)) # save as dataframe
# change by the correct names
tax_mat$Phylum[tax_mat$Phylum == "Acidobacteria"] <- "Acidobacteriota"
tax_mat$Phylum[tax_mat$Phylum == "Actinobacteriota"] <- "Actinobacteria"
tax_mat$Phylum[tax_mat$Phylum == "Armatimonadetes"] <- "Armatimonadota"
tax_mat$Phylum[tax_mat$Phylum == "Bacteroidetes"] <- "Bacteroidota"
tax_mat$Phylum[tax_mat$Phylum == "Elusimicrobia"] <- "Elusimicrobiota"
tax_mat$Phylum[tax_mat$Phylum == "Fibrobacteres"] <- "Fibrobacterota"
tax_mat$Phylum[tax_mat$Phylum == "Fusobacteria"] <- "Fusobacteriota"
tax_mat$Phylum[tax_mat$Phylum == "Gemmatimonadetes"] <- "Gemmatimonadota"
tax_mat$Phylum[tax_mat$Phylum == "Nitrospirae"] <- "Nitrospirota"
tax_mat$Phylum[tax_mat$Phylum == "Planctomycetes"] <- "Planctomycetota"
tax_mat$Phylum[tax_mat$Phylum == "Spirochaetes"] <- "Spirochaetota"
tax_mat$Phylum[tax_mat$Phylum == "Synergistetes"] <- "Synergistota"
tax_mat$Phylum[tax_mat$Phylum == "Verrucomicrobia"] <- "Verrucomicrobiota"
# re-built the new ps oject with the corrected taxa table
tax_mat<- read_csv("Gen_tax.csv")
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("otu")
# Transform into matrix tax tables (sample table can be left as data frame)
tax_mat <- as.matrix(tax_mat)
# Transform to phyloseq objects
OTU = otu_table(ps_ALL_phy)
TAX = tax_table(as.matrix(tax_mat))
samples = sample_data(ps_ALL_phy)
ps_ALL_phy2 <- phyloseq(OTU, TAX, samples)
# SECOND tax_glom to genus level (dropping NAs)
ps_ALL_phy2 <- tax_glom(ps_ALL_phy2, "Phylum", bad_empty = c("NA", NA, "", " ", "\t"))
View(tax_table(ps_ALL_phy2))
#Getting the upsetda information for the upset plot
# library("MicrobiotaProcess")
upsetda_phy <- get_upset(ps_ALL_phy2, factorNames="Method") ## ASV
colnames(upsetda_phy) <- c("16S-V4", "SMS cs0","SMS cs1")
```

## SUP FIGURE 2 
```{r}
F2A <- upset(upsetda_phy, sets=c("16S-V4", "SMS cs0","SMS cs1"),text.scale = 2,sets.bar.color = c( "#17b12b", "#5086ff", "#f35e5a"), order.by = "freq", empty.intersections = "on")
```

## 3.4. Genera in common (UpSet plot)
Agglomerating and formatting ps object before "get_upset" information
```{r}
# create ps objects with modified number of samples and OTU/ASV: DO IT FOR 16S, K0, K1
otu_mat<- read_excel("List of common taxa.xlsx", sheet = "OTU_K0")
tax_mat<- read_excel("List of common taxa.xlsx", sheet = "Tax_K0")
samples_df<- read_excel("List of common taxa.xlsx", sheet = "Samples_K0")
# formatting and re-leveling categorical variables 
samples_df$collection_date= factor(samples_df$collection_date, levels=c("9.21.21", "10.18.21", "10.20.21"))
samples_df$collection_day= factor(samples_df$collection_day, levels=c("Pre_weaning", "At_weaning", "Post_weaning"))
samples_df$castration_group= factor(samples_df$castration_group, levels=c("Birth", "Turnout", "Pre_weaning","Weaning"))
samples_df$castration= factor(samples_df$castration, levels=c("Birth", "Turnout", "Pre_weaning","Weaning", "Not_castrated"))
samples_df$castration_status= factor(samples_df$castration_status, levels = c("Castrated", "Not_castrated"))
samples_df$weaning= factor(samples_df$weaning, levels=c("Fence_line", "Truck","Not_weaned"))
samples_df$weaning_group=factor(samples_df$weaning_group, levels=c("Fence_line", "Truck"))
samples_df$weaning_status=factor(samples_df$weaning_status, levels=c("Weaned", "Not_weaned"))
samples_df$cow_ID= as.factor(samples_df$cow_ID)
samples_df$extraction_date= as.factor(samples_df$extraction_date)
samples_df$extraction_run= as.factor(samples_df$extraction_run)
# Define the row names from the otu column in otu_mat and tax_mat
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("otu") 
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("otu")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("samples")
# Transform into matrixes otu and tax tables (sample table can be left as data frame)
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
# Transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
ps.object <- phyloseq(OTU, TAX, samples)
# Change the left part of the equivalence 
ps_K0 = ps.object

## 1.3. Dropping viruses, viroids, eukaryota reads from kraken classification
ps_K0 <- subset_taxa(ps_K0, !Kingdom %in% c("Viruses","Viroids","Eukaryota")) # Raw 11588 OTU, clean 4706 OTU
ps_K1 <- subset_taxa(ps_K1, !Kingdom %in% c("Viruses","Viroids","Eukaryota")) # Raw 2796 OTU, clean 2772 OTU
# REVIEW PREVIOUS STEPS
# Clean phyloseq object
# ps_16S: 7653 taxa
# ps_K0: 4706 taxa
# ps_K1: 2772 taxa
# Subset the outliers from all the phyloseq objects: R0015 and R0127
ps_16S <- prune_samples(!(sample_names(ps_16S) %in% c("R0015", "R0127")), ps_16S)
ps_K0 <- prune_samples(!(sample_names(ps_K0) %in% c("K0_R0015", "K0_R0127")), ps_K0)
ps_K1 <- prune_samples(!(sample_names(ps_K1) %in% c("K1_R0015", "K1_R0127")), ps_K1)

# YOU WILL OBTAIN: ps_16S, ps_K0, ps_K1. THEY ARE EXACTLY THE SAME AS R_16S_Sample_clean, R_K0_Sample_clean AND R_K1_Sample_clean BUT WITH MODIFIED SAMPLE NAMES AND OTU NUMBERS SO THEY DON'T GET COMBINED WHEN PS OBJECT ARE MERGED.
```
merge ps objects to a full single ps object with a new variable "method" that account for the 3 methods
```{r}
ps_16S_K1 <-merge_phyloseq(ps_16S,ps_K1)
ps_ALL <- merge_phyloseq(ps_16S_K1,ps_K0)
# take only samples
ps_ALL <- subset_samples(ps_ALL, sample_type =="Sample")
#ps_ALL = 15131 taxa, 279 samples, 73 variables, 7 taxonomic ranks. PERFECT!
# categorizing our new variable "METHOD" 
sample_data(ps_ALL)$Method <- as.factor(sample_data(ps_ALL)$Method)
# The new phyloseq object will sum up the number of taxa and the number of samples because 
# they have different numbers (the taxa) and different names (the samples) and 
# they will not be recognized as the same thing. 
# Also I have categorized the variable Method in 3 levels: 16S, Kraken_0, Kraken_1
```
Agglomerate and format tables before obtaining information for upset plot (work around duplicate genera because of different taxonomic-tree classification in higher taxonomic ranks)
```{r}
# FIRST tax_glom to genus level (dropping NAs)
ps_ALL_gen <- tax_glom(ps_ALL, "Genus", bad_empty = c("NA", NA, "", " ", "\t"))
# View(phyloseq::tax_table(ps_ALL_gen2))
# WE GOT DUPLICATED GENERA BECAUSE OF TAXONOMIC INFORMATION AT HIGHER TAXONOMIC RANKS

# WORK AROUND: Manually moved genus column to the first column (where kingdom is) and then re-agglomerate
# write the tax_table to fix it manually with excel
write.csv(phyloseq::tax_table(ps_ALL_gen), "Gen_tax.csv")
# re-built the new ps oject
tax_mat<- read_csv("Gen_tax.csv")
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("otu")
# Transform into matrix tax tables (sample table can be left as data frame)
tax_mat <- as.matrix(tax_mat)
# Transform to phyloseq objects
OTU = otu_table(phyloseq::otu_table(ps_ALL_gen))
TAX = phyloseq::tax_table(tax_mat)
samples = sample_data(ps_ALL_gen)
ps.object <- phyloseq(OTU, TAX, samples)
# Change the left part of the equivalence 
ps_ALL_gen2 = ps.object

# SECOND tax_glom to genus level (dropping NAs)
ps_ALL_gen2 <- tax_glom(ps_ALL_gen2, "Genus", bad_empty = c("NA", NA, "", " ", "\t"))
# View(phyloseq::tax_table(ps_ALL_gen))
```
Obtaining information for upset plot and plotting
```{r}
# obtain information for upset plot
# install Microbiota process
#BiocManager::install("MicrobiotaProcess")
# library("MicrobiotaProcess")
upsetda <- get_upset(ps_ALL_gen2, factorNames="Method") ## ASV
View(phyloseq::tax_table(ps_ALL_gen2))
colnames(upsetda) <- c("16S-V4", "SMS cs0","SMS cs1")
```
## FIGURE 3
```{r}
upset_genus <- upset(upsetda, sets=c("16S-V4", "SMS cs0","SMS cs1"), sets.bar.color = c( "#17b12b", "#5086ff", "#f35e5a"),text.scale = 1,
                     order.by = "freq", empty.intersections = "on")

######### CODE THAT FINALLY WORKED AS 04-09-2023 ###########
# 16S=403 genera
# K0=1145 genera
# K1=861 genera

# write the same table for tax and otu tables of 16S, K1, K0
write.csv(phyloseq::tax_table(RK1_Gen), "K1Gen_tax.csv")
write.csv(phyloseq::otu_table(RK1_Gen), "K1Gen_otu.csv")

# create ps objects with modified number of samples and OTU/ASV: DO IT FOR 16S, K0, K1
otu_mat<- read_excel("List of common taxa_2.xlsx", sheet = "OTU_K1Gen")
tax_mat<- read_excel("List of common taxa_2.xlsx", sheet = "Tax_K1Gen")
samples_df<- read_excel("List of common taxa.xlsx", sheet = "Samples_K1")
# Define the row names from the otu column in otu_mat and tax_mat
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("otu") 
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("otu")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("samples")
# Transform into matrixes otu and tax tables (sample table can be left as data frame)
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
# Transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
ps.object <- phyloseq(OTU, TAX, samples)
# Change the left part of the equivalence 
ps_gen_K1 = ps.object

# you will get
ps_gen_16S
ps_gen_K0
ps_gen_K1

#merge ps objects to a full single ps object with a new variable "method" that account for the 3 methods
ps_gen_16S_K1 <-merge_phyloseq(ps_gen_16S,ps_gen_K1)
ps_gen_ALL <- merge_phyloseq(ps_gen_16S_K1,ps_gen_K0)

# categorizing our new variable "METHOD" 
sample_data(ps_gen_ALL)$Method <- as.factor(sample_data(ps_gen_ALL)$Method)
# The new phyloseq object will sum up the number of taxa and the number of samples because 
# they have different numbers (the taxa) and different names (the samples) and 
# they will not be recognized as the same thing. 
# Also I have categorized the variable Method in 3 levels: 16S, Kraken_0, Kraken_1
ps_gen_ALL <- tax_glom(ps_gen_ALL, "Genus", bad_empty = c("NA", NA, "", " ", "\t"))

upsetda2 <- get_upset(ps_gen_ALL, factorNames="Method") ## ASV
View(phyloseq::tax_table(ps_gen_ALL))
colnames(upsetda2) <- c("16S-V4", "SMS cs0","SMS cs1")

F2B <- upset(upsetda2, sets=c("16S-V4", "SMS cs0","SMS cs1"), sets.bar.color = c( "#17b12b", "#5086ff", "#f35e5a"),
                     order.by = "freq", empty.intersections = "on")


```

# 4. Relative abundance
We will analyze to Phylum and Genus level 
## 4.1. Phylum level
```{r}
# Change the right part of the equivalence 
  ps = RK1_Phylum_clean
# We will use R16S_Phylum_clean, RK0_Phylum_clean, RK1_Phylum_clean
```
Obtain the relative_long for ggplot:
```{r}
# 1. Converting to Percentage and transforming Phylum tax_glom object to long-format table
  relative <- transform_sample_counts(ps, function(x) (x / sum(x))*100 )
  relative_long <- psmelt(relative)
# 2. Creating a new column (variable) for the  long-format table, called "mean_relative_abund" 
  relative_long <- relative_long %>%
  group_by(Phylum) %>%
  mutate(mean_relative_abund = mean(Abundance))
# 3. Formatting the "Phylum" column variable as names and the "mean_relative_abund" column variable as continuous numbers
  relative_long$Phylum <- as.character(relative_long$Phylum)
  relative_long$mean_relative_abund <- as.numeric(relative_long$mean_relative_abund)
# 4. Giving the name "Others (< 1%)" to all the OTUs in the "Genus" column that have less than 1 in the "Abundace" column in the long-format table 
  relative_long$Phylum[relative_long$Abundance < 1] <- "Others (< 1%)" 
```
Exploratory figure: Relative abundance per collection day
```{r}
relative_long %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day+weaning, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20")
```

## SUP FIGURE 1. Classification at phylum level
```{r}
#For 16S:
#relative_long$Phylum = factor(relative_long$Phylum, levels=c("NA",	"Bdellovibrionota",	"Patescibacteria",	"Planctomycetota",	"Spirochaetota",	"Thermoplasmatota",	"Desulfobacterota",	"Verrucomicrobiota", "Fusobacteriota",	"Cyanobacteria",	"Proteobacteria",	"Actinobacteriota",	"Fibrobacterota",	"Euryarchaeota",	"Firmicutes",	"Bacteroidota",	"Others (< 1%)")) # Change order in plot HERE
# For K0
#relative_long$Phylum <- relative_long$Phylum %>% replace_na('NA')
#relative_long$Phylum = factor(relative_long$Phylum, levels=c("NA","Cyanobacteria",	"Proteobacteria",	"Actinobacteria",	"Fibrobacteres",	"Euryarchaeota",	"Firmicutes",	"Bacteroidetes",	"Others (< 1%)"))
#K0P = c("#3f7aab","#86BCB6",	"#E15759",	"#FF9D9A",	"#79706E",	"#BAB0AC",	"#D37295",	"#FABFD2",	"#B07AA1")
# For K1
#relative_long$Phylum <- relative_long$Phylum %>% replace_na('NA')
#relative_long$Phylum = factor(relative_long$Phylum, levels=c("NA","Fusobacteria",	"Synergistetes",	"Proteobacteria",	"Actinobacteria",	"Fibrobacteres",	"Euryarchaeota",	"Firmicutes",	"Bacteroidetes",	"Others (< 1%)"))
#K1P = c("#3f7aab","#499894",	"#86BCB6",	"#E15759",	"#FF9D9A",	"#79706E",	"#BAB0AC",	"#D37295",	"#FABFD2",	"#B07AA1")

SF2A <- relative_long %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity",width = 1,color="black") +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day+weaning, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20")
SF2B <- relative_long %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) +
  #ggplot(aes(x = Sample, y = Abundance, fill = fct_reorder(Genus, Abundance))) +
  scale_fill_manual(values = K0P)+
  geom_bar(stat = "identity",width = 1,color="black") +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day+weaning, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")
  #ggthemes::scale_fill_tableau("Tableau 20")
SF2C <- relative_long %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) +
  #ggplot(aes(x = Sample, y = Abundance, fill = fct_reorder(Genus, Abundance))) +
  scale_fill_manual(values = K1P)+
  geom_bar(stat = "identity",width = 1,color="black") +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day+weaning, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")
  #ggthemes::scale_fill_tableau("Tableau 20")
```
Obtain phylum-level relative abundance for Bacteria
```{r}
# Overall Phyla abundance in bacteria
  bac_Phy<-table1(~ Abundance | Phylum, data=relative_long)
```


## 4.2. Genus level
```{r}
# Change the right part of the equivalence 
  ps = RK1_Genus_clean
# We will use R16S_Genus_clean, RK0_Genus_clean, RK1_Genus_clean
```
Obtain the relative_long for ggplot:
```{r}
# 1. Converting to Percentage and transforming Genus tax_glom object to long-format table
  relative2 <- transform_sample_counts(ps, function(x) (x / sum(x))*100 )
  relative_long2 <- psmelt(relative2)
# 2. Creating a new column (variable) for the  long-format table, called "mean_relative_abund" 
  relative_long2 <- relative_long2 %>%
  group_by(Genus) %>%
  mutate(mean_relative_abund = mean(Abundance))
# 3. Formatting the "Genus" column variable as names and the "mean_relative_abund" column variable as continuous numbers
  relative_long2$Genus <- as.character(relative_long2$Genus)
  relative_long2$mean_relative_abund <- as.numeric(relative_long2$mean_relative_abund)
# 4. Giving the name "Others (< 1%)" to all the OTUs in the "Genus" column that have less than 1 in the "Abundace" column in the long-format table 
  relative_long2$Genus[relative_long2$Abundance < 3] <- "Others (< 3%)" 
```

## FIGURE 2. Genus level
```{r}
# BACTERIA
#Sixteen = c("#4E79A7","#A0CBE8","#F28E2B","#FFBE7D","#59A14F","#8CD17D","#B6992D","#F1CE63","#499894","#86BCB6","#E15759","#FF9D9A","#79706E","#BAB0AC","#D37295","#FABFD2","#B07AA1","#D4A6C8","#D7B5A6")
#Sixteen2 = c("#D7B5A6", "#9D7660", "#D4A6C8", "#B07AA1", "#FABFD2","#D37295","#BAB0AC","#79706E","#FF9D9A","#E15759","#86BCB6", "#499894","#F1CE63","#B6992D","#8CD17D","#59A14F","#FFBE7D","#F28E2B","#A0CBE8","#4E79A7")
K0 = c("#3f7aab","#F7F0A8",	"#BAB0AC",	"#D37295",	"#FABFD2",	"#B07AA1",	"#D4A6C8",	"#9D7660")
#K0 = c("#4E79A7","#a2ceaa","#f47942", "#B07AA1", "#59A14F","#E15759","#D7B5A6")
#K02 = c("#D7B5A6", "#a2ceaa", "#f47942","#FFBE7D", "#FABFD2","#86BCB6","#A0CBE8")
K1 = c("#3f7aab","#B7F1A5",	"#D5F7C4",	"#FFF3F0",	"#DCFFFB",	"#C1E7E3",	"#F7F0A8",	"#BAB0AC",	"#FABFD2",	"#B07AA1",	"#D4A6C8",	"#9D7660")
#K1 = c("#4E79A7","#79706E","#B07AA1","#7e756d","#b9aa97","#d7ce9f","#b66353","#fbb04e","#f47942","#bfbb60", "#59A14F", "#E15759", "#D7B5A6")
#K12 = c("#D7B5A6","#79706E","#FFBE7D", "#7e756d", "#b9aa97", "#d7ce9f", "#b66353", "#fbb04e", "#f47942", "#bfbb60", "#FABFD2", "#86BCB6", "#A0CBE8")

# For 16S Genus:
#relative_long2$Genus = factor(relative_long2$Genus, levels=c("NA",	"UCG-002",	"Olsenella",	"Lachnospiraceae NK3A20 group",	"Succiniclasticum",	"Anaeroplasma",	"Rikenellaceae RC9 gut group",	"Ruminococcus",	"Selenomonas",	"Pseudobutyrivibrio",	"Treponema",	"NK4A214 group",	"Clostridium sensu stricto 1",	"Christensenellaceae R-7 group",	"Methanobrevibacter",	"Fibrobacter",	"Butyrivibrio",	"Prevotella",	"Others (< 3%)")) # Change order in plot HERE
# For K0:
#relative_long2$Genus <- relative_long2$Genus %>% replace_na('NA')
#relative_long2$Genus = factor(relative_long2$Genus, levels=c("NA","Bacteroides",	"Eubacterium",	"Methanobrevibacter",	"Fibrobacter",	"Butyrivibrio",	"Prevotella",	"Others (< 3%)"))
# For K1:
relative_long2$Genus <- relative_long2$Genus %>% replace_na('NA')
relative_long2$Genus = factor(relative_long2$Genus, levels=c("NA","Chryseobacterium",	"Listeria",	"Xanthomonas",	"Selenomonas",	"Trueperella",	"Fusobacterium",	"Eubacterium",	"Fibrobacter",	"Butyrivibrio",	"Prevotella",	"Others (< 3%)"))

F3A1 <- relative_long2 %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
  #ggplot(aes(x = Sample, y = Abundance, fill = fct_reorder(Genus, Abundance))) +
  #scale_fill_manual(values = Sixteen)+
  geom_bar(stat = "identity",width = 1,color="black") +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day+weaning, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20")+
  guides(fill=guide_legend(title="Genus"))
F3A2 <- relative_long2 %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
  #ggplot(aes(x = Sample, y = Abundance, fill = fct_reorder(Genus, Abundance))) +
  scale_fill_manual(values = K0)+
  geom_bar(stat = "identity",width = 1,color="black") +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day+weaning, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")
  #ggthemes::scale_fill_tableau("Tableau 20")+
  #guides(fill=guide_legend(title="Genus"))
F3A3 <- relative_long2 %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
  #ggplot(aes(x = Sample, y = Abundance, fill = fct_reorder(Genus, Abundance))) +
  scale_fill_manual(values = K1)+
  geom_bar(stat = "identity",width = 1,color="black") +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day+weaning, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")
  #ggthemes::scale_fill_tableau("Tableau 20")+
  #guides(fill=guide_legend(title="Genus"))
```
Obtain tables for relative abundance
```{r}
 # Overall Genera abundance in bacteria
  arc_Gen <- table1(~ Abundance | Genus, data=relative_long2)
  bac_Gen <- table1(~ Abundance | Genus, data=relative_long2)
```

## 4.3. Genera in common
Formatting dataframes to produce scatterplot + marginal histogram for genera in common
```{r}
# subsetting from ps object that have all in raw counts
df_16S<- subset_samples(ps_ALL_gen2, Method == "16S")
df_K0<- subset_samples(ps_ALL_gen2, Method == "Kraken_0")
df_K1<- subset_samples(ps_ALL_gen2, Method == "Kraken_1")

# Summing up the number of reads per OTU
#df_16S
df_16S <-as.data.frame(taxa_sums(df_16S))
colnames(df_16S)[1] ="R16S" 
df_16S$R16S_only <- ifelse(df_16S$R16S > 0, "1", "0")
df_16S$R16S_log2 <- log2(df_16S$R16S)
df_16S$R16S_log10 <- log10(df_16S$R16S)

#df_K0
df_K0 <-as.data.frame(taxa_sums(df_K0))
colnames(df_K0)[1] ="RK0" 
df_K0$RK0_only <- ifelse(df_K0$RK0 > 0, "1", "0")
df_K0$RK0_log2 <- log2(df_K0$RK0)
df_K0$RK0_log10 <- log10(df_K0$RK0)

#df_K1
df_K1 <-as.data.frame(taxa_sums(df_K1))
colnames(df_K1)[1] ="RK1" 
df_K1$RK1_only <- ifelse(df_K1$RK1 > 0, "1", "0")
df_K1$RK1_log2 <- log2(df_K1$RK1)
df_K1$RK1_log10 <- log10(df_K1$RK1)

#creating dataframe for the figures
df_K0_16S <-cbind(df_K0, df_16S)
df_K1_16S <-cbind(df_K1, df_16S)
df_K0_K1 <-cbind(df_K0, df_K1)
#df_all <- cbind(df_K0, df_K1, df_16S)

# Final dataframe for the figures
df_K0_16S2 <-df_K0_16S %>%
  mutate(Genera_identified_by = case_when(RK0_only == 1 & R16S_only == 1 ~ "16S-V4 and SMS cs0",
                                    RK0_only == 1 & R16S_only == 0 ~ "Only SMS cs0",
                                    RK0_only == 0 & R16S_only == 1 ~ "Only 16S-V4")) %>% drop_na()
df_K1_16S2 <-df_K1_16S %>%
  mutate(Genera_identified_by = case_when(RK1_only == 1 & R16S_only == 1 ~ "16S-V4 and SMS cs1",
                                          RK1_only == 1 & R16S_only == 0 ~ "Only SMS cs1",
                                          RK1_only == 0 & R16S_only == 1 ~ "Only 16S-V4")) %>% drop_na()
df_K0_K12 <-df_K0_K1 %>%
  mutate(Genera_identified_by = case_when(RK0_only == 1 & RK1_only == 1 ~ "SMS cs0 and cs1",
                                          RK0_only == 1 & RK1_only == 0 ~ "Only SMS cs0",
                                          RK0_only == 0 & RK1_only == 1 ~ "Only SMS cs1")) %>% drop_na()
```

## FIGURE X. Genus distribution
Plot it
```{r}

 # Kraken cs0 vs 16S
 p1 <- df_K0_16S2 %>%
   ggplot(aes(x=R16S_log2, RK0_log2, color=Genera_identified_by))+
   geom_point() +
   theme(legend.position="none")+
   theme_bw()+
   labs(title = "Genera distribution (16S-V4 vs SMS cs0)", x = "Genera by 16S-V4 (log2 scale)", y = "Genera by SMS cs0 (log2 scale)")

 F4A<- ggMarginal(p1, type="histogram", groupColour = FALSE, groupFill = TRUE, size = 3)
 #cor.testp1 <-cor.test(df_K0_16S2$R16S, df_K0_16S2$RK0, method = "pearson", conf.level = 0.95)

  # Kraken cs1 vs 16S
 p2 <- df_K1_16S2 %>%
   ggplot(aes(x=R16S_log2, RK1_log2, color=Genera_identified_by))+
   geom_point() +
   theme(legend.position="none")+
   theme_bw()+
   labs(title = "Genera distribution (16S-V4 vs SMS cs1)", x = "Genera by 16S-V4 (log2 scale)", y = "Genera by SMS cs1 (log2 scale)")
 
  F4B<- ggMarginal(p2, type="histogram", groupColour = FALSE, groupFill = TRUE, size = 3) 
# cor.testp2 <-cor.test(df_K1_16S2$R16S, df_K1_16S2$RK1, method = "pearson", conf.level = 0.95)

 # Kraken cs1 vs cs0
 df_K0_K12$Genera_identified_by = factor(df_K0_K12$Genera_identified_by, levels=c("SMS cs0 and cs1", "Only SMS cs0", "Only SMS cs1"))
 p3 <- df_K0_K12 %>%
   ggplot(aes(x=RK1_log2, RK0_log2, color=Genera_identified_by))+
   geom_point() +
   theme(legend.position="none")+
   theme_bw() +
   labs(title = "Genera distribution (SMS cs1 vs cs0)", x = "Genera by SMS cs1 (log2 scale)", y = "Genera by SMS cs0 (log2 scale)")

  F4C<- ggMarginal(p3, type="histogram", groupColour = FALSE, groupFill = TRUE, size = 3)
 #cor.testp3 <-cor.test(df_K0_K12$RK1,df_K0_K12$RK0, method = "pearson", conf.level = 0.95)
```

# 5. Genera Abundance correlation
We only will analyze the genera in common per pairwise
```{r}
# formatting
phyloseq::otu_table(ps_ALL_gen2)
upsetda2 <- tibble::rownames_to_column(upsetda, var = "otu")
colnames(upsetda2)[2] ="R16S" # because name 16S was troublesome 
upsetda2 <- upsetda2 %>% filter(R16S == 1 & Kraken_0 == 1 & Kraken_1 == 1) # to keep only common otus
# subsetting only common otus
ps_ALL_gen_common <- prune_taxa(upsetda2$otu, ps_ALL_gen2)
df <- as.data.frame(otu_table(ps_ALL_gen_common))

# correlation plot
# 16S(x) vs Kraken_0(y)
cor1 <- ggplot(df, aes(R0025, K0_R0025))+ # Pre-weaning
  geom_point() + 
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  labs(x = "16S(log10 scale)", y = "Kraken_0(log10 scale)") + 
  geom_smooth(method="lm")
cor2 <- ggplot(df, aes(R0001, K0_R0001))+ # At-weaning
  geom_point() + 
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  labs(x = "16S(log10 scale)", y = "Kraken_0(log10 scale)") + 
  geom_smooth(method="lm")
cor3 <- ggplot(df, aes(R0003, K0_R0003))+ # Post-weaning TRUCK WEANING
  geom_point() + 
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  labs(x = "16S(log10 scale)", y = "Kraken_0(log10 scale)") + 
  geom_smooth(method="lm")
cor4 <- ggplot(df, aes(R0020, K0_R0020))+ # Post-weaning FENCE-LINE WEANING
  geom_point() + 
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  labs(x = "16S(log10 scale)", y = "Kraken_0(log10 scale)") + 
  geom_smooth(method="lm")

# 16S(x) vs Kraken_1(y)
cor5 <- ggplot(df, aes(R0025, K1_R0025))+
  geom_point() + 
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  labs(x = "16S(log10 scale)", y = "Kraken_1(log10 scale)") + 
  geom_smooth(method="lm")
cor6 <- ggplot(df, aes(R0001, K1_R0001))+
  geom_point() + 
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  labs(x = "16S(log10 scale)", y = "Kraken_1(log10 scale)") + 
  geom_smooth(method="lm")
cor7 <- ggplot(df, aes(R0003, K1_R0003))+
  geom_point() + 
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  labs(x = "16S(log10 scale)", y = "Kraken_1(log10 scale)") + 
  geom_smooth(method="lm")
cor8 <- ggplot(df, aes(R0020, K1_R0020))+
  geom_point() + 
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  labs(x = "16S(log10 scale)", y = "Kraken_1(log10 scale)") + 
  geom_smooth(method="lm")
# Kraken_1(x) vs Kraken_0(y)
cor9 <- ggplot(df, aes(K1_R0025, K0_R0025))+
  geom_point() + 
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  labs(x = "Kraken_1(log10 scale)", y = "Kraken_0(log10 scale)") + 
  geom_smooth(method="lm")
cor10 <- ggplot(df, aes(K1_R0001, K0_R0001))+
  geom_point() + 
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  labs(x = "Kraken_1(log10 scale)", y = "Kraken_0(log10 scale)") + 
  geom_smooth(method="lm")
cor11 <- ggplot(df, aes(K1_R0003, K0_R0003))+
  geom_point() + 
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  labs(x = "Kraken_1(log10 scale)", y = "Kraken_0(log10 scale)") + 
  geom_smooth(method="lm")
cor12 <- ggplot(df, aes(K1_R0020, K0_R0020))+
  geom_point() + 
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  labs(x = "Kraken_1(log10 scale)", y = "Kraken_0(log10 scale)") + 
  geom_smooth(method="lm")

# Pearson correlation and p-value
cor.test1 <-cor.test(df$R0025, df$K0_R0025, method = "pearson", conf.level = 0.95)
cor.test2 <-cor.test(df$R0001, df$K0_R0001, method = "pearson", conf.level = 0.95)
cor.test3 <-cor.test(df$R0003, df$K0_R0003, method = "pearson", conf.level = 0.95)
cor.test4 <-cor.test(df$R0020, df$K0_R0020, method = "pearson", conf.level = 0.95)
cor.test5 <-cor.test(df$R0025, df$K1_R0025, method = "pearson", conf.level = 0.95)
cor.test6 <-cor.test(df$R0001, df$K1_R0001, method = "pearson", conf.level = 0.95)
cor.test7 <-cor.test(df$R0003, df$K1_R0003, method = "pearson", conf.level = 0.95)
cor.test8 <-cor.test(df$R0020, df$K1_R0020, method = "pearson", conf.level = 0.95)
cor.test9 <-cor.test(df$K1_R0025, df$K0_R0025,  method = "pearson", conf.level = 0.95)
cor.test10 <-cor.test(df$K1_R0001, df$K0_R0001, method = "pearson", conf.level = 0.95)
cor.test11 <-cor.test(df$K1_R0003, df$K0_R0003, method = "pearson", conf.level = 0.95)
cor.test12 <-cor.test(df$K1_R0020, df$K0_R0020, method = "pearson", conf.level = 0.95)
```

Make the cor.test in a pairwise way and save the output
Define your "df". Use R16S_K0, R16S_K1, RK1_K0
```{r}
# 16S vs K0: R16S_K0
#df <- R16S_K0
# 16S vs K1: R16S_K1
#df <- R16S_K1
# K1 vs K0: RK1_K0
#df <- RK1_K0
```

Run loop
```{r}
results_df <- data.frame()
# loop through column pairs
for(i in 1:(ncol(df)-1)){
  for(j in (i+1):ncol(df)){
    if(substr(colnames(df)[i], nchar(colnames(df)[i])-2, nchar(colnames(df)[i])) == substr(colnames(df)[j], nchar(colnames(df)[j])-2, nchar(colnames(df)[j]))){
      # perform cor.test
      cor_test <- cor.test(df[,i], df[,j], method="pearson", conf.level=0.95)
      # extract relevant information from cor.test output
      result <- data.frame(Column1=colnames(df)[i], Column2=colnames(df)[j], Correlation=cor_test$estimate, P_value=cor_test$p.value)
      # add result to results_df
      results_df <- rbind(results_df, result)
    }
  }
}

# print results dataframe
print(results_df)
```

Recover your output
```{r}
# change the left part of the equation
# 16S vs K0: R16S_K0
#R2_R16S_K0 <- results_df
# 16S vs K1: R16S_K1
#R2_R16S_K1 <- results_df
# K1 vs K0: RK1_K0
# R2_RK1_K0 <- results_df
```

## FIGURE 4

```{r}
# adding new column to each dataframe giving their classification
R2_R16S_K0$comparison <- paste("16S vs K0")
R2_R16S_K1$comparison <- paste("16S vs K1")
R2_RK1_K0$comparison <- paste("K1 vs K0")
# merging all dataframes
merged_R2 <- rbind(R2_R16S_K0, R2_R16S_K1, R2_RK1_K0)
merged_R2$comparison<-as.factor(merged_R2$comparison) # factoring 
# testing ANOVA
ANOVA<-aov(merged_R2$Correlation~merged_R2$comparison)
summary(ANOVA)
pairwise.t.test(merged_R2$Correlation,merged_R2$comparison, p.adj = "bonferroni")
# testing kruskal wallis
kruskal.test(merged_R2$Correlation~merged_R2$comparison)
pairwise.wilcox.test(merged_R2$Correlation,merged_R2$comparison, p.adjust.method = "BH")

# plotting
levels(merged_R2$comparison) <- c("16S-V4 vs SMS cs0", "16S-V4 vs SMS cs1","SMS cs0 vs SMS cs1")

F5<-ggplot(merged_R2, aes(x = comparison, y= Correlation, fill= factor(comparison)))+
  geom_point(position = position_jitter(width = 0.2), alpha = 0.7)+
  geom_boxplot(position=position_dodge(0.8),lwd=1, alpha= 0.7)+
  theme_classic()+ 
  theme(axis.text.x = element_blank())+
  labs(title = "R2 values for genera abundance correlation", x = "Comparison type", y = "R2")+  
  scale_fill_discrete(name = "Comparison type") 

```


# 6. Alpha diversity
We will estimate the alpha diversity indices and plot them. Finally we will modelate properly our indices as outcomes and the predictors of our study.
## 6.1. Estimate alpha diversity
Define your raw phyloseq object to be subsetted and the Domains you want to subset. Also define your level for taxa_glom (Usually genus)
```{r}
# Change the right part of the equivalence 
  ps = RK1_div_gen
```
Estimate the richness, Simpson, Shannon and Evenness indices
```{r}
  div1= estimate_richness(ps, measures=c("Observed", "InvSimpson", "Shannon"))
  div2= evenness(ps, 'pielou')
# estimate diversity metrics at the phylum level
  div_16S= estimate_richness(R16S_div_phy, measures=c("Observed", "InvSimpson", "Shannon"))
  div_16S= cbind(div_16S,(evenness(R16S_div_phy, 'pielou')))
  names(div_16S) <- paste(names(div_16S), "16S", sep = "_")

  div_K0= estimate_richness(RK0_div_phy, measures=c("Observed", "InvSimpson", "Shannon"))
  div_K0= cbind(div_K0,(evenness(RK0_div_phy, 'pielou')))
  names(div_K0) <- paste(names(div_K0), "K0", sep = "_")

  div_K1= estimate_richness(RK1_div_phy, measures=c("Observed", "InvSimpson", "Shannon"))
  div_K1= cbind(div_K1, (evenness(RK1_div_phy, 'pielou')))
  names(div_K1) <- paste(names(div_K1), "K1", sep = "_")
# get phylum-level diversity indices table for all of them
  div_phy = cbind(div_16S,div_K0,div_K1)
```
Save the results in a table
```{r}
# Change the name of the csv file. div1 is for richness, simpson, shannon) and div2 is for evenness
  write.csv(div1, "RK1_div1_genus.csv")
  write.csv(div2, "RK1_div2_genus.csv")
```

## 6.2. Plotting alpha diversity
Re-level the categorical variables
```{r}
samples_2$castration_group=factor(samples_2$castration_group, levels=c("Birth", "Turnout", "Pre_weaning", "Weaning"))
samples_2$castration=factor(samples_2$castration, levels=c("Birth", "Turnout", "Pre_weaning", "Weaning", "Not_castrated"))
samples_2$weaning_group=factor(samples_2$weaning_group, levels=c("Fence_line", "Truck"))
samples_2$weaning=factor(samples_2$weaning, levels=c("Fence_line", "Truck", "Not_weaned"))
samples_2$collection_day=factor(samples_2$collection_day, levels=c("Pre_weaning", "At_weaning", "Post_weaning"))
```

Define the variables for the plot:
y_axis can be   = Kra0.1_Obs_Genus, Kra0.1_Shan_Genus, Kra0.1_Simp_Genus, Kra0.1_Even_Genus
y_label can be  = Richness / Shannon Index / Inv-Simpson Index / Pielou's evenness index

FOR colorby = collection_day 
  color = c("brown", "#D98326", "#E8B600")
FOR colorby = weaning_group or weaning
  color = c("#967bba", "#505359", "#c9c8c8")
FOR colorby = castration_group / castration
  color = c("#f54545", "#85a9ca", "#e0b54e", "#73a054", "#c9c8c8")
```{r}
# Change the right part of the formula
  x_axis = "collection_day"
  x_label = "Collection day"
  y_axis = "Shan_K_1_gen"
  y_label = "Shannon Index"
  colorby = "weaning_group"
  color = c("#967bba", "#505359", "#c9c8c8")
```
Exploratory figure. Note that you will need to change the dataframe "samples_2" according to your needs
```{r}
# plotting richness, shannon and Pielou by WEANING and Collection day
  ggplot(samples_2_clean,
         aes_string(x_axis,y_axis, color=colorby))+
  geom_boxplot(position=position_dodge(0.8),lwd=1)+
  scale_fill_manual(values = color)+
  scale_color_manual(values = color)+
  theme_classic()+ 
  #ylim(0,1)+
  geom_jitter(position=position_dodge(0.8),size=2)+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  #theme(axis.text.x = element_blank())+
  ylab(y_label)+
  xlab(x_label)
```
Filtering out "outliers" from 16S dataset
```{r}
# I will filter out a sample from pre-weaning (R0015) and post_weaning (R0127) based on what I saw in Shan_16S_gen results
  #samples_2_clean <- samples_2[!(samples_2_clean$samples %in% c("R0015","R0127")),]
```

## FIGURE 6A: Alpha div
```{r}
F3A3 <-   ggplot(samples_2_clean,
         aes_string(x_axis,y_axis, color=colorby))+
  geom_boxplot(position=position_dodge(0.8),lwd=1)+
  scale_fill_manual(values = color)+
  scale_color_manual(values = color)+
  theme_classic()+ 
  #ylim(0,5.5)+
  scale_y_continuous(limits = c(0,5.5), breaks = seq(0,5.5, by = 0.5))+
  #scale_y_continuous(n.breaks=2)+
  #scale_y_continuous(breaks=c(0,3,6))+
  geom_jitter(position=position_dodge(0.8),size=1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab(y_label)+
  xlab(x_label)
```

## 6.3. Statistics 
We used a linear mixed-effect model to fit random (cow_ID) and fixed effects (castration / castration_group, weaning / weaning_group, age_days, weight_ADGUpToDate)
* Our outcome (response variable) will be genus-level alpha indices: richness and shannon 
* Our predictor (explanatory variable) will be the interventions and time: castration, weaning and collection_days
```{r}
# Define data frames
#md= samples_2
md_clean=samples_2_clean
md1= samples_D1
md2= samples_D2
md3= samples_D3
```
Run model for WEANING
```{r}
  weaning_16S = lmer(Shan_16S_gen ~ weaning_group*collection_day +castration_group + (1|cow_ID), data=md_clean, REML = T)
  summary(weaning_16S)
  #confint(weaning_16S)
  emmeans(weaning_16S,pairwise~weaning_group|collection_day)
  Anova(weaning_16S, type = "III")
  #Anova(weaning_16S, type = "III", test.statistic = "F")
  
  weaning_K0 = lmer(Shan_K_0_gen ~ weaning_group*collection_day +castration_group+ (1|cow_ID), data=md_clean, REML = T)
  summary(weaning_K0)
  #confint(weaning_K0)
  emmeans(weaning_K0,pairwise~weaning_group|collection_day)
  Anova(weaning_K0, type = "III")
  #Anova(weaning_K0, type = "III", test.statistic = "F")

  weaning_K1 = lmer(Shan_K_1_gen ~ weaning_group*collection_day +castration_group + (1|cow_ID), data=md_clean, REML = T)
  summary(weaning_K1)
  #confint(weaning_K1)
  emmeans(weaning_K1,pairwise~weaning_group|collection_day)
  Anova(weaning_K1, type = "III")
  #Anova(weaning_K1, type = "III", test.statistic = "F")
  
  # weaning_obs = lmer(Kra0.1_Obs_Genus ~ weaning_group*collection_day +castration_group + age_days + weight_ADGUpToDate+ (1|cow_ID), data=md, REML = T)
  # summary(weaning_obs)
  # confint(weaning_obs)
  # emmeans(weaning_obs,pairwise~weaning_group|collection_day)
  # Anova(weaning_obs, type = "III")
  # Anova(weaning_obs, type = "III", test.statistic = "F")
```
Run model for CASTRATION
```{r}
# CASTRATION
  castration_shan = lmer(Kra0.1_Shan_Genus~castration_group*collection_day+weaning_group + (1|cow_ID), data=md, REML = T)
  # anova(castration_shan1, castration_shan2)
  summary(castration_shan)
  emmeans(castration_shan,pairwise~castration_group|collection_day)
  Anova(castration_shan, type = "III")
  Anova(castration_shan, type = "III", test.statistic = "F")
  
  castration_obs = lmer(Kra0.1_Obs_Genus~castration_group*collection_day+weaning_group + (1|cow_ID), data=md, REML = T)
  summary(castration_obs)
  emmeans(castration_obs,pairwise~castration_group|collection_day)
  Anova(castration_obs, type = "III")
  Anova(castration_obs, type = "III", test.statistic = "F")
  
# CASTRATION PER DAYS
  castration_shan1 = lm(Kra0.1_Shan_Genus~castration + age_days + weight_ADGUpToDate, data=md1)
  summary(castration_shan1)
  emmeans(castration_shan1,pairwise~castration) 
  Anova(castration_shan1, type = "III")
 
  castration_shan2 = lm(Kra0.1_Shan_Genus~castration + age_days + weight_ADGUpToDate, data=md2)
  summary(castration_shan2)
  emmeans(castration_shan2,pairwise~castration) 
  Anova(castration_shan2, type = "III")
  
  castration_shan3 = lm(Kra0.1_Shan_Genus~castration*weaning_group + age_days + weight_ADGUpToDate, data=md3)
  summary(castration_shan3)
  emmeans(castration_shan3,pairwise~castration) 
  Anova(castration_shan3, type = "III")
  
# WEANING PER DAY
    weaning_shan1 = lm(Kra0.1_Shan_Genus~weaning_group + age_days + weight_ADGUpToDate, data=md1)
    weaning_shan2 = lm(Kra0.1_Shan_Genus~weaning_group + age_days + weight_ADGUpToDate, data=md2)
    weaning_shan3 = lm(Kra0.1_Shan_Genus~weaning_group + age_days + weight_ADGUpToDate, data=md3)
    summary(weaning_shan3)
```

## 6.4. Correlation alpha div
At the Genus level
```{r}
# 16S vs SMS cs0
ggplot(samples_2_clean, aes(Obs_16S_gen, Obs_K_0_gen))+
  geom_point() + 
  #scale_y_continuous(trans='log10')+
  #scale_x_continuous(trans='log10')+
  labs(x = "SMS cs0 (Shannon's Index)", y = "16S-V4 (Shannon's Index)") + 
  geom_smooth(method="lm")

cor.rich1 <-cor.test(samples_2_clean$Obs_16S_gen, samples_2_clean$Obs_K_0_gen, method = "pearson", conf.level = 0.95)
cor.even1 <-cor.test(samples_2_clean$Even_16S_gen, samples_2_clean$Even_K_0_gen, method = "pearson", conf.level = 0.95)
cor.shan1 <-cor.test(samples_2_clean$Shan_16S_gen, samples_2_clean$Shan_K_0_gen, method = "pearson", conf.level = 0.95)

# 16S vs SMS cs1
cor.rich2 <-cor.test(samples_2_clean$Obs_16S_gen, samples_2_clean$Obs_K_1_gen, method = "pearson", conf.level = 0.95)
cor.even2 <-cor.test(samples_2_clean$Even_16S_gen, samples_2_clean$Even_K_1_gen, method = "pearson", conf.level = 0.95)
cor.shan2 <-cor.test(samples_2_clean$Shan_16S_gen, samples_2_clean$Shan_K_1_gen, method = "pearson", conf.level = 0.95)

# SMS cs0 vs SMS cs1
cor.rich3 <-cor.test(samples_2_clean$Obs_K_0_gen, samples_2_clean$Obs_K_1_gen, method = "pearson", conf.level = 0.95)
cor.even3 <-cor.test(samples_2_clean$Even_K_0_gen, samples_2_clean$Even_K_1_gen, method = "pearson", conf.level = 0.95)
cor.shan3 <-cor.test(samples_2_clean$Shan_K_0_gen, samples_2_clean$Shan_K_1_gen, method = "pearson", conf.level = 0.95)
```

At the phylum level
```{r}
#dataframe with diversity metrics for phylum level
div_phy
# 16S vs SMS cs0
cor.rich4 <-cor.test(div_phy$Observed_16S, div_phy$Observed_K0, method = "pearson", conf.level = 0.95)
cor.even4 <-cor.test(div_phy$Shannon_16S, div_phy$Shannon_K0, method = "pearson", conf.level = 0.95)
cor.shan4 <-cor.test(div_phy$pielou_16S, div_phy$pielou_K0, method = "pearson", conf.level = 0.95)

# 16S vs SMS cs1
cor.rich5 <-cor.test(div_phy$Observed_16S, div_phy$Observed_K1, method = "pearson", conf.level = 0.95)
cor.even5 <-cor.test(div_phy$Shannon_16S, div_phy$Shannon_K1, method = "pearson", conf.level = 0.95)
cor.shan5 <-cor.test(div_phy$pielou_16S, div_phy$pielou_K1, method = "pearson", conf.level = 0.95)

# SMS cs0 vs SMS cs1
cor.rich6 <-cor.test(div_phy$Observed_K0, div_phy$Observed_K1, method = "pearson", conf.level = 0.95)
cor.even6 <-cor.test(div_phy$Shannon_K0, div_phy$Shannon_K1, method = "pearson", conf.level = 0.95)
cor.shan6 <-cor.test(div_phy$pielou_K0, div_phy$pielou_K1, method = "pearson", conf.level = 0.95)

```

# 7. Beta diversity
We will estimate the *ordination* for the plot and the *distance matrix* for the statistics
## 7.1. Estimate beta diversity

Define your CSS-normalized genus-agglomerate phyloseq object from which we will get the ordination.
```{r}
# Change the right part of the equivalence 
  ps = RK0_Genus_clean
  set.seed(143) # it will set a fixed number for all the random-based permutational tests, so it can be reproducible
  #you can use: R16S_Genus_clean, RK0_Genus_clean, RK1_Genus_clean
```
Get the ordination using NMDS and bray curtis distance matrix. FOR PLOT
```{r}
  # R16S_Genus_clean <- prune_samples(!(sample_names(R16S_Genus) %in% c("R0015", "R0127")), R16S_Genus)
  ord <- ordinate(ps, "NMDS", "bray")
```
Recover object -- Give a proper name to the ordination
```{r}
# change the left part of the equivalence 
  RK0_ordination = ord
```

Define your CSS-normalized genus-agglomerate phyloseq object from which we will get the distance matrix (dist)
```{r}
# Change the right part of the equivalence 
  ps = RK1_Genus_clean
  #you can use: R16S_Genus_clean, RK0_Genus_clean, RK1_Genus_clean
```
Get the distance matrix. FOR STATISTICS (permanova)
```{r}
  dist = vegdist(t(otu_table(ps)), method ="bray" )
```
Recover object -- Give a proper name to the distance matrix (e.g., Rumen_dist, Rumen_dist1, etc)
```{r}
# change the left part of the equivalence 
  RK1_dist = dist
```

## 7.2. Plotting beta diversity
Define the variables for the plot:
FOR colorby = collection_day 
  color = c("brown", "#D98326", "#E8B600")
FOR colorby = weaning_group or weaning
  color = c("#967bba", "#505359", "#c9c8c8")
FOR colorby = castration_group / castration
  color = c("#f54545", "#85a9ca", "#e0b54e", "#73a054", "#c9c8c8")
FOR shapeby = collection_day
  shape = c(15,16,17)
```{r}
# Change the right part of the formula
  ps = RK0_Genus_clean
  ord = RK0_ordination
  colorby = "weaning"
  color = c("#967bba", "#505359", "#c9c8c8")
  shapeby = "collection_day" # OPTIONAL
  shape = c(15,16,17) # OPTIONAL
```
Exploratory figure. Note that you will need to change the dataframe "samples_2" according to your needs
```{r}
p<-plot_ordination(ps, ord, type="sample", color=colorby,
                   shape="collection_day")+
  geom_point(size=4.5, alpha=0.8)+
  scale_fill_manual(values = color)+
  scale_shape_manual(values = shape)+
  scale_color_manual(values = color)+
  stat_ellipse()+
  theme_classic() 

p$layers
p$layers <- p$layers[-1]
p
```
## FIGURE 6B
```{r}
levels(sample_data(ps)$weaning) <- c("Fence_line", "Truck","Before_weaning") # Re-labeling data

F6B3 <-plot_ordination(ps, ord, type="sample", color=colorby,
                   shape="collection_day")+
  geom_point(size=4.5, alpha=0.8)+
  scale_fill_manual(values = color)+
  scale_shape_manual(values = shape)+
  scale_color_manual(values = color)+
  stat_ellipse()+
  theme_classic() 

F6B3$layers
F6B3$layers <- F6B3$layers[-1]

#F6B1
#F6B2
#F6B3
```

## 7.3. Statistics
We used a use PERMANOVA (adonis2) to determine the variation explained by our predictors (castration / castration_group, weaning / weaning_group, collection_day)
* Our outcome (response variable) will be distance matrix 
* Our predictor (explanatory variable) will be the interventions and time: castration, weaning and collection_days

Define data frames (metadata)
```{r}
md_16S = data.frame(sample_data(R16S_Genus_clean))
md_K0 = data.frame(sample_data(RK0_Genus_clean))
md_K1 = data.frame(sample_data(RK1_Genus_clean))

set.seed(143) # it will set a fixed number for all the random-based permutational tests, so it can be reproducible
```
Run model for WEANING
```{r}
weaning_16S= adonis2(R16S_dist ~ weaning_group  + collection_day + castration, data=md_16S , by="margin")
weaning_K0= adonis2(RK0_dist ~ weaning_group  + collection_day + castration, data=md_K0 , by="margin")
weaning_K1= adonis2(RK1_dist ~ weaning_group  + collection_day + castration, data=md_K1 , by="margin")
# Rumen_dist_clean = removed R0048 because did not allow to get ordination for collection_day 1 and improved the figure but DID NOT CHANGE THE DIST MATRIX NOR STATISTICS 

```

## 7.4. Procrustes anlysis
```{r}
# more details in: https://john-quensen.com/tutorials/procrustes-analysis/
# and: https://stackoverflow.com/questions/30325739/ggplot2-for-procrustes-rotation-in-vegan
# decent interpretation for protest: https://github.com/vegandevs/vegan/issues/335

# transform bray curtis to euclidean distance required for analysis. We can use pcoa instead of NMDS
add_16S <-  !(is.euclid(R16S_dist))
add_K0 <-  !(is.euclid(RK0_dist))
add_K1 <-  !(is.euclid(RK1_dist))
pcoa.16S <- cmdscale(R16S_dist, k = nrow(t(otu_table(R16S_Genus_clean)))-1, eig = TRUE, add = add_16S)
pcoa.K0 <- cmdscale(RK0_dist, k = nrow(t(otu_table(RK0_Genus_clean)))-1, eig = TRUE, add = add_K0)
pcoa.K1 <- cmdscale(RK1_dist, k = nrow(t(otu_table(RK1_Genus_clean)))-1, eig = TRUE, add = add_K1)
# run procrustes analysis in pairs
pro1 <- procrustes(pcoa.16S,pcoa.K0,  symmetric = TRUE, scores = "sites",choices=c(1,2))
pro2 <- procrustes(pcoa.16S,pcoa.K1,  symmetric = TRUE, scores = "sites",choices=c(1,2))
pro3 <- procrustes(pcoa.K1,pcoa.K0,  symmetric = TRUE, scores = "sites",choices=c(1,2))
# run statistical test
protest1<-protest(X = pcoa.16S, Y = pcoa.K0, scores = "sites", permutations = 999)
protest2<-protest(X = pcoa.16S, Y = pcoa.K1, scores = "sites", permutations = 999)
protest3<-protest(X = pcoa.K1, Y = pcoa.K0, scores = "sites", permutations = 999)
# classic plot
plot(pro3, kind=1, type = "points")
# personalized plot
# create a new data frame for ggplot
ctest1<- data.frame(rda1=pro1$Yrot[,1],rda2=pro1$Yrot[,2],xrda1=pro1$X[,1],xrda2=pro1$X[,2])
ctest1<-merge(ctest1, md_16S, by ='row.names', all = TRUE)
ctest2<- data.frame(rda1=pro2$Yrot[,1],rda2=pro2$Yrot[,2],xrda1=pro2$X[,1],xrda2=pro2$X[,2])
ctest2<-merge(ctest2, md_16S, by ='row.names', all = TRUE)
ctest3<- data.frame(rda1=pro3$Yrot[,1],rda2=pro3$Yrot[,2],xrda1=pro3$X[,1],xrda2=pro3$X[,2])
ctest3<-merge(ctest3, md_K0, by ='row.names', all = TRUE)
```
## FIGURE 7
```{r}
# ggplot it
F7A <- ggplot(ctest1) +
geom_point(aes(x=rda1, y=rda2, colour=weaning)) +
geom_point(aes(x=xrda1, y=xrda2, colour=weaning)) +
geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2,colour=weaning),arrow=arrow(length=unit(0.2,"cm")))+
scale_color_manual(values=c("#967bba", "#505359", "#c9c8c8"))+
  theme_bw()

F7B <- ggplot(ctest2) +
geom_point(aes(x=rda1, y=rda2, colour=weaning)) +
geom_point(aes(x=xrda1, y=xrda2, colour=weaning)) +
geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2,colour=weaning),arrow=arrow(length=unit(0.2,"cm")))+
scale_color_manual(values=c("#967bba", "#505359", "#c9c8c8"))+
  theme_bw()

F7C <- ggplot(ctest3) +
geom_point(aes(x=rda1, y=rda2, colour=weaning)) +
geom_point(aes(x=xrda1, y=xrda2, colour=weaning)) +
geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2,colour=weaning),arrow=arrow(length=unit(0.2,"cm")))+
scale_color_manual(values=c("#967bba", "#505359", "#c9c8c8"))+
  theme_bw()
```

