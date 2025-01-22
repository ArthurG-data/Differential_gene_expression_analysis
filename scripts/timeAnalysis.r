##The goal is to compare runtime for different tools

library(DESeq2)
library(dplyr)
library(tibble)
library(edgeR)
library(limma)
library(RColorBrewer)
library(caret)

library(patchwork)
library(tidyverse)
library(NOISeq)
library(ggplot2)
library(readr)
library(tidyr)
library(data.table)
library(Matrix)
library(ggrepel)
library(gridExtra)
############################################################################################################################
#time efficiency of deseq2
time_evaluation_desq2<-function(matrixcount,metaFile)
{
  #take one file and return the time to do the DEF
  start <- Sys.time()
  ##start with DESeq2
  data <- read.table(paste("Data/",matrixcount,sep = ""), header=T, row.names=1)
  meta <- read.table(paste("MetaFull/",metaFile,sep = ""), header=T, row.names=1)
  dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~Condition)
  # when have to perform the median of ratio of normalization
  dds <- estimateSizeFactors(dds)
  #We then retrieve the normalize count matrix
  normalized_counts <- counts(dds, normalized=TRUE)
  dds <- DESeq(dds)
  contrast_ce <- c("Condition", "cond2", "cond1")
  res_tableCE_unshrunken <- results(dds, contrast=contrast_ce, alpha = 0.05)
  res_tableCE <- lfcShrink(dds, contrast=contrast_ce, res=res_tableCE_unshrunken,type="normal")
  #we use filter to get only the significant
  res_tableCE_tb <- res_tableCE %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()

  padj.cutoff <- 0.05
  lfc.cutoff <-0.58
  res_tableCE_tb <- res_tableCE_tb %>%
    mutate(threshold_CE = padj <  padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
  return  (Sys.time() - start)
}
################################################################################################################################
#time efficiency of edgeR
time_evaluation_edgeR<-function(matrixcount,metaFile)
{
  start <- Sys.time()
  data <-read.table(paste("Data/",matrixcount,sep = ""), header=T, row.names=1)
  meta <- read.table(paste("MetaFull/",metaFile,sep = ""), header=T, row.names=1)
  #create the DGEList object, takes at leat a count matrix
  samples <- substr(matrixcount, 1, nchar(matrixcount) - 4)
  split_values <- strsplit(matrixcount, "_")[[1]]
  #extract the forst number
  num_samples <- as.integer(split_values[1])
  num_conditions <- 2
  groups <-c(rep(1, num_samples), rep(2, num_samples))
  y<- DGEList(counts = data, group=groups)
  #filter genes with low counts, remove rows a small count(the smallest group size)
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  #normalization by the size of the library
  y <- normLibSizes((y))
  #estimate dispersion
  y <- estimateDisp(y)
  et <- exactTest(y)
  padj.cutoff <- 0.05
  lfc.cutoff <-0.58
  et_table = as.data.frame(et$table)
  filtered_genes <- et_table[et_table$padj < padj.cutoff & abs(et_table$logFC) > lfc.cutoff ]
  return (Sys.time() - start)
}
###########################################################################################################
#time efficiency limma
time_evaluation_limma<-function(matrixcount,metaFile)
{
  start <- Sys.time()
  data <-read.table(paste("Data/",matrixcount,sep = ""), header=T, row.names=1)
  meta <- read.table(paste("MetaFull/",metaFile,sep = ""), header=T, row.names=1)
  
  dg <- DGEList(data)

  #Creating groups
  if (substring(metaFile, 1, 1) == "3"){
    group <- as.factor(c("Condition1", "Condition1", "Condition1", 
                         "Condition2", "Condition2", "Condition2"))  
  }
  else if (substring(metaFile, 1, 1) == "6"){
    group <- as.factor(c("Condition1", "Condition1", "Condition1", 
                         "Condition1", "Condition1", "Condition1",
                         "Condition2", "Condition2", "Condition2",
                         "Condition2", "Condition2", "Condition2"))
  }
  else if (substring(metaFile, 1, 1) == "9"){
    group <- as.factor(c("Condition1", "Condition1", "Condition1", 
                         "Condition1", "Condition1", "Condition1",
                         "Condition1", "Condition1", "Condition1",
                         "Condition2", "Condition2", "Condition2",
                         "Condition2", "Condition2", "Condition2",
                         "Condition2", "Condition2", "Condition2"))
  }
  
  dg$samples$group <- group
  ## Data preprocessing ####
  
  # Removing genes that are lowly expressed
  keep.exprs <- filterByExpr(dg, group=group)
  dg <- dg[keep.exprs,, keep.lib.sizes=FALSE]
  dim(dg)
  
  #normalisation
  dg <- calcNormFactors(dg, method = "TMM")
  dg$samples$norm.factors
  
  #Design matrix
  design <- model.matrix(~0+group)
  colnames(design) <- gsub("group", "", colnames(design))
  design
  
  
  contr.matrix <- makeContrasts(
    Cnd2vsCnd1 = Condition2-Condition1, 
    levels = colnames(design))
  contr.matrix
  
  vdg <- voom(dg, design, plot=FALSE)
  vdg
  
  #Fitting linear model
  vfit <- lmFit(vdg, design)
  vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit)
  
  #Table with pvalue and mapping differentially expressed genes
  top.table <- topTable(efit, sort.by = "P", n = Inf) %>%
    mutate(upregulation=ifelse(adj.P.Val <0.05  & logFC > 0.58, 1, 0)) %>%
    mutate(downregulation=ifelse(adj.P.Val <0.05  & logFC < 0.58, 1, 0)) %>%
    mutate(differential.expression=ifelse(upregulation ==1|
                                            downregulation ==1, 1, 0)) 
  
  top.table <- top.table[,!names(top.table) %in% 
                           c("logFC","AveExpr","t","P.Value","B")]
  
  
  return (Sys.time() - start)
}
###########################################################################################################
#time study NOISeq
#ER_boot_list_fold[[j]] <- rownames(ER_res[ER_res$log2FC > 2,])
time_evaluation_NOISeq<-function(matrixcount,metaFile)
{
  start <- Sys.time()
  data <-read.table(paste("Data/",matrixcount,sep = ""), header=T, row.names=1)
  meta <- read.table(paste("MetaFull/",metaFile,sep = ""), header=T, row.names=1)
  NoisData<-readData(data=data, factors=meta)

  mynoiseq = noiseqbio(NoisData ,norm = "uqua", factor = "Condition", lc=0)

  res_tb<-mynoiseq@results %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble() 
  #q : threshold for q ,set to 0.95 for Noiseqbio
  q<-0.8
  mynoiseq.deg = degenes(mynoiseq, q = q, M = NULL)


  mynoiseq.deg1 = degenes(mynoiseq, q = q, M = "up")
  mynoiseq.deg2 = degenes(mynoiseq, q = q, M = "down")
  head(mynoiseq@results[[1]])
  mynoiseq.deg<-mynoiseq.deg%>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  print(mynoiseq)
  return (Sys.time() - start)
}

###########################################################################################################
##

time_n3_deseq2<-0
time_n6_deseq2<-0
time_n9_deseq2<-0

time_n6_edgeR<-0
time_n3_edgeR<-0
time_n9_edgeR<-0

time_n3_limma<-0
time_n6_limma<-0
time_n9_limma<-0

time_n3_NOISeq<-0
time_n6_NOISeq<-0
time_n9_NOISeq<-0


for(i in 1:10)
{
  print(i)
  time_n3_deseq2<-time_n3_deseq2+ time_evaluation_desq2("3_500_500.tsv","3_500_500_full.txt")+time_evaluation_desq2("3_750_250.tsv","3_750_250_full.txt")+time_evaluation_desq2("3_1000_0.tsv","3_1000_0_full.txt")
  time_n6_deseq2<- time_n6_deseq2 + time_evaluation_desq2("6_500_500.tsv","6_500_500_full.txt")+time_evaluation_desq2("6_750_250.tsv","6_750_250_full.txt")+time_evaluation_desq2("6_1000_0.tsv","6_1000_0_full.txt")
  time_n9_deseq2<-time_n9_deseq2 + time_evaluation_desq2("9_500_500.tsv","9_500_500_full.txt")+time_evaluation_desq2("9_750_250.tsv","9_750_250_full.txt")+time_evaluation_desq2("9_1000_0.tsv","9_1000_0_full.txt")
  #time_n3_edgeR<-time_n3_edgeR + time_evaluation_edgeR("3_500_500.tsv","3_500_500_full.txt")+time_evaluation_edgeR("3_750_250.tsv","3_750_250_full.txt")+time_evaluation_edgeR("3_1000_0.tsv","3_1000_0_full.txt")
  #time_n6_edgeR<-time_n6_edgeR + time_evaluation_edgeR("6_500_500.tsv","6_500_500_full.txt")+time_evaluation_edgeR("6_750_250.tsv","6_750_250_full.txt")+time_evaluation_edgeR("6_1000_0.tsv","6_1000_0_full.txt")
  #time_n9_edgeR<-time_n9_edgeR + time_evaluation_edgeR("9_500_500.tsv","9_500_500_full.txt")+time_evaluation_edgeR("9_750_250.tsv","9_750_250_full.txt")+time_evaluation_edgeR("9_1000_0.tsv","9_1000_0_full.txt")
  #time_n3_limma<-time_n3_limma + time_evaluation_limma("3_500_500.tsv","3_500_500_full.txt")+time_evaluation_limma("3_750_250.tsv","3_750_250_full.txt")+time_evaluation_limma("3_1000_0.tsv","3_1000_0_full.txt")
  #time_n6_limma<-time_n6_limma + time_evaluation_limma("6_500_500.tsv","6_500_500_full.txt")+time_evaluation_limma("6_750_250.tsv","6_750_250_full.txt")+time_evaluation_limma("6_1000_0.tsv","6_1000_0_full.txt")
  #time_n9_limma<-time_n9_limma + time_evaluation_limma("9_500_500.tsv","9_500_500_full.txt")+time_evaluation_limma("9_750_250.tsv","9_750_250_full.txt")+time_evaluation_limma("9_1000_0.tsv","9_1000_0_full.txt")
  #time_n3_NOISeq<-time_n3_NOISeq + time_evaluation_NOISeq("3_500_500.tsv","3_500_500_full.txt")+time_evaluation_NOISeq("3_750_250.tsv","3_750_250_full.txt")+time_evaluation_NOISeq("3_1000_0.tsv","3_1000_0_full.txt")
  #time_n6_NOISeq<-time_n6_NOISeq + time_evaluation_NOISeq("6_500_500.tsv","6_500_500_full.txt")+time_evaluation_NOISeq("6_750_250.tsv","6_750_250_full.txt")+time_evaluation_NOISeq("6_1000_0.tsv","6_1000_0_full.txt")
  #time_n9_NOISeq<-time_n9_NOISeq + time_evaluation_NOISeq("9_500_500.tsv","9_500_500_full.txt")+time_evaluation_NOISeq("9_750_250.tsv","9_750_250_full.txt")+time_evaluation_NOISeq("9_1000_0.tsv","9_1000_0_full.txt")
}

your_data_frame <- data.frame(
  Tool = rep(c("DESeq2", "edgeR", "limma", "NOISeq"), each = 3),
  NumberOfSamples = rep(c(3, 6, 9), times = 4),
  Time = c(
    time_n3_deseq2, time_n6_deseq2, time_n9_deseq2,
    time_n3_edgeR, time_n6_edgeR, time_n9_edgeR,
    time_n3_limma, time_n6_limma, time_n9_limma,
    time_n3_NOISeq, time_n6_NOISeq, time_n9_NOISeq
  )/30
)
custom_colors <- c("edgeR" = "#7CAE00", "DESeq2" = "#F8766D", "NOISeq" = "#00BFC4","limma" = "#C77CFF" )
print(your_data_frame)
ggplot(your_data_frame, aes(x = NumberOfSamples, y = Time, color = Tool)) +
  geom_point(size = 3) +
  geom_line() +
  labs(
    x = "No of Samples",
    y = "Run time(sec)",
    color = NULL
  ) +
  scale_color_manual(values = custom_colors)+
  ggtitle("Runtime Analysis") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))
