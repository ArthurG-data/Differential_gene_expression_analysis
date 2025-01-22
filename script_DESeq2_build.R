## Gene-level differential expression analysis using DESeq
## Setup
### Bioconductor and CRAN libraries used
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(ashr)
library(ggplot2)
library(ggrepel)
library(DEGreport)
##directory toward the folder
setwd("C:/Users/lorie/OneDrive/Bureau/Semester_2/IFN646_BiomedicalDataScience/assessments/project/Project")
## Load in data
MetaCreator <- function(input_string)
{
  #split string on underscore
  string_output <- substr(input_string, 1, nchar(input_string) - 4)
  split_values <- strsplit(input_string, "_")[[1]]
  #extract the forst number
  num_samples <- as.integer(split_values[1])
  num_conditions <- 2
  total_row <- num_conditions*num_samples
  #generate content
  sample_info <- data.frame(
    sampletype = paste0("sample", 1:total_row),
    Condition = rep(paste0("cond", ifelse(1:total_row <= num_samples,1,2)))
  )
  #Specify the output file name
  output_file <-  paste0("MetaFull/", string_output, "_full.txt")
  
  # Write the sample information to a tab-separated text file
  write.table(sample_info, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Print a message indicating that the file has been created
  cat("Sample information has been written to", output_file, "\n")
}

create_result_file <-function(resTable, resultFile, lfcCutoff=0,outputname)
{
  
  padj.cutoff <- 0.05
  lfc.cutoff <- lfcCutoff
  trueResult <-read.table(paste("Meta/",resultFile,sep = ""), header=T, row.names=1)
  
  resTable <- as.data.frame(resTable)
  res_table_filtered <- resTable %>%
    filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
  
  result_log_filtered <- data.frame(
    gene = rownames(trueResult),
    upregulation = rep(0, nrow(trueResult)),
    downregulation = rep(0, nrow(trueResult)),
    differential.expression = rep(0, nrow(trueResult))
  )
  #check results in res_tableCE_TB and update result_meta_df
  for (i in 1:nrow(res_table_filtered)) 
  {
    
    if (!is.na(res_table_filtered$padj[i])) 
    {
      numeric_value <- as.numeric(sub("^g([0-9]+).*", "\\1",rownames(res_table_filtered[i,])))
      result_log_filtered$differential.expression[numeric_value] <- 1
      if(res_table_filtered$log2FoldChange[i]>0)
      {
        result_log_filtered$upregulation[numeric_value] <- 1
      }
      else 
      {
        result_log_filtered$downregulation[numeric_value] <- 1
      }
    }

  }
  string_output <- substr(resultFile, 1, nchar(resultFile) - 4)
  #save the result_file
  output_file <-  paste0("Results/DESeq2/DESeq2_", string_output,'_',outputname, "_Results.txt")
  write.table(  result_log_filtered, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

DESeqEvaluator <-function(dataFile, metaFile,resultFile)
  {
  min_count <- 10

  data <- read.table(paste("Data/",dataFile,sep = ""), header=T, row.names=1)
  meta <- read.table(paste("MetaFull/",metaFile,sep = ""), header=T, row.names=1)
  max_group <- ncol(data)/2

  ## Check that sample names match in both files
  all(colnames(data) %in% rownames(meta))
  all(colnames(data) == rownames(meta))
  ## Create DESeq2Dataset object
  dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~Condition)
  # View the original counts matrix
  
  # when have to perform the median of ratio of normalization
  dds <- estimateSizeFactors(dds)
  #filter low gene with total count less than specified accross all genes
  keep= rowSums( counts(dds) >= min_count ) >= max_group
  filtered_dds <-dds[keep,]
  dds <- DESeq(dds)
 
  dds_filtered <-DESeq(filtered_dds)
  #### Define contrasts, extract results table, and shrink the log2 fold changes
  contrast_ce <- c("Condition", "cond2", "cond1")
  res_tableCE_unshrunken <- results(dds, contrast=contrast_ce, alpha = 0.05)
  res_tableCE_unshrunken_filtered <-results(dds_filtered, contrast=contrast_ce, alpha = 0.05)
  
  res_tableCE <- lfcShrink(dds, contrast=contrast_ce, res=res_tableCE_unshrunken,type="normal")
  res_tableCE_filtered <- lfcShrink(dds_filtered, contrast=contrast_ce, res=res_tableCE_unshrunken_filtered,type="normal")
  
  create_result_file(res_tableCE_filtered,resultFile, 0, 'low_filtered')
  create_result_file(res_tableCE_filtered,resultFile, 0.58, 'low_log_filtered')
  create_result_file(res_tableCE,resultFile, 0.58, 'log_filtered')
  create_result_file(res_tableCE,resultFile, 0, 'no_filtered')
  print("Over")

}
##here wa create the necessary files for the meta
MetaCreator("3_500_500.tsv")
MetaCreator("3_750_250.tsv")
MetaCreator("3_1000_0.tsv")
MetaCreator("6_500_500.tsv")
MetaCreator("6_750_250.tsv")
MetaCreator("6_1000_0.tsv")
MetaCreator("9_500_500.tsv")
MetaCreator("9_750_250.tsv")
MetaCreator("9_1000_0.tsv")
#here we create vectors with the necessary statistical values
DESeqEvaluator("3_500_500.tsv","3_500_500_full.txt","3_500_500_meta.tsv")
DESeqEvaluator("3_750_250.tsv","3_750_250_full.txt","3_750_250_meta.tsv")
DESeqEvaluator("3_1000_0.tsv","3_1000_0_full.txt","3_1000_0_meta.tsv")
DESeqEvaluator("6_500_500.tsv","6_500_500_full.txt","6_500_500_meta.tsv")
DESeqEvaluator("6_750_250.tsv","6_750_250_full.txt","6_750_250_meta.tsv")
DESeqEvaluator("6_1000_0.tsv","6_1000_0_full.txt","6_1000_0_meta.tsv")
DESeqEvaluator("9_500_500.tsv","9_500_500_full.txt","9_500_500_meta.tsv")
DESeqEvaluator("9_750_250.tsv","9_750_250_full.txt","9_750_250_meta.tsv")
DESeqEvaluator("9_1000_0.tsv","9_1000_0_full.txt","9_1000_0_meta.tsv")


sessionInfo()