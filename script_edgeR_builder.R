library(tidyverse)
library(RColorBrewer)
library(edgeR)
library(pheatmap)
library(DEGreport)
library(ashr)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(dplyr)


create_result_file <-function(resTable, resultFile, lfcCutoff=0,outputname)
{
  #takes the table produced after DGE, the true result file to which it will be compared, and the lfc cutoff to be used
  p_cutoff <- 0.05
  lfc_cutoff <- lfcCutoff
  trueResult <-read.table(paste("Meta/",resultFile,sep = ""), header=T, row.names=1)
  #filter by p value and lfc
  resTable <- as.data.frame(resTable)
  res_table_filtered <- resTable %>%
    filter(PValue < p_cutoff & abs(logFC) > lfc_cutoff)
  #create final output file
  result_log_filtered <- data.frame(
    gene = rownames(trueResult),
    upregulation = rep(0, nrow(trueResult)),
    downregulation = rep(0, nrow(trueResult)),
    differential.expression = rep(0, nrow(trueResult))
  )
  #check results in res_tableCE_TB and update result_meta_df
  for (i in 1:nrow(res_table_filtered)) 
  {
    if (!is.na(res_table_filtered$PValue[i])) 
    {
      #The following line are to make sure the loop works even if the genes are not ordered
      numeric_value <- as.numeric(sub("^g([0-9]+).*", "\\1",rownames(res_table_filtered[i,])))
      result_log_filtered$differential.expression[numeric_value] <- 1
      if(res_table_filtered$logFC[i]>0)
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
  output_file <-  paste0("Results/edger/edgeR_", string_output,'_',outputname, "_Results.txt")
  write.table(  result_log_filtered, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}


EdgeREvaluator <-function(dataFile, metaFile, resultFile)
{
  data <-read.table(paste("Data/",dataFile,sep = ""), header=T, row.names=1)
  meta <- read.table(paste("MetaFull/",metaFile,sep = ""), header=T, row.names=1)
  
  #create the DGEList object, takes at leat a count matrix
  samples <- substr(dataFile, 1, nchar(dataFile) - 4)
  split_values <- strsplit(dataFile, "_")[[1]]
  #extract the forst number
  num_samples <- as.integer(split_values[1])
  num_conditions <- 2
  groups <-c(rep(1, num_samples), rep(2, num_samples))
  y<- DGEList(counts = data, group=groups)
  #filter genes with low counts, remove rows a small count(the smallest group size)
  keep <- filterByExpr(y)
  y_filterd_low_count <- y 
  y_filtered_low_count <- y[keep, , keep.lib.sizes=FALSE]
  #normalization by the size of the library, TMM
  y <- normLibSizes(y)
  y_filtered_low_count<-normLibSizes(y_filtered_low_count)
  #estimate dispersion
  y <- estimateDisp(y)
  y_filtered_low_count<-estimateDisp(y_filtered_low_count)
  #Testing for DE genes, if a single factor exact test is possible
  et <- exactTest(y)
  et_filtered_low_count<-exactTest(y_filtered_low_count)

  create_result_file(et$table,resultFile, 0.58,'log_filtered')
  create_result_file(et_filtered_low_count$table,resultFile, 0.58, 'low_log_filtered')
  create_result_file(et$table,resultFile, 0,'no_filtered')
  create_result_file(et_filtered_low_count$table,resultFile, 0,'low_filtered')
}
##############################################################################
EdgeREvaluator("3_500_500.tsv","3_500_500_full.txt","3_500_500_meta.tsv")
EdgeREvaluator("3_750_250.tsv","3_750_250_full.txt","3_750_250_meta.tsv")
EdgeREvaluator("3_1000_0.tsv","3_1000_0_full.txt", "3_1000_0_meta.tsv")
EdgeREvaluator("6_500_500.tsv","6_500_500_full.txt", "6_500_500_meta.tsv")
EdgeREvaluator("6_750_250.tsv","6_750_250_full.txt", "6_750_250_meta.tsv")
EdgeREvaluator("6_1000_0.tsv","6_1000_0_full.txt", "6_1000_0_meta.tsv")
EdgeREvaluator("9_500_500.tsv","9_500_500_full.txt", "9_500_500_meta.tsv")
EdgeREvaluator("9_750_250.tsv","9_750_250_full.txt","9_750_250_meta.tsv")
EdgeREvaluator("9_1000_0.tsv","9_1000_0_full.txt", "9_1000_0_meta.tsv")
##############################################################################

evaluate_DE <- function(calculated_data, results, method)
{
  ##takes the file with the result, the file with the calculated results, the method used 
  result_file <- file.path("Results", method, calculated_data)
  answer_file <- file.path("Meta", results)
  results <- read.table(result_file, header=T, row.names=1)
  answers <- read.table(answer_file, header=T, row.names=1)
 
  truePositive <- 0
  trueNegative <- 0
  falsePositive <-0
  falseNegative <- 0
  
  for(i in 1:nrow(results))
  {
    row1 <- results[i,]
    row2 <- answers[i,]
    {
      if (row1$differential.expression == 1 && row2$differential.expression == 1) {
        truePositive <- truePositive + 1
      } else if (row1$differential.expression == 0 && row2$differential.expression == 0) {
        trueNegative <- trueNegative + 1
      } else if (row1$differential.expression == 1 && row2$differential.expression == 0) {
        falsePositive <- falsePositive + 1
      } else if (row1$differential.expression == 0 && row2$differential.expression == 1) {
        falseNegative <- falseNegative + 1
      }
    }
  }
  recall <- truePositive/(truePositive+falseNegative)
  accuracy <- (truePositive+trueNegative)/(trueNegative+truePositive+falseNegative+falsePositive)
  precision <- truePositive/(truePositive+falsePositive)
  falsePositiveRate <- falsePositive/(falsePositive+trueNegative)
  
  
  result_vector <- c("True Positives:" = truePositive,
                     "True Negatives:"= trueNegative,
                     "False Positives:" = falsePositive,
                     "False Negatives:" = falseNegative,
                     "Recall:" = recall,
                     "Accuracy:" =accuracy,
                     "False Positive Rate:" = falsePositiveRate)
  return(result_vector)
  
}
results_list <- list(evaluate_DE("3_500_500_Results.txt", "3_500_500_meta.tsv", "edgeR"),
                     evaluate_DE("3_750_250_Results.txt", "3_750_250_meta.tsv", "edgeR"),
                     evaluate_DE("3_1000_0_Results.txt", "3_1000_0_meta.tsv", "edgeR"),
                     evaluate_DE("6_500_500_Results.txt", "6_500_500_meta.tsv", "edgeR"),
                     evaluate_DE("6_750_250_Results.txt", "6_750_250_meta.tsv", "edgeR"),
                     evaluate_DE("6_1000_0_Results.txt", "6_1000_0_meta.tsv", "edgeR"),
                     evaluate_DE("9_500_500_Results.txt", "9_500_500_meta.tsv", "edgeR"),
                     evaluate_DE("9_750_250_Results.txt", "9_750_250_meta.tsv", "edgeR"),
                     evaluate_DE("9_1000_0_Results.txt", "9_1000_0_meta.tsv", "edgeR")
)

accuracy_values <- c()
recall_values <- c()

colors <- c("black", "orange", "red")
#extract the appropiate values
for (result in results_list) {
  accuracy_values <- c(accuracy_values, result[["Accuracy:"]])
  recall_values <- c(recall_values, result[["Recall:"]])
}
# Add vertical dashed lines at specific x-values (n=3 and n=6)
x_values <- c(0.915, 0.95)  # Example x-values where you want to add vertical lines

plot(accuracy_values, recall_values, type = "p", pch = 19, col = colors,
     xlab = "Accuracy", ylab = "Recall", main = "Accuracy vs. Recall")
#add legend
legend("topleft", legend = c("500_500", "750_250", "1000_0"), fill = colors)
#make it pretty
annotations <- c("n=3 ", "n=6", "n=9")  # Annotations corresponding to the vertical lines
text(x = c(0.910, 0.925, 0.934), y = c(0.55, 0.55, 0.55), labels = annotations, pos = 1, col = "gray")
abline(v = x_values, col = "gray", lty = 2)

EdgeREvaluator("3_500_500.tsv","3_500_500_full.txt")
EdgeREvaluator("3_750_250.tsv","3_750_250_full.txt")
EdgeREvaluator("3_1000_0.tsv","3_1000_0_full.txt")
EdgeREvaluator("6_500_500.tsv","6_500_500_full.txt")
EdgeREvaluator("6_750_250.tsv","6_750_250_full.txt")
EdgeREvaluator("6_1000_0.tsv","6_1000_0_full.txt")
EdgeREvaluator("9_500_500.tsv","9_500_500_full.txt")
EdgeREvaluator("9_750_250.tsv","9_750_250_full.txt")
EdgeREvaluator("9_1000_0.tsv","9_1000_0_full.txt")
sessionInfo()