
library(NOISeq)
library(dplyr)
library(tibble)
create_result_file <-function(resTable, resultFile, lfcCutoff=0,outputname)
{

  lfc_cutoff <- lfcCutoff
  trueResult <-read.table(paste("Meta/",resultFile,sep = ""), header=T, row.names=1)
  
  res_table_filtered <- resTable %>%
    filter(prob >= 0.95 & abs(log2FC) > lfc_cutoff)
  
  result_log_filtered <- data.frame(
    gene = rownames(trueResult),
    upregulation = rep(0, nrow(trueResult)),
    downregulation = rep(0, nrow(trueResult)),
    differential.expression = rep(0, nrow(trueResult))
  )
  #check results in res_tableCE_TB and update result_meta_df
  for (i in 1:nrow(res_table_filtered)) 
  {
    if (!is.na(res_table_filtered$prob[i])) 
    {
      numeric_value <- as.numeric(sub("^g([0-9]+).*", "\\1",res_table_filtered$gene[i]))
      result_log_filtered$differential.expression[numeric_value] <- 1
      if(res_table_filtered$log2FC[i]<0)
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
  output_file <-  paste0("Results/NOISeq/NOISeq_", string_output,'_',outputname, "_Results.txt")
  write.table(  result_log_filtered, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

NOISeqEvaluator<-function(dataFile, metaFile, resultFile)
 {
  data <-read.table(paste("Data/",dataFile,sep = ""), header=T, row.names=1)
  meta <- read.table(paste("MetaFull/",metaFile,sep = ""), header=T, row.names=1)

  NoisData<-readData(data=data, factors=meta)
  
  mynoiseq_low_filtered = noiseqbio(NoisData ,k = 0.5,norm = "rpkm", factor = "Condition", lc=1, filter=1, cpm=1, r = 20, adj = 1.5, nclust = 15,a0per = 0.9, random.seed = 12345)
  mynoiseq = noiseqbio(NoisData ,norm = "rpkm",k = 0.5, factor = "Condition", lc = 1, r = 20, adj = 1.5, filter=0,nclust = 15,plot = FALSE, a0per = 0.9, random.seed = 12345)
  
  #q : threshold for q ,set to 0.95 for Noiseqbio, the treshold will be applied in the create rseulf file function
  q<-0

  mynoiseq.deg = degenes(mynoiseq, q = q, M = NULL)
  mynoiseq_low_filtered.deg= degenes(mynoiseq_low_filtered, q = q, M = NULL)
  
  mynoiseq_low_filtered<-mynoiseq_low_filtered.deg%>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  mynoiseq_no<-mynoiseq.deg%>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  create_result_file(mynoiseq_no,resultFile,lfcCutoff=0,'no_filtered')
  create_result_file(mynoiseq_no,resultFile,lfcCutoff=0.58,'log_filtered')
  create_result_file(mynoiseq_low_filtered,resultFile,lfcCutoff=0,'low_filtered')
  create_result_file(mynoiseq_low_filtered,resultFile,lfcCutoff=0.58,'low_log_filtered')
}
NOISeqEvaluator("3_500_500.tsv","3_500_500_full.txt", "3_500_500_meta.tsv")
NOISeqEvaluator("3_750_250.tsv","3_750_250_full.txt", "3_750_250_meta.tsv")
NOISeqEvaluator("3_1000_0.tsv","3_1000_0_full.txt", "3_1000_0_meta.tsv")
NOISeqEvaluator("6_500_500.tsv","6_500_500_full.txt", "6_500_500_meta.tsv")
NOISeqEvaluator("6_750_250.tsv","6_750_250_full.txt", "6_750_250_meta.tsv")
NOISeqEvaluator("6_1000_0.tsv","6_1000_0_full.txt", "6_1000_0_meta.tsv")
NOISeqEvaluator("9_500_500.tsv","9_500_500_full.txt", "9_500_500_meta.tsv")
NOISeqEvaluator("9_750_250.tsv","9_750_250_full.txt", "9_750_250_meta.tsv")
NOISeqEvaluator("9_1000_0.tsv","9_1000_0_full.txt", "9_1000_0_meta.tsv")

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
  print(result_vector)
  return(result_vector)
}
results_list <- list(evaluate_DE("3_500_500_Results.txt", "3_500_500_meta.tsv", "NOISeq"),
                     evaluate_DE("3_750_250_Results.txt", "3_750_250_meta.tsv", "NOISeq"),
                     evaluate_DE("3_1000_0_Results.txt", "3_1000_0_meta.tsv", "NOISeq"),
                     evaluate_DE("6_500_500_Results.txt", "6_500_500_meta.tsv", "NOISeq"),
                     evaluate_DE("6_750_250_Results.txt", "6_750_250_meta.tsv", "NOISeq"),
                     evaluate_DE("6_1000_0_Results.txt", "6_1000_0_meta.tsv", "NOISeq"),
                     evaluate_DE("9_500_500_Results.txt", "9_500_500_meta.tsv", "NOISeq"),
                     evaluate_DE("9_750_250_Results.txt", "9_750_250_meta.tsv", "NOISeq"),
                     evaluate_DE("9_1000_0_Results.txt", "9_1000_0_meta.tsv", "NOISeq")
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
x_values <- c(0.915, 0.932)  # Example x-values where you want to add vertical lines

plot(accuracy_values, recall_values, type = "p", pch = 19, col = colors,
     xlab = "Accuracy", ylab = "Recall", main = "Accuracy vs. Recall")
#add legend
legend("topleft", legend = c("500_500", "750_250", "1000_0"), fill = colors)
#make it pretty
annotations <- c("n=3 ", "n=6", "n=9")  # Annotations corresponding to the vertical lines
text(x = c(0.910, 0.925, 0.934), y = c(0.55, 0.55, 0.55), labels = annotations, pos = 1, col = "gray")
abline(v = x_values, col = "gray", lty = 2)

v3_500_500EdgeR =NOISeqEvaluator("3_500_500.tsv","3_500_500_full.txt")
v3_750_250EdgeR=NOISeqEvaluator("3_750_250.tsv","3_750_250_full.txt")
v3_1000_0EdgeR=NOISeqEvaluator("3_1000_0.tsv","3_1000_0_full.txt")
v6_500_500EdgeR=NOISeqEvaluator("6_500_500.tsv","6_500_500_full.txt")
v6_750_250EdgeR=NOISeqEvaluator("6_750_250.tsv","6_750_250_full.txt")
v6_1000_0EdgeR=NOISeqEvaluator("6_1000_0.tsv","6_1000_0_full.txt")
v9_500_500EdgeR=NOISeqEvaluator("9_500_500.tsv","9_500_500_full.txt")
v9_750_250EdgeR=NOISeqEvaluator("9_750_250.tsv","9_750_250_full.txt")
v9_1000_0EdgeR=NOISeqEvaluator("9_1000_0.tsv","9_1000_0_full.txt")

