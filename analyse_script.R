#creates the necessary graph from the files

setwd("C:/Users/lorie/OneDrive - Queensland University of Technology/Semester_2/IFN646_BiomedicalDataScience/assessments/project/Project")
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(dplyr)
library(ggVennDiagram)
library(dplyr)
library(ggpubr)
#work with the highest number samples, with 500 500

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
  
  for(i in 1:nrow(answers))
  {
    row1 <- results[i,]
    row2 <- answers[i,]
    {
      if(i< nrow(results)){
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
  }
  recall <- truePositive/(truePositive+falseNegative)
  accuracy <- (truePositive+trueNegative)/(trueNegative+truePositive+falseNegative+falsePositive)
  precision <- truePositive/(truePositive+falsePositive)
  falsePositiveRate <- falsePositive/(falsePositive+trueNegative)
  trueNegativeRate <- trueNegative/(trueNegative+falsePositive)
  sensitivity <- truePositive/(truePositive+falseNegative)
  specificity<-trueNegative/(trueNegative+falsePositive)

  
  
  result_vector <- c("True Positives:" = truePositive,
                     "True Negatives:"= trueNegative,
                     "False Positives:" = falsePositive,
                     "False Negatives:" = falseNegative,
                     "Recall:" = recall,
                     "Accuracy:" =accuracy,
                     "False Positive Rate:" = falsePositiveRate,
                     "True Negative Rate:"= trueNegativeRate,
                    "Sensitivity:"=sensitivity,
                    "specificity:"=specificity)
  print(result_vector)
  return(result_vector)
}

transform_rowname <- function(rowname) {
  first_char <- substr(rowname, 1, 1)  # Extract the first character
  rest_chars <- as.integer(substr(rowname, 2, nchar(rowname)))  # Convert the rest to integers
  
  return(rest_chars)
}
####################################################################################
#DESeq_low_log_filtered
v3_500_500DESeq =evaluate_DE("DESeq2_3_500_500_meta_low_log_filtered_Results.txt", "3_500_500_meta.tsv","DESeq2")
v3_750_250DESeq=evaluate_DE("DESeq2_3_750_250_meta_low_log_filtered_Results.txt", "3_750_250_meta.tsv","DESeq2")
v3_1000_0DESeq=evaluate_DE("DESeq2_3_1000_0_meta_low_log_filtered_Results.txt", "3_1000_0_meta.tsv","DESeq2")
v6_500_500DESeq=evaluate_DE("DESeq2_6_500_500_meta_low_log_filtered_Results.txt", "6_500_500_meta.tsv","DESeq2")
v6_750_250DESeq=evaluate_DE("DESeq2_6_750_250_meta_low_log_filtered_Results.txt", "6_750_250_meta.tsv","DESeq2")
v6_1000_0DESeq=evaluate_DE("DESeq2_6_1000_0_meta_low_log_filtered_Results.txt", "6_1000_0_meta.tsv","DESeq2")
v9_500_500DESeq=evaluate_DE("DESeq2_9_500_500_meta_low_log_filtered_Results.txt", "9_500_500_meta.tsv","DESeq2")
v9_750_250DESeq=evaluate_DE("DESeq2_9_750_250_meta_low_log_filtered_Results.txt", "9_750_250_meta.tsv","DESeq2")
v9_1000_0DESeq=evaluate_DE("DESeq2_9_1000_0_meta_low_log_filtered_Results.txt", "9_1000_0_meta.tsv","DESeq2")
#NOISeq_low_filtered
v3_500_500noiseq = evaluate_DE("NOISeq_3_500_500_meta_low_log_filtered_Results.txt", "3_500_500_meta.tsv", "NOISeq")
v3_750_250noiseq = evaluate_DE("NOISeq_3_750_250_meta_low_log_filtered_Results.txt", "3_750_250_meta.tsv", "NOISeq")
v3_1000_0noiseq = evaluate_DE("NOISeq_3_1000_0_meta_low_log_filtered_Results.txt", "3_1000_0_meta.tsv", "NOISeq")
v6_500_500noiseq = evaluate_DE("NOISeq_6_500_500_meta_low_log_filtered_Results.txt", "6_500_500_meta.tsv", "NOISeq")
v6_750_250noiseq = evaluate_DE("NOISeq_6_750_250_meta_low_log_filtered_Results.txt", "6_750_250_meta.tsv", "NOISeq")
v6_1000_0noiseq = evaluate_DE("NOISeq_6_1000_0_meta_low_log_filtered_Results.txt", "6_1000_0_meta.tsv", "NOISeq")
v9_500_500noiseq = evaluate_DE("NOISeq_9_500_500_meta_low_log_filtered_Results.txt", "9_500_500_meta.tsv", "NOISeq")
v9_750_250noiseq = evaluate_DE("NOISeq_9_750_250_meta_low_log_filtered_Results.txt", "9_750_250_meta.tsv", "NOISeq")
v9_1000_0noiseq = evaluate_DE("NOISeq_9_1000_0_meta_low_log_filtered_Results.txt", "9_1000_0_meta.tsv", "NOISeq")
####################################################################################
#limma log low filtered
v3_500_500limma = evaluate_DE("limmaAG_3_500_500_meta_low_log_filtered_Results.txt", "3_500_500_meta.tsv", "limmaAG")
v3_750_250limma = evaluate_DE("limmaAG_3_750_250_meta_low_log_filtered_Results.txt", "3_750_250_meta.tsv", "limmaAG")
v3_1000_0limma = evaluate_DE("limmaAG_3_1000_0_meta_low_filtered_Results.txt", "3_1000_0_meta.tsv", "limmaAG")
v6_500_500limma = evaluate_DE("limmaAG_6_500_500_meta_low_filtered_Results.txt", "6_500_500_meta.tsv", "limmaAG")
v6_750_250limma = evaluate_DE("limmaAG_6_750_250_meta_low_filtered_Results.txt", "6_750_250_meta.tsv", "limmaAG")
v6_1000_0limma = evaluate_DE("limmaAG_6_1000_0_meta_low_filtered_Results.txt", "6_1000_0_meta.tsv", "limmaAG")
v9_500_500limma = evaluate_DE("limmaAG_9_500_500_meta_low_filtered_Results.txt", "9_500_500_meta.tsv", "limmaAG")
v9_750_250limma = evaluate_DE("limmaAG_9_750_250_meta_low_filtered_Results.txt", "9_750_250_meta.tsv", "limmaAG")
v9_1000_0limma = evaluate_DE("limmaAG_9_1000_0_meta_low_filtered_Results.txt", "9_1000_0_meta.tsv", "limmaAG")
#edgeR low filtered
v3_500_500edgeR = evaluate_DE("edgeR_3_500_500_meta_low_log_filtered_Results.txt", "3_500_500_meta.tsv", "edgeR")
v3_750_250edgeR = evaluate_DE("edgeR_3_750_250_meta_low_log_filtered_Results.txt", "3_750_250_meta.tsv", "edgeR")
v3_1000_0edgeR = evaluate_DE("edgeR_3_1000_0_meta_low_log_filtered_Results.txt", "3_1000_0_meta.tsv", "edgeR")
v6_500_500edgeR = evaluate_DE("edgeR_6_500_500_meta_low_log_filtered_Results.txt", "6_500_500_meta.tsv", "edgeR")
v6_750_250edgeR = evaluate_DE("edgeR_6_750_250_meta_low_log_filtered_Results.txt", "6_750_250_meta.tsv", "edgeR")
v6_1000_0edgeR = evaluate_DE("edgeR_6_1000_0_meta_low_log_filtered_Results.txt", "6_1000_0_meta.tsv", "edgeR")
v9_500_500edgeR = evaluate_DE("edgeR_9_500_500_meta_low_log_filtered_Results.txt", "9_500_500_meta.tsv", "edgeR")
v9_750_250edgeR = evaluate_DE("edgeR_9_750_250_meta_low_log_filtered_Results.txt", "9_750_250_meta.tsv", "edgeR")
v9_1000_0edgeR = evaluate_DE("edgeR_9_1000_0_meta_low_log_filtered_Results.txt", "9_1000_0_meta.tsv", "edgeR")
#########################################################################################
#create a graph with the results
DESeq_results_list <- list(v3_500_500DESeq,
                     v3_750_250DESeq,
                     v3_1000_0DESeq,
                     v6_500_500DESeq,
                     v6_750_250DESeq,
                     v6_1000_0DESeq,
                     v9_500_500DESeq,
                     v9_750_250DESeq,
                     v9_1000_0DESeq
)
edgeR_results_list <- list(v3_500_500edgeR,
                           v3_750_250edgeR,
                           v3_1000_0edgeR,
                           v6_500_500edgeR,
                           v6_750_250edgeR,
                           v6_1000_0edgeR,
                           v9_500_500edgeR,
                           v9_750_250edgeR,
                           v9_1000_0edgeR
)
limma_results_list <- list(v3_500_500limma,
                           v3_750_250limma,
                           v3_1000_0limma,
                           v6_500_500limma,
                           v6_750_250limma,
                           v6_1000_0limma,
                           v9_500_500limma,
                           v9_750_250limma,
                           v9_1000_0limma
)
noiseq_results_list <- list(v3_500_500noiseq,
                           v3_750_250noiseq,
                           v3_1000_0noiseq,
                           v6_500_500noiseq,
                           v6_750_250noiseq,
                           v6_1000_0noiseq,
                           v9_500_500noiseq,
                           v9_750_250noiseq,
                           v9_1000_0noiseq
)
#extract the necessary values, create a dataframe
specificity_values <- c()
sensitivity_values <- c()

for (result in DESeq_results_list) {
  specificity_values <- c(specificity_values, result[["specificity:"]])
  sensitivity_values <- c(sensitivity_values, result[["Sensitivity:"]])
 
}
for (result in edgeR_results_list) {
  specificity_values <- c(specificity_values, result[["specificity:"]])
  sensitivity_values <- c(sensitivity_values, result[["Sensitivity:"]])
}
for (result in limma_results_list) {
  specificity_values <- c(specificity_values, result[["specificity:"]])
  sensitivity_values <- c(sensitivity_values, result[["Sensitivity:"]])
}
for (result in noiseq_results_list) {
  specificity_values <- c(specificity_values, result[["specificity:"]])
  sensitivity_values <- c(sensitivity_values, result[["Sensitivity:"]])
}
df <-data.frame(
  x= specificity_values,
  y= sensitivity_values,
  group1 = rep(c(500,750,1000), times=3),
  group2 = rep(c(3,6,9), each=3,times=1),
  group3 = rep(c("DESeq2", "edgeR", "limma", "NOIseq"),each=9,times=1)
)
print(df)
# Create a column for symbols based on group3
df$symbol <- factor(df$group3, levels = c("DESeq2", "edgeR", "limma", "NOIseq"), labels = c("DESeq2", "edgeR", "limma", "NOIseq"))
##############################################################################################################################################
#figure overall
average_data <- df %>%
  group_by(group3) %>%
  summarize(avg_specificity = mean(x)*100, avg_sensitivity = mean(y)*100)
print(average_data)
# Create the plot
custom_colors <- c("edgeR" = "#7CAE00", "DESeq2" = "#F8766D", "NOIseq" = "#00BFC4","limma" = "#C77CFF" )
ggplot(data = average_data, aes(x = avg_specificity, y = avg_sensitivity, color = group3)) +
  geom_point(size = 3) +
  labs(title = "Overall Performance: Filtered by  FoldChange",
       x = "Specificity(%)",
       y = "Sensitivity(%)",
       color = NULL) +
  scale_color_manual(values = custom_colors)+
  theme_minimal()

##############################################################################################################################################
scatterplot <- ggplot(df, aes(x = x, y = y, shape = symbol)) +
  geom_point() +
  labs(title = "DEG Pipeline Comparations",
       x = "Specificity",
       y="Sensitivity") +
  
  scale_shape_manual(values = c("DESeq2" = 1, "edgeR" = 2, "limma"=3, "NOISeq"=4))  # Customize symbols

# Create a 3 by 3 grid using facet_grid
grid_plot <- scatterplot +
  facet_grid(group1 ~ group2)
################################################
result <- df %>%
  group_by(group1, group3) %>%
  summarise(mean_sensitivity = mean(y))
result
#version 2 of effect of ratio
specificity_plot <- ggplot(df, aes(x = group1, y = x, color = group3)) +
  stat_summary(fun.y = "mean", geom = "point", size = 3) +
  labs(
    title = "Effect of Ratio of Increased DEG on Specificity",
    x = "Number of Overexpressed Genes",
    y = "Specificity",
    color = "Tool"
  ) +
  theme_minimal()

# Create a scatter plot for Sensitivity
sensitivity_plot <- ggplot(df, aes(x = group1, y = y, color = group3)) +
  stat_summary(fun.y = "mean", geom = "point", size = 3) +
  labs(
    title = "Effect of Ratio of Increased DEG on Sensitivity",
    x = "Number of Overexpressed Genes",
    y = "Sensitivity",
    color = "Tool"
  ) +
  theme_minimal()

# Display the Specificity and Sensitivity plots side by side
grid.arrange(specificity_plot, sensitivity_plot, ncol = 2)
################################################
#table 1
###################################################
# Display the grid with the legend outside the graphs
print(grid_plot + theme(legend.position = "bottom"))
#create a table with the averaged values

# Create a matrix from the vectors
DESeq_result_matrix <- t(do.call(rbind, DESeq_results_list))
DESeq_row_averages <- rowMeans(DESeq_result_matrix)
DESeq_result_table <- data.frame(Parameters = names(DESeq_row_averages), Average = DESeq_row_averages)
DESeq_result_table$Average <- sprintf("%.2f", DESeq_result_table$Average)

edgeR_result_matrix <- t(do.call(rbind, edgeR_results_list))
edgeR_row_averages <- rowMeans(edgeR_result_matrix)
edgeR_result_table <- data.frame(Parameters = names(edgeR_row_averages), Average = edgeR_row_averages)
edgeR_result_table$Average <- sprintf("%.2f", edgeR_result_table$Average)

# Create a matrix from the vectors
limma_result_matrix <- t(do.call(rbind, limma_results_list))
limma_row_averages <- rowMeans(limma_result_matrix)
limma_result_table <- data.frame(Parameters = names(limma_row_averages), Average = limma_row_averages)
limma_result_table$Average <- sprintf("%.2f", limma_result_table$Average)

# Create a matrix from the vectors
noiseq_result_matrix <- t(do.call(rbind, noiseq_results_list))
noiseq_row_averages <- rowMeans(noiseq_result_matrix)
noiseq_result_table <- data.frame(Parameters = names(noiseq_row_averages), Average = noiseq_row_averages)
noiseq_result_table$Average <- sprintf("%.2f", noiseq_result_table$Average)



##########################################################
#figure 1
#############################################################
#select the appropriate values
specificity_DESeq<-c()
sensitivity_DESeq<-c()

specificity_edgeR<-c()
sensitivity_edgeR<-c()

specificity_limma<-c()
sensitivity_limma<-c()

specificity_noiseq<-c()
sensitivity_noiseq<-c()

for (values in DESeq_results_list){
  specificity_DESeq <- c(specificity_DESeq, values[["True Negative Rate:"]])
  sensitivity_DESeq <- c(sensitivity_DESeq, values[["Recall:"]])
}

for (values in edgeR_results_list){
  specificity_edgeR <- c(specificity_edgeR, values[["True Negative Rate:"]])
  sensitivity_edgeR <- c(sensitivity_edgeR, values[["Recall:"]])
}

for (result in limma_results_list) {
  specificity_limma <- c(specificity_limma, result[["True Negative Rate:"]])
  sensitivity_limma <- c(sensitivity_limma, result[["Recall:"]])
}
for (result in noiseq_results_list) {
  specificity_noiseq <- c(specificity_noiseq, result[["True Negative Rate:"]])
  sensitivity_noiseq <- c(sensitivity_noiseq, result[["Recall:"]])
}
#average the values by number of samples
result_table_edgeR <- data.frame(
  Group = c("n=3", "n=6", "n=9"),
  specificity = c(mean(specificity_edgeR[1:3]), mean(specificity_edgeR[4:6]), mean(specificity_edgeR[7:9])),
  sensitivity = c(mean(sensitivity_edgeR[1:3]), mean(sensitivity_edgeR[4:6]), mean(sensitivity_edgeR[7:9]))
)


result_table_DESeq <- data.frame(
  Group = c("n=3", "n=6", "n=9"),
  specificity = c(mean(specificity_DESeq[1:3]), mean(specificity_DESeq[4:6]), mean(specificity_DESeq[7:9])),
  sensitivity = c(mean(sensitivity_DESeq[1:3]), mean(sensitivity_DESeq[4:6]), mean(sensitivity_DESeq[7:9]))
)

result_table_limma <- data.frame(
  Group = c("n=3", "n=6", "n=9"),
  specificity = c(mean(specificity_limma[1:3]), mean(specificity_limma[4:6]), mean(specificity_limma[7:9])),
  sensitivity = c(mean(sensitivity_limma[1:3]), mean(sensitivity_limma[4:6]), mean(sensitivity_limma[7:9]))
)

result_table_noiseq <- data.frame(
  Group = c("n=3", "n=6", "n=9"),
  specificity = c(mean(specificity_noiseq[1:3]), mean(specificity_noiseq[4:6]), mean(specificity_noiseq[7:9])),
  sensitivity = c(mean(sensitivity_noiseq[1:3]), mean(sensitivity_noiseq[4:6]), mean(sensitivity_noiseq[7:9]))
)
print(result_table_noiseq)
#create a figure with 2 correlation plot

library(patchwork)
##correlation for DESeq2
correlation_plot_DESeq2 <- ggplot(result_table_DESeq, aes(x = specificity, y = sensitivity, color = Group)) +
  geom_point(size = 3) +
  geom_text(
    aes(label = paste("r =", round(cor(specificity, sensitivity), 2))),
    x = Inf, y = -Inf,
    hjust = 1, vjust = 0,
    size = 4
  ) +
  labs(
   
    title = "DESeq2"
  ) +
  theme_minimal()+
  guides(color = FALSE)
# Display the plot

##correlation for edgeR
correlation_plot_edgeR <- ggplot(result_table_edgeR, aes(x = specificity, y = sensitivity, color = Group)) +
  geom_point(size = 3) +
  geom_text(
    aes(label = paste("r =", round(cor(specificity, sensitivity), 2))),
    x = Inf, y = -Inf,
    hjust = 1, vjust = 0,
    size = 4
  ) +
  labs(
 
    title = "edgeR"
  ) +
  theme_minimal()+
  guides(color = FALSE)

##correlation for limma
correlation_plot_limma <- ggplot(result_table_limma, aes(x = specificity, y = sensitivity, color = Group)) +
  geom_point(size = 3) +
  geom_text(
    aes(label = paste("r =", round(cor(specificity, sensitivity), 2))),
    x = Inf, y = -Inf,
    hjust = 1, vjust = 0,
    size = 4
  ) +
  labs(

    title = "limma"
  ) +
  theme_minimal()+
  guides(color = FALSE)

##correlation for noiseq
correlation_plot_noiseq <- ggplot(result_table_noiseq, aes(x = specificity, y = sensitivity, color = Group)) +
  geom_point(size = 3) +
  geom_text(
    aes(label = paste("r =", round(cor(specificity, sensitivity), 2))),
    x = Inf, y = -Inf,
    hjust = 1, vjust = 0,
    size = 4
  ) +
  labs(
 
    title = "noiseq"
  ) +
  theme_minimal()

combined_plots <- correlation_plot_DESeq2 + correlation_plot_edgeR + correlation_plot_limma+correlation_plot_noiseq+ plot_layout(ncol = 4)

# Display the combined plot
print(combined_plots)
#############################################################################################
#create a plot for specificity over number of samples and one with sensitivity over samples
combined_specificity_plot <- ggplot() +
  geom_point(data = result_table_edgeR, aes(x = Group, y = specificity, color = "edgeR"), size = 3) +
  geom_point(data = result_table_DESeq, aes(x = Group, y = specificity, color = "DESeq2"), size = 3) +
  geom_point(data = result_table_limma, aes(x = Group, y = specificity, color = "limma"), size = 3) +
  geom_point(data = result_table_noiseq, aes(x = Group, y = specificity, color = "noiseq"), size = 3) +
  labs(title = "Effect of Sample Size on Specificity", x = "Samples", y = "Specificity",  color = "Tool") +
  scale_color_manual(values = c("edgeR" = "darkturquoise", "DESeq2" = "coral1","limma" = "black", "noiseq" = "orange")) +
  theme_minimal()

# Create a scatter plot for sensitivity with both DESeq2 and edgeR data
combined_sensitivity_plot <- ggplot() +
  geom_point(data = result_table_edgeR, aes(x = Group, y = sensitivity, color = "edgeR"), size = 3) +
  geom_point(data = result_table_DESeq, aes(x = Group, y = sensitivity, color = "DESeq2"), size = 3) +
  geom_point(data = result_table_limma, aes(x = Group, y = sensitivity, color = "limma"), size = 3) +
  geom_point(data = result_table_noiseq, aes(x = Group, y = sensitivity, color = "noiseq"), size = 3) +
  labs(title = "Effect of Sample Size on Sensitivity", x = "Samples", y = "Sensitivity", color = "Tool") +
  scale_color_manual(values = c("edgeR" = "darkturquoise", "DESeq2" = "coral1", "limma" = "black", "noiseq" = "orange")) +
  theme_minimal()

combined_plots <- combined_specificity_plot + combined_sensitivity_plot + plot_layout(ncol = 4)
print(combined_plots)
#############################################################################################
#look at the pairwise aggrement between tools
result_file_edgeR_6_500_500_filtered <- file.path("Results/edgeR/edgeR_6_500_500_meta_low_log_filtered_Results.txt")
result_file_deseq2_6_500_500_filtered <- file.path("Results/DESeq2/DESeq2_6_500_500_meta_low_log_filtered_Results.txt")
result_file_limma_6_500_500_filtered <- file.path("Results/limmaAG/limmaAG_6_500_500_meta_low_log_filtered_Results.txt")
result_file_noiseq_6_500_500_filtered <- file.path("Results/NOISeq/NOISeq_6_500_500_meta_low_log_filtered_Results.txt")

result_file_edgeR_6_500_500 <- file.path("Results/edgeR/edgeR_6_500_500_meta_low_filtered_Results.txt")
result_file_deseq2_6_500_500 <- file.path("Results/DESeq2/DESeq2_6_500_500_meta_low_filtered_Results.txt")
result_file_limma_6_500_500<- file.path("Results/limmaAG/limmaAG_6_500_500_meta_low_filtered_Results.txt")
result_file_noiseq_6_500_500 <- file.path("Results/NOISeq/NOISeq_6_500_500_meta_low_filtered_Results.txt")
meta_file_6_500_500 <- file.path("Meta/6_500_500_meta.tsv")
#############################
result_file_edgeR_3_500_500_filtered <- file.path("Results/edgeR/edgeR_3_500_500_meta_low_log_filtered_Results.txt")
result_file_deseq2_3_500_500_filtered <- file.path("Results/DESeq2/DESeq2_3_500_500_meta_low_log_filtered_Results.txt")
result_file_limma_3_500_500_filtered <- file.path("Results/limmaAG/limmaAG_3_500_500_meta_low_log_filtered_Results.txt")
result_file_noiseq_3_500_500_filtered <- file.path("Results/NOISeq/NOISeq_3_500_500_meta_low_log_filtered_Results.txt")

result_file_edgeR_3_500_500 <- file.path("Results/edgeR/edgeR_3_500_500_meta_low_filtered_Results.txt")
result_file_deseq2_3_500_500 <- file.path("Results/DESeq2/DESeq2_3_500_500_meta_low_filtered_Results.txt")
result_file_limma_3_500_500 <- file.path("Results/limmaAG/limmaAG_3_500_500_meta_low_filtered_Results.txt")
result_file_noiseq_3_500_500 <- file.path("Results/NOISeq/NOISeq_3_500_500_meta_low_filtered_Results.txt")
###########################################
result_file_edgeR_9_500_500_filtered <- file.path("Results/edgeR/edgeR_9_500_500_meta_low_log_filtered_Results.txt")
result_file_deseq2_9_500_500_filtered <- file.path("Results/DESeq2/DESeq2_9_500_500_meta_low_log_filtered_Results.txt")
result_file_limma_9_500_500_filtered <- file.path("Results/limmaAG/limmaAG_9_500_500_meta_low_log_filtered_Results.txt")
result_file_noiseq_9_500_500_filtered <- file.path("Results/NOISeq/NOISeq_9_500_500_meta_low_log_filtered_Results.txt")

result_file_edgeR_9_500_500<- file.path("Results/edgeR/edgeR_9_500_500_meta_low_filtered_Results.txt")
result_file_deseq2_9_500_500 <- file.path("Results/DESeq2/DESeq2_9_500_500_meta_low_filtered_Results.txt")
result_file_limma_9_500_500 <- file.path("Results/limmaAG/limmaAG_9_500_500_meta_low_filtered_Results.txt")
result_file_noiseq_9_500_500 <- file.path("Results/NOISeq/NOISeq_9_500_500_meta_low_filtered_Results.txt")
###########################################
result_edgeR_9_500_500_filtered <- read.table(result_file_edgeR_9_500_500_filtered, header = TRUE, row.names = 1)
result_deseq2_9_500_500_filtered <- read.table(result_file_deseq2_9_500_500_filtered, header = TRUE, row.names = 1)
result_limma_9_500_500_filtered <- read.table(result_file_limma_9_500_500_filtered, header = TRUE, row.names = 1)
result_noiseq_9_500_500_filtered <- read.table(result_file_noiseq_9_500_500_filtered, header = TRUE, row.names = 1)

result_edgeR_9_500_500<- read.table(result_file_edgeR_9_500_500, header = TRUE, row.names = 1)
result_deseq2_9_500_500 <- read.table(result_file_deseq2_9_500_500, header = TRUE, row.names = 1)
result_limma_9_500_500 <- read.table(result_file_limma_9_500_500, header = TRUE, row.names = 1)
result_noiseq_9_500_500 <- read.table(result_file_noiseq_9_500_500, header = TRUE, row.names = 1)
###########################################
result_edgeR_3_500_500_filtered <- read.table(result_file_edgeR_3_500_500_filtered, header = TRUE, row.names = 1)
result_deseq2_3_500_500_filtered <- read.table(result_file_deseq2_3_500_500_filtered, header = TRUE, row.names = 1)
result_limma_3_500_500_filtered <- read.table(result_file_limma_3_500_500_filtered, header = TRUE, row.names = 1)
result_noiseq_3_500_500_filtered <- read.table(result_file_noiseq_3_500_500_filtered, header = TRUE, row.names = 1)

result_edgeR_3_500_500<- read.table(result_file_edgeR_3_500_500, header = TRUE, row.names = 1)
result_deseq2_3_500_500 <- read.table(result_file_deseq2_3_500_500, header = TRUE, row.names = 1)
result_limma_3_500_500 <- read.table(result_file_limma_3_500_500, header = TRUE, row.names = 1)
result_noiseq_3_500_500 <- read.table(result_file_noiseq_3_500_500, header = TRUE, row.names = 1)
###########################################



##################################################################################################################
#create summarary input 4 df with the result, one for each tool, and output 1 df
create_summary_result_nsamples <- function(edge, deseq2, limma, Noiseq)
{
  number_row <-10000
  tool_columns <- c("edgeR", "DESeq2", "limma", "NOIseq")
  gene_tool_matrix <- data.frame(matrix(0, nrow = number_row, ncol = 4))
  colnames(gene_tool_matrix) <- tool_columns
  # Set the row names to "g1" to "g10000"
  row.names(gene_tool_matrix) <- paste("g", 1:number_row, sep = "")
  
  for (i in 1:number_row) {
    gene_name <- paste("g", i, sep = "")
    gene_tool_matrix[gene_name, "edgeR"] <- ifelse( edge[gene_name, "differential.expression"] == 1, 1, 0)
    gene_tool_matrix[gene_name, "DESeq2"] <- ifelse(deseq2[gene_name, "differential.expression"] == 1, 1, 0)
    gene_tool_matrix[gene_name, "limma"] <- ifelse( limma[gene_name, "differential.expression"] == 1, 1, 0)
    gene_tool_matrix[gene_name, "NOIseq"] <- ifelse( Noiseq[gene_name, "differential.expression"] == 1, 1, 0)
  }
  return(gene_tool_matrix)
  
}
##################################################
#create a df for N3,6,9, filtered by log and low

n_9_result_df_filtered = create_summary_result_nsamples(result_edgeR_9_500_500_filtered, result_deseq2_9_500_500_filtered,result_limma_9_500_500_filtered,result_noiseq_9_500_500_filtered)
n_3_result_df_filtered = create_summary_result_nsamples(result_edgeR_3_500_500_filtered, result_deseq2_3_500_500_filtered,result_limma_3_500_500_filtered,result_noiseq_3_500_500_filtered)
n_9_result_df = create_summary_result_nsamples(result_edgeR_9_500_500, result_deseq2_9_500_500,result_limma_9_500_500,result_noiseq_9_500_500)
n_3_result_df = create_summary_result_nsamples(result_edgeR_3_500_500, result_deseq2_3_500_500,result_limma_3_500_500,result_noiseq_3_500_500)
######################################################

count_deseq2<- sum(gene_tool_matrix[,"DESeq2"])
count_edgeR<-sum(gene_tool_matrix[,"edgeR"])
count_limma <-sum(gene_tool_matrix[,"limma"])
count_noiseq <-sum(gene_tool_matrix[,"NOIseq"])

both_deseq2_edgeR <- rowSums(gene_tool_matrix[, c("DESeq2", "edgeR")]) >0 
both_deseq2_lima <- rowSums(gene_tool_matrix[, c("DESeq2", "limma")]) >0 
both_deseq2_noiseq <- rowSums(gene_tool_matrix[, c("DESeq2", "NOIseq")]) >0 
both_lima_noiseq <- rowSums(gene_tool_matrix[, c("limma", "NOIseq")]) >0 
both_lima_edgeR <- rowSums(gene_tool_matrix[, c("limma", "edgeR")]) >0 
both_edgeR_noiseq <- rowSums(gene_tool_matrix[, c("NOIseq", "edgeR")]) >0 

both_deseq2_edgeR_limma <- rowSums(gene_tool_matrix[, c("DESeq2", "edgeR","limma")]) >0 
both_deseq2_edgeR_noiseq <- rowSums(gene_tool_matrix[, c("DESeq2", "edgeR","NOIseq")]) >0 
both_limma_edgeR_noiseq <- rowSums(gene_tool_matrix[, c("limma", "edgeR","NOIseq")]) >0 
all<- rowSums(gene_tool_matrix[, c("DESeq2", "edgeR","NOIseq","limma")]) >0 

count_both_deseq2_edgeR <- sum(both_deseq2_edgeR)
count_both_deseq2_lima <- sum(both_deseq2_lima)
count_both_deseq2_noiseq <- sum(both_deseq2_noiseq)
count_both_lima_noiseq <- sum(both_lima_noiseq)
count_both_lima_edgeR <- sum(both_lima_edgeR)
count_both_edgeR_noiseq <- sum(both_edgeR_noiseq)
count_both_deseq2_edgeR_limma<-sum(both_deseq2_edgeR_limma)
count_both_deseq2_edgeR_noiseq<-sum(both_deseq2_edgeR_noiseq )
count_both_limma_edgeR_noiseq<-sum(both_limma_edgeR_noiseq)
count_all <-sum(all)

counts <- data.frame(
  Method = c("DESeq2", "edgeR", "limma", "NOIseq","DESeq2+edgeR","DESeq2+limma","DESeq2+NOIseq","limma+NOIseq","limma+edgeR","NOIseq+edgeR","DESeq2+edgeR+limma", "DESeq2+edgeR+NOIseq","limma+edgeR+NOIseq","all"),
  Count = c(count_deseq2, count_edgeR, count_limma, count_noiseq,count_both_deseq2_edgeR,count_both_deseq2_lima, count_both_deseq2_noiseq,count_both_lima_noiseq, count_both_lima_edgeR,count_both_edgeR_noiseq,count_both_deseq2_edgeR_limma,count_both_deseq2_edgeR_noiseq,count_both_deseq2_edgeR_limma,count_all)
)
counts$Method <- reorder(counts$Method, -counts$Count)
print(counts)

library(reshape2)
custom_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", 
  "#1a9850", "#d73027", "#66a61e", "#e6ab02", "#a6761d"
)

df<-melt(counts)
ggplot(counts, aes(x = Method, y = Count, fill = Method)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Counts of Combinations, n=3",
    x = "Method Combination",
    y = "Count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = custom_colors) 

#####################################################################################################
#now look at concordance
##venn Diagram
create_venn<- function(result4tools)
{
  edgeR_genes <- rownames(result4tools[result4tools$edgeR == 1, ],)
  Deseq2_genes <- rownames(result4tools[result4tools$DESeq2 == 1, ])
  limma_genes <- rownames(result4tools[result4tools$limma == 1, ])
  rows_with_na <- grepl("NA", limma_genes)
  # Use the logical vector to subset and remove the row names with "NA"
  limma_genes <- limma_genes[!rows_with_na]
 
  Noiseq_genes <- rownames(result4tools[result4tools$NOIseq == 1, ])
  rows_with_na <- grepl("NA", Noiseq_genes)
  # Use the logical vector to subset and remove the row names with "NA"
  Noiseq_genes <- Noiseq_genes[!rows_with_na]
  X<-list(edgeR_genes,Deseq2_genes,limma_genes,Noiseq_genes)
  
  ggVennDiagram(X,label = "count",label_size = 2,set_size = 3,label_alpha = 0,color = 1, lwd = 0.7,category.names = c("edgeR","DESeq2","limma", "NOISeq"))+
    theme(legend.position = "none")+
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
    theme(plot.margin = unit(c(2,2,2,2), "lines") )+
    scale_x_continuous(expand = expansion(mult = .2))
    
}
venn_n_3 = create_venn(n_3_result_df)
venn_n_3_filtered = create_venn(n_3_result_df_filtered)
venn_n_9 = create_venn(n_9_result_df)
venn_n_9_filtered = create_venn(n_9_result_df_filtered)

combined_plots <- ggarrange(venn_n_3, venn_n_9, venn_n_3_filtered ,venn_n_9_filtered + rremove("x.text"), 
                            labels = c("3 Samples, No LFC Filtering", "9 Samples, No LFC Change Filtering", "3 Samples, Filtered Logfold 0.58","9 Samples, Filtered Logfold 0.58"),
                            ncol = 2, nrow = 2,vjust=2.5, hjust=-0.1)
combined_plots
ggsave("figures/vennFInalD.svg", plot = combined_plots, device = "svg", dpi = 300)
#####################################################################################################
#plot roc curves we want 1 for each
library(pROC)
#define object to plot
rocobjedgeR <- roc(result_edgeR$differential.expression, meta$differential.expression)
rocobjdeseq2<- roc(result_deseq2$differential.expression, meta$differential.expression)

#create ROC plot
ggroc(list(EdgeR = rocobjedgeR, DESeq2 = rocobjdeseq2))

