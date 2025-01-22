setwd("C:/Users/lorie/OneDrive - Queensland University of Technology/Semester_2/IFN646_BiomedicalDataScience/assessments/project/Project")
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(dplyr)
library(ggVennDiagram)
library(dplyr)
library(ggpubr)
########################################################################
#folders
########################################################################
##figure 1 , will plot filtered vs unfiltered for low cout, without log fold
########################################################################

extract_number_deg<-function(result_file)
{
  #taxe a result file and count the number of 1 in DEG, does not compare with gound truth
  
  data<-read.table(paste("Results/",result_file,sep = ""), header=T, row.names=1)
  number_DEG <- sum(data$differential.expression)
  return (number_DEG)
}

tools <- c("edgeR", "DESeq2", "NOIseq","limmaAG")
ratios <- c("500_500", "750_250", "1000_0")
methods <-c("low_filtered","log_filtered","low_log_filtered","no_filtered" )
sample_sizes<-c(3,6,9)

# Create an empty list to store tables for each sample size
result_df <- data.frame(tools=c("edgeR", "DESeq2", "NOIseq","limmaAG","edgeR", "DESeq2", "NOIseq","limmaAG","edgeR", "DESeq2", "NOIseq","limmaAG"))
result_df$sample <- c(3,3,3,3,6,6,6,6,9,9,9,9)
result_df$low_filtered <-c(0,0,0,0,0,0,0,0,0,0,0,0)
result_df$log_filtered <-c(0,0,0,0,0,0,0,0,0,0,0,0)
result_df$low_log_filtered <-c(0,0,0,0,0,0,0,0,0,0,0,0)
result_df$no_filtered <-c(0,0,0,0,0,0,0,0,0,0,0,0)

# Loop through each sample size
for (sample_size in sample_sizes) 
  {
  # Initialize an empty data frame to store results for this sample size
  # Loop through each tool
  for (tool in tools) 
    {
    # Loop through each method
    for (method in methods) 
      {
     
      number_DEG <- 0
      
      for (ratio in ratios) 
        {
        # Create the formatted string
        formatted_string <- paste(tool,paste(tool,sample_size, ratio,"meta", method, "Results.txt", sep = "_"), sep="/")
        # Call the extract_number_deg function and store the result
        number_DEG = number_DEG +extract_number_deg(formatted_string)
        }
      
      result_df[result_df$tools == tool & result_df$sample == sample_size, method] <- number_DEG
    }
  }
}
print(result_df)
###################################################################
#filtered low count vs unfiltered low count

custom_colors <- c("edgeR" = "#7CAE00", "DESeq2" = "#F8766D", "NOIseq" = "#00BFC4","limmaAG" = "#C77CFF" )
filtered_low_plot = ggplot(result_df, aes(x = sample,group=tools, color=tools)) +
  geom_line(aes(y = low_filtered)) +
  geom_line(aes(y = no_filtered),linetype="twodash")+
  scale_x_continuous(breaks = c(3, 6, 9))+# Use geom_line to create a line plot
  labs(
    title = "Effect of Low Count Filtering on the Number of Detected DEG",
    x = "Number of Samples",
    y = "Number of Detected Genes")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),  
        panel.grid.minor = element_blank(),
        legend.title.align = 0.5)+
  scale_color_manual(values = custom_colors)+
  guides(color = guide_legend(title = NULL))
  
print(filtered_low_plot)
ggsave("figures/countFiltering.svg", plot = filtered_low_plot, device = "svg", dpi = 300)
###############################
#filtered logfold count vs unfiltered log
filtered_log_plot = ggplot(result_df, aes(x = sample,group=tools, color=tools)) +
  geom_line(aes(y = log_filtered)) +
  geom_line(aes(y = no_filtered),linetype="twodash")+
  scale_x_continuous(breaks = c(3, 6, 9)) +# Use geom_line to create a line plot
  labs(
    title = "Effect of logFoldchange Filtering on the Number of Detected DEG. lfc=0.58",
    x = "Number of Samples",
    y = "Number of Detected Genes")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        legend.title.align = 0.5)+
  scale_color_manual(values = custom_colors)+
  guides(color = guide_legend(title = NULL)) 

print(filtered_log_plot)
ggsave("figures/FOldFiltering.svg", plot = filtered_log_plot, device = "svg", dpi = 300)


