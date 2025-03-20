##get genes of interest per substance in one file 
###genes of interest are: p < 0.05 (not adjusted) and Fold regulation < -1.5 or >1.5
###graph with how many upregulated genes and downregulated genes per substance
#read GeneGlobe files with fold regulation and p-values

library(readr)
library(tidyverse)
library(ggplot2)

setwd("3.2 Changes in gene transcription levels/")

Results_HepaRG <- read_delim("HepaRG_IPA_Upload.csv", delim = ";", escape_double = FALSE, locale = locale(), trim_ws = TRUE, n_max = 370)
                                                  
Results_RPTEC <- read_delim("RPTEC_IPA_Upload.csv", delim = ";", escape_double = FALSE, locale = locale(), trim_ws = TRUE, n_max = 84)

#genes of interest
#find interesting genes

interesting_genes_function <- function(df_results){
  for(i in seq(2, 11, 2)){
    interesting_genes <- mutate(df_results[, c(i, i+1)],
                                interesting_candidate_fc=ifelse(df_results[i] <= -1.5, "DOWN", ifelse(df_results[i] >= 1.5, "UP", "NO")),
                                interesting_candidate_pval=ifelse(df_results[i+1] >= 0.05,"NO","YES"),
                                interesting_candidate=ifelse(interesting_candidate_fc== "DOWN" & interesting_candidate_pval=="YES", "DOWN", ifelse(interesting_candidate_fc== "UP" & interesting_candidate_pval=="YES","UP","NO")),
                                .keep = "none")
    names(interesting_genes) <- gsub("interesting_candidate", paste0("GO_genes_", gsub("_Fold_Regulation", "", names(df_results[i]))), names(interesting_genes))
    interesting_genes <- apply(interesting_genes, 2, unlist)
    df_results <- data.frame(df_results, as.data.frame(interesting_genes))
  }
  return(df_results)
}

interesting_genes_HepaRG <- interesting_genes_function(Results_HepaRG)
interesting_genes_RPTEC <- interesting_genes_function(Results_RPTEC)

#create Barchart with up and down regulated genes
##function to get number of up and down regulated genes per substance
sum_up_down <- function(interesting_genes){
  substances <- c("PFOS", "PFOA", "AB1", "CdCl2", "Las")
  vector_up <- numeric(0)
  vector_down <- numeric(0)
  for (i in seq(14, 26, 3)){
    up <- length(grep("UP", interesting_genes[,i]))
    vector_up <- c(vector_up, up)
    down <- length(grep("DOWN", interesting_genes[,i]))
    vector_down <- c(vector_down, down)
  } 
  df <- data.frame(substances = substances, up = vector_up, down = vector_down, row.names = substances)
  return(df)
}

sum_up_down_HepaRG <- sum_up_down(interesting_genes_HepaRG)
sum_up_down_RPTEC <- sum_up_down(interesting_genes_RPTEC)

##mutate df for ggplot2
df1_long <- sum_up_down_HepaRG[, c(2,3)] %>%
  mutate(Category = factor(rownames(sum_up_down_HepaRG), levels = c("CdCl2", "Las", "AB1", "PFOA", "PFOS")), 
    Dataset = "HepaRG") %>%
  pivot_longer(cols = -c(Category, Dataset), names_to = "Variable", values_to = "Value")

df2_long <- sum_up_down_RPTEC[, c(2,3)] %>%
  mutate(Category = factor(rownames(sum_up_down_RPTEC), levels = c("CdCl2", "Las", "AB1", "PFOA", "PFOS")), Dataset = "RPTEC") %>%
  pivot_longer(cols = -c(Category, Dataset), names_to = "Variable", values_to = "Value")

combined_df <- bind_rows(df1_long, df2_long)
combined_df$Variable <- factor(combined_df$Variable, c("up", "down"))

##barchart with HepaRG and RPTEC next to each other
ggplot(combined_df, aes(x = Category, y = Value, fill = Variable)) +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.7)) +
  geom_text(aes(label = Value), hjust = -0.3, size = 5, position = position_dodge(width = 0.7))+
  facet_wrap(~Dataset, scales = "free", nrow = 2, strip.position = "right", labeller = as_labeller(c("HepaRG" = "HepaRG \n(total = 370)", "RPTEC" = "RPTEC/TERT1 \n(total = 84)"))) +
  theme_bw() +
  scale_fill_manual(values = c("up" = "red2", "down" = "dodgerblue2"))+
  scale_x_discrete(labels = c(expression(CdCl[2]), "Las", expression(AB[1]), "PFOA", "PFOS"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  coord_flip()+
  labs(x = element_blank(),
       y = "number of DEG",
       fill = "direction") +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        panel.spacing = unit(1, "cm"),
        legend.position = "bottom")

ggsave("number_of_DEG.png", width = unit(12, "cm"), height = unit(7, "cm"))

