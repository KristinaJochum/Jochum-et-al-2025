#Array
#get data from GeneGlobe

library(readr)

setwd("3.2 Changes in gene transcription levels/")
HepaRG_IPA_Upload <- read_delim("HepaRG_IPA_Upload.csv", 
                               delim = ";", escape_double = FALSE, locale = locale(), 
                               trim_ws = TRUE)

RPTEC_IPA_Upload <- read_delim("D:/BfR/Labor/Results/RPTEC/Array/.txt_data/Contaminats/RPTEC_IPA_Upload.csv", 
                               delim = ";", escape_double = FALSE, locale = locale(), 
                               trim_ws = TRUE)

HepaRG_Genes_Molecular_Toxicology_Pathway_Finder <- read_csv("HepaRG_Genes_Molecular_Toxicology_Pathway_Finder.csv")

RPTEC_Array_GeneList <- read_csv("D:/BfR/Labor/Results/RPTEC/Array/.txt_data/Contaminats/RPTEC_Array_GeneList.csv")

## perform fdr for all genes

file_1 <- HepaRG_IPA_Upload
file_2 <- RPTEC_IPA_Upload

#run once with HepaRG file and once with RPTEC
#for HepaRG
IPA_upload <- file_1[1:370,]
IPA_upload$genes <- HepaRG_Genes_Molecular_Toxicology_Pathway_Finder$Symbol[1:370]

#for RPTEC
IPA_upload <- file_2[1:84,]
IPA_upload$genes <- RPTEC_Array_GeneList$`Gene Symbol`[1:84]

#get names and number of DEG for no fdr and pvalue below 0.05; perform p-value adjustment
no_fdr <- numeric(0)
no_fdr_names <- numeric(0)
for (i in c(3,5,7,9,11)) {
  fdr <- p.adjust(IPA_upload[[i]], "fdr")
  IPA_upload <- cbind(IPA_upload, fdr)
  no_fdr_names <- c(no_fdr_names, paste(IPA_upload$genes[which(IPA_upload[[i]]<= 0.05& (IPA_upload[[i-1]] <= -1.5 | IPA_upload[[i-1]] >= 1.5))], collapse = ", "))
  no_fdr <- c(no_fdr, length(which(IPA_upload[[i]]<= 0.05 & (IPA_upload[[i-1]] <= -1.5 | IPA_upload[[i-1]] >= 1.5))))
}

no_fdr_names

names(IPA_upload)[13:17] <- c("PFOS_adj_pvalue", "PFOA_adj_pvalue", "Aflatoxin B1_adj_pvalue", "Cadmium chloride_adj_pvalue", "Lasiocarpine_adj_pvalue")

new_order <- c("Refseq", "genes",
               "PFOS_Fold_Regulation", "PFOS_pvalue", "PFOS_adj_pvalue",
               "PFOA_Fold_Regulation", "PFOA_pvalue", "PFOA_adj_pvalue",
               "Aflatoxin B1_Fold_Regulation", "Aflatoxin B1_pvalue", "Aflatoxin B1_adj_pvalue",
               "Cadmium chloride_Fold_Regulation", "Cadmium chloride_pvalue", "Cadmium chloride_adj_pvalue",
               "Lasiocarpine_Fold_Regulation", "Lasiocarpine_pvalue", "Lasiocarpine_adj_pvalue")

IPA_upload <- IPA_upload[, new_order]

#names and number of DEGs for different fdr thresholds
fdr_05 <- numeric(0)
fdr_1 <- numeric(0)
fdr_2 <- numeric(0)

fdr_05_names <- numeric(0)
fdr_1_names <- numeric(0)
fdr_2_names <- numeric(0)
for (i in c(seq(5, 17, 3))) {
  fdr_05 <- c(fdr_05, length(which(IPA_upload[[i]]<= 0.05 & (IPA_upload[[i-2]] <= -1.5 | IPA_upload[[i-2]] >= 1.5))))
  fdr_05_names <- c(fdr_05_names, paste(IPA_upload$genes[which(IPA_upload[[i]]<= 0.05 & (IPA_upload[[i-2]] <= -1.5 | IPA_upload[[i-2]] >= 1.5))], collapse = ", "))
  fdr_1 <- c(fdr_1, length(which(IPA_upload[[i]]<= 0.1 & (IPA_upload[[i-2]] <= -1.5 | IPA_upload[[i-2]] >= 1.5))))
  fdr_1_names <- c(fdr_1_names, paste(IPA_upload$genes[which(IPA_upload[[i]]<= 0.1 & (IPA_upload[[i-2]] <= -1.5 | IPA_upload[[i-2]] >= 1.5))], collapse = ", "))
  fdr_2 <- c(fdr_2, length(which(IPA_upload[[i]]<= 0.2 & (IPA_upload[[i-2]] <= -1.5 | IPA_upload[[i-2]] >= 1.5))))
  fdr_2_names <- c(fdr_2_names, paste(IPA_upload$genes[which(IPA_upload[[i]]<= 0.2 & (IPA_upload[[i-2]] <= -1.5 | IPA_upload[[i-2]] >= 1.5))], collapse = ", "))
}

DEGs <- data.frame(rbind(fdr_05, fdr_1, fdr_2, no_fdr))
DEGs_names <- data.frame(rbind(fdr_05_names, fdr_1_names, fdr_2_names, no_fdr_names))

names(DEGs) <- c("PFOS", "PFOA", "Aflatoxin B1", "Cadmium chloride", "Lasiocarpine")
names(DEGs_names) <- c("PFOS", "PFOA", "Aflatoxin B1", "Cadmium chloride", "Lasiocarpine")
DEGs
#write.csv2(DEGs, "DEG_FDR_HepaRG.csv")

DEGs_names
#write.csv2(DEGs_names, "DEG_names_FDR_RPTEC.csv")

##Heatmap with complexHeatmap
###HepaRG 

library(reshape)

HepaRG_IPA_Upload <- IPA_upload
names(HepaRG_IPA_Upload)[2] <- "Gene"

#load list assorting genes to effects
HepaRG_Array_genes_per_effect <- read_delim("HepaRG_genes_per_effect_MolToxPathFind.csv", 
                                            delim = ";", 
                                            escape_double = FALSE, 
                                            trim_ws = TRUE)

#correct spelling
HepaRG_IPA_Upload$Gene <- gsub("HLA.DRB1", "HLA-DRB1", HepaRG_IPA_Upload$Gene)

#merge to file per effects to get results per effect
HepaRG_IPA_Upload_effect <- merge(
  x = HepaRG_Array_genes_per_effect, 
  y = HepaRG_IPA_Upload, 
  by.x = "Gene", 
  by.y = "Gene", 
  all = T)

#order alphabetically by effect and gene
HepaRG_IPA_Upload_effect <- HepaRG_IPA_Upload_effect[order(
  HepaRG_IPA_Upload_effect$effect, 
  HepaRG_IPA_Upload_effect$Gene, 
  decreasing = c(F,F), 
  method = "radix"),]

row.names(HepaRG_IPA_Upload_effect) <- paste(
  HepaRG_IPA_Upload_effect$effect, 
  HepaRG_IPA_Upload_effect$Gene, 
  sep = "_")

#safe file, this will be adapted a lot now
safe <- HepaRG_IPA_Upload_effect
HepaRG_IPA_Upload_effect <- safe

library(ComplexHeatmap)

#inlcude levels for fold regulation
for(i in grep("_Fold_Regulation", names(HepaRG_IPA_Upload_effect), ignore.case = T)){
  levels <- cut(HepaRG_IPA_Upload_effect[,i], 
                breaks = c(-Inf, -2, -1.5, 1.50, 2.00, Inf), 
                labels = c("< 0.50", "0.50-0.67", "0.67-1.50", "1.50-2.00", "> 2.00"),
                right = FALSE)
  HepaRG_IPA_Upload_effect <-  cbind(HepaRG_IPA_Upload_effect, levels)
  names(HepaRG_IPA_Upload_effect)[ncol(HepaRG_IPA_Upload_effect)] <- paste0(names(HepaRG_IPA_Upload_effect)[i], "_levels")
}

#change adj_pvalue to adjpvalue in names
names(HepaRG_IPA_Upload_effect) <- gsub("adj_p", "adjp", names(HepaRG_IPA_Upload_effect))

#asterisks for pvalues below 0.05
for(i in grep("_pvalue", names(HepaRG_IPA_Upload_effect), ignore.case = T)){
  asterisks <- cut(HepaRG_IPA_Upload_effect[,i], 
                   breaks = c(0, 0.05, Inf), 
                   labels = c("*", ""),
                   right = FALSE)
  HepaRG_IPA_Upload_effect <-  cbind(HepaRG_IPA_Upload_effect, asterisks)
  names(HepaRG_IPA_Upload_effect)[ncol(HepaRG_IPA_Upload_effect)] <- paste0(names(HepaRG_IPA_Upload_effect)[i], "_asterisks")
}

#circles for adj pvalues below 0.05
for(i in grep("adjpvalue", names(HepaRG_IPA_Upload_effect), ignore.case = T)){
  asterisks <- cut(HepaRG_IPA_Upload_effect[,i], 
                breaks = c(0, 0.05, Inf), 
                labels = c("°", ""),
                right = FALSE)
  HepaRG_IPA_Upload_effect <-  cbind(HepaRG_IPA_Upload_effect, asterisks)
  names(HepaRG_IPA_Upload_effect)[ncol(HepaRG_IPA_Upload_effect)] <- paste0(names(HepaRG_IPA_Upload_effect)[i], "_asterisks")
}

#combine asterisks and circles to be displayed in heatmap
for (i in grep("_pvalue_asterisks", names(HepaRG_IPA_Upload_effect), ignore.case = T)) {
  col_display <- paste0(HepaRG_IPA_Upload_effect[,i], HepaRG_IPA_Upload_effect[,i+5])
  HepaRG_IPA_Upload_effect <-  cbind(HepaRG_IPA_Upload_effect, col_display)
  names(HepaRG_IPA_Upload_effect)[ncol(HepaRG_IPA_Upload_effect)] <- paste0(names(HepaRG_IPA_Upload_effect)[i], "_col_display")
}

#set colors for heatmap based on fold regulation levels
colors <- structure(
  c("dodgerblue2", "lightskyblue1", "white", "lightpink", "red2"), 
  names = c("< 0.50", "0.50-0.67", "0.67-1.50", "1.50-2.00", "> 2.00"))

#adapt row and col titles 
HepaRG_IPA_Upload_effect$effect <-  gsub("%", "\n", HepaRG_IPA_Upload_effect$effect)
HepaRG_IPA_Upload_effect <- HepaRG_IPA_Upload_effect[c(1:21, 23,22, 24:26, 28, 27, 29:31, 33, 32, 34:36, 38, 37)]
col_title <- c("PFOS 120 µM", "PFOA 120 µM", expression(AB[1]*" 0.25 µM"), "Las 30 µM", expression(CdCl[2]*" 2.5 µM"))

#Heatmap: 3 single heatmaps will be put together in the end
HT_1 <- Heatmap(HepaRG_IPA_Upload_effect[c(1:137), c(grep("_levels", names(HepaRG_IPA_Upload_effect), ignore.case = T))],
                width = unit(20, "mm"),
                show_heatmap_legend = F,
                name = "fold change",
                col = colors,
                column_names_rot = 90,
                row_labels = HepaRG_IPA_Upload_effect$Gene[c(1:137)],
                row_names_gp = gpar(fontsize = 6),
                split = HepaRG_IPA_Upload_effect$effect[c(1:137)],
                row_title_rot = 360,
                row_title_gp = gpar(fontsize = 7, just = "center"),
                column_labels = col_title,
                column_names_centered = F,
                column_names_gp = gpar(fontsize = 7),
                cluster_columns = F, 
                cluster_rows = F, 
                row_names_side = "right",
                row_title_side = "right",
                border = T,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize =7)),
                cell_fun = function(j, i, x, y, width, height, fill){
                  grid.text(HepaRG_IPA_Upload_effect[i, j+33], x, y - height*0.3, gp= gpar(fontsize = 10))
                }
                )

HT_1

HT_2 <- Heatmap(HepaRG_IPA_Upload_effect[c(138:267), c(grep("_levels", names(HepaRG_IPA_Upload_effect), ignore.case = T))],
                width = unit(20, "mm"),
                name = "fold change",
                col = colors,
                column_names_rot = 90,
                row_labels = HepaRG_IPA_Upload_effect$Gene[c(138:267)],
                row_names_gp = gpar(fontsize = 6),
                row_title_rot = 360,
                row_title_gp = gpar(fontsize = 7, just = "center"),
                column_labels = col_title,
                column_names_centered = F,
                column_names_gp = gpar(fontsize = 7),
                cluster_columns = F, 
                cluster_rows = F, 
                row_names_side = "right", 
                row_title_side = "right",
                split = HepaRG_IPA_Upload_effect$effect[c(138:267)],
                border = T,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize =7)),
                cell_fun = function(j, i, x, y, width, height, fill){
                  grid.text(HepaRG_IPA_Upload_effect[i+137,j+33], x, y - height*0.3, gp= gpar(fontsize = 10))
                },
                show_heatmap_legend = F)

HT_2

HT_3 <- Heatmap(HepaRG_IPA_Upload_effect[c(268:398), c(grep("_levels", names(HepaRG_IPA_Upload_effect), ignore.case = T))],
                width = unit(20, "mm"),
                name = "fold change",
                col = colors,
                column_names_rot = 90,
                row_labels = HepaRG_IPA_Upload_effect$Gene[c(268:398)],
                row_names_gp = gpar(fontsize = 6),
                row_title_rot = 360,
                row_title_gp = gpar(fontsize = 7, just = "center"),
                column_labels = col_title,
                column_names_centered = F,
                column_names_gp = gpar(fontsize = 7),
                cluster_columns = F, 
                cluster_rows = F, 
                row_names_side = "right", 
                row_title_side = "right",
                split = HepaRG_IPA_Upload_effect$effect[c(268:398)],
                border = T,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize =7)),
                cell_fun = function(j, i, x, y, width, height, fill){
                  grid.text(HepaRG_IPA_Upload_effect[i+267,j+33], x, y - height*0.3, gp= gpar(fontsize = 10))
                })

HT_3

#grid of heatmaps

png(filename = "HepaRG_Array.png", 
    width = 220,
    height = 280,
    units = "mm",
    res = 900)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 3,
                                           widths = unit(c(0.9, 0.9, 1.2), "null"),
                                           heights = unit(1, "null"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(HT_1, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(HT_2, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
draw(HT_3, newpage = FALSE)
upViewport()

dev.off()

#percentage of DEGs per effect

df <- HepaRG_IPA_Upload_effect[, order(names(HepaRG_IPA_Upload_effect))]
effects <- unique(HepaRG_IPA_Upload_effect$effect)
new_df <- data.frame(0)

#make bar graph
#count times per effect when FR>1.5 for up or FR<-1.5 for down + significant

for(k in grep("_Fold_regulation$", names(df), ignore.case = T)){
  affected <- numeric(0)
  for(m in 1:nrow(df)){
    if(df[m,k] > 1.5 & df[m, k+3] == "*"){
      affected <- c(affected, df$effect[m])
    }
    if(df[m,k] < -1.5 & df[m, k+3] == "*"){
      affected <- c(affected, df$effect[m])
    }
  }
  percent_vector <- numeric(0)
  for (i in effects) {
    Summe <- numeric(0)
    Summe <- length(grep(i, affected))
    total <- nrow(df[c(grep(i, df$effect)),])
    Percent <- round(Summe*100/total)
    percent_vector <- c(percent_vector, Percent)
  }
  new_df <- cbind(new_df, percent_vector)
}

names(new_df)[2:6] <- names(df[grep("_Fold_regulation$", names(df), ignore.case = T)])
row.names(new_df) <- effects
new_df$effects <- effects

new_df[,1] <- c(rep("HepaRG", nrow(new_df)))
names(new_df)[1] <- "Cell Line"
percentage_HepaRG <- new_df[, c(7,1,6,5,2,4,3)]

#write_excel_csv2(x = new_df, file = "percentage_HepaRG_IPA_1_5.csv")

#Supplement
#Ct values per replicate
HepaRG_Array_rep1 <- as.data.frame(read_delim("HepaRG_Array_rep1_data.txt", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE))

HepaRG_Array_rep2 <- as.data.frame(read_delim("HepaRG_Array_rep2_data.txt", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE))

HepaRG_Array_rep3 <- as.data.frame(read_delim("HepaRG_Array_rep3_data.txt", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE))

HepaRG_Array_rep4 <- as.data.frame(read_delim("HepaRG_Array_rep4_data.txt", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE))

#set cut-off to 35 -> all NA values and values above 35 transformed to 35
cutoff <- function(dataframe_column){
  for(i in 1:length(dataframe_column)){
    if(is.na(dataframe_column[i]) == TRUE){
      dataframe_column[i] <- 35
    }
    if(dataframe_column[i] > 35){
      dataframe_column[i] <- 35
    }
  }
  return(dataframe_column[c(1:length(dataframe_column))])
}

HepaRG_Array_rep1 <- as.data.frame(apply(HepaRG_Array_rep1, 2, cutoff))
HepaRG_Array_rep2 <- as.data.frame(apply(HepaRG_Array_rep2, 2, cutoff))
HepaRG_Array_rep3 <- as.data.frame(apply(HepaRG_Array_rep3, 2, cutoff))
HepaRG_Array_rep4 <- as.data.frame(apply(HepaRG_Array_rep4, 2, cutoff))

#sheet including gene names in same order as on array

HepaRG_Array_Genes <- HepaRG_Genes_Molecular_Toxicology_Pathway_Finder

#assign gene names as rownames

row.names(HepaRG_Array_rep1) <- make.names(HepaRG_Array_Genes$Symbol, unique = T)
row.names(HepaRG_Array_rep2) <- make.names(HepaRG_Array_Genes$Symbol, unique = T)
row.names(HepaRG_Array_rep3) <- make.names(HepaRG_Array_Genes$Symbol, unique = T)
row.names(HepaRG_Array_rep4) <- make.names(HepaRG_Array_Genes$Symbol, unique = T)

#function (applicable to all reps, per rep) for -ddct
##HKG mean NC: mean of HKG of solvent control
##dCt NC: substract HKG_mean_NC from all genes of interest
##HKG mean Substances: per substance: calculate mean of HKGs
##dCt_substance_gene: subtract mean of HKGs from genes of interest
##ddCt_substance: subtract deltaCt of gene of interest of solvent control from deltaCt of gene of interest of test substance incubation

calculate_ddCt <- function(df_Ct){
  HKG_mean_NC <- mean(df_Ct[c(371:375),1])
  dCt_NC <- numeric(0)
  for(k in 1:nrow(df_Ct)){
    dCt_NC_gene <- df_Ct[k,1] - HKG_mean_NC
    dCt_NC <- c(dCt_NC, dCt_NC_gene)
  }
  df_ddCt <- numeric(0)
  for(i in 2:ncol(df_Ct)){
    HKG_mean_substance <- mean(df_Ct[c(371:375),i])
    ddCt_substance <- numeric(0)
    for(m in 1:370){
      dCt_substance_gene <- df_Ct[m,i] - HKG_mean_substance
      ddCt_substance_gene <- dCt_substance_gene - dCt_NC[m]
      ddCt_substance <- c(ddCt_substance, ddCt_substance_gene)
    }
    df_ddCt <- as.data.frame(cbind(df_ddCt, ddCt_substance))
    colnames(df_ddCt)[i-1] <- colnames(df_Ct)[i]
  }
  row.names(df_ddCt) <- row.names(df_Ct)[1:370]
  return(df_ddCt)
}

#apply to all rep
HepaRG_Array_ddCt_rep1 <- calculate_ddCt(HepaRG_Array_rep1)
HepaRG_Array_ddCt_rep2 <- calculate_ddCt(HepaRG_Array_rep2)
HepaRG_Array_ddCt_rep3 <- calculate_ddCt(HepaRG_Array_rep3)
HepaRG_Array_ddCt_rep4 <- calculate_ddCt(HepaRG_Array_rep4)

#change colnames to indicate replicate number
colnames(HepaRG_Array_ddCt_rep1) <- paste(colnames(HepaRG_Array_ddCt_rep1), "rep1", sep = "_")
colnames(HepaRG_Array_ddCt_rep2) <- paste(colnames(HepaRG_Array_ddCt_rep2), "rep2", sep = "_")
colnames(HepaRG_Array_ddCt_rep3) <- paste(colnames(HepaRG_Array_ddCt_rep3), "rep3", sep = "_")
colnames(HepaRG_Array_ddCt_rep4) <- paste(colnames(HepaRG_Array_ddCt_rep4), "rep4", sep = "_")

#bind all reps together
HepaRG_Array_ddCt_allrep <- cbind(HepaRG_Array_ddCt_rep1, HepaRG_Array_ddCt_rep2, HepaRG_Array_ddCt_rep3, HepaRG_Array_ddCt_rep4)

#graph with ggplot
library(ggplot2)

#function applicable to all substances to get one boxplot per substance
boxplot_Array_HepaRG <- function(Substance){
  Sub <- grep(Substance, names(HepaRG_Array_ddCt_allrep))
  datamelt_boxplot <- melt(cbind(Genes = row.names(HepaRG_Array_ddCt_allrep), HepaRG_Array_ddCt_allrep[, Sub]))
  datamelt_boxplot$variable <- gsub("_rep[0-9]+", "", datamelt_boxplot$variable)
  print(unique(datamelt_boxplot$variable))
  datamelt_boxplot$value <- (-1)*datamelt_boxplot$value
  col_title <- c("PFOS 120 µM", "PFOA 120 µM", "Aflatoxin B1 0.25 µM", "Cadmium chloride 2.5 µM", "Lasiocarpine 30 µM")
  datamelt_boxplot$variable <- gsub(Substance, col_title[grep(Substance, col_title)],datamelt_boxplot$variable)
  print(unique(datamelt_boxplot$variable))
  Plot <- ggplot(datamelt_boxplot, aes(Genes, value))+
    geom_point(size = 1)+
    geom_boxplot(width = 0.5, linewidth = 0.3)+
    stat_summary(
      fun.data = "mean_sdl",
      fun.args = list(mult = 1),
      geom = "errorbar",
      width = 0.2,
      linewidth = 0.3
    )+
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size = 1
    )+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)+
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.5)+
    geom_hline(yintercept = -1, linetype = "dashed", color = "red", linewidth = 0.5)+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_blank(),
      strip.text = element_text(size = 6),
      strip.background = element_blank(), 
      strip.text.x = element_blank())+
    ylab(expression("relative gene transcription "*log[2]*"FC"))+
    facet_wrap(~ ifelse(as.numeric(as.factor(Genes)) <= 185, "Oben", "Unten"), scales = "free", ncol = 1) +
    ggtitle(paste("HepaRG", "\u2012", unique(datamelt_boxplot$variable), " "))
  return(Plot)
}

names(HepaRG_Array_ddCt_allrep)
Substances <- c("PFOS", "PFOA", "Aflatoxin.B1", "Cadmium.chloride", "Lasiocarpine")

boxplot_Array_HepaRG(Substances[5])

for(i in Substances){
  p <- boxplot_Array_HepaRG(i)
  ggsave(
    filename = paste0(i,"_Sup_Fig_HepaRG_Array", ".png"),
    plot = p,
    device = "png",
    width = 297,
    height = 210,
    units = "mm",
    dpi = 900)
}

###RPTEC

#get IPA_upload with all adjpvalues from code in the beginning
RPTEC_IPA_Upload <- IPA_upload
names(RPTEC_IPA_Upload)[2] <- "Gene"

#list with genes per effect
RPTEC_Array_genes_per_effect <- read_delim("RPTEC_genes_per_effect_nephotoxicity.csv", 
                                            delim = ";", 
                                            escape_double = FALSE, 
                                            trim_ws = TRUE)
#merge results to genes per effect
RPTEC_IPA_Upload_effect <- merge(
  x = RPTEC_Array_genes_per_effect, 
  y = RPTEC_IPA_Upload, 
  by.x = "Gene", 
  by.y = "Gene", 
  all = T)

#sort alphabetically by effect and gene
RPTEC_IPA_Upload_effect <- RPTEC_IPA_Upload_effect[order(
  RPTEC_IPA_Upload_effect$effect, 
  RPTEC_IPA_Upload_effect$Gene, 
  decreasing = c(F,F), 
  method = "radix"),]

#set rownames
row.names(RPTEC_IPA_Upload_effect) <- paste(
  RPTEC_IPA_Upload_effect$effect, 
  RPTEC_IPA_Upload_effect$Gene, 
  sep = "_")

#safe, this file will be adapted a lot
RPTEC_safe <- RPTEC_IPA_Upload_effect
RPTEC_IPA_Upload_effect <- RPTEC_safe

# graph with complexHeatmap

###inlcude levels for fold regulation in dataframe to be displayed as foldchange in heatmap
for(i in grep("_Fold_Regulation", names(RPTEC_IPA_Upload_effect), ignore.case = T)){
  levels <- cut(RPTEC_IPA_Upload_effect[,i], 
                breaks = c(-Inf, -2, -1.5, 1.50, 2.00, Inf), 
                labels = c("< 0.50", "0.50-0.67", "0.67-1.50", "1.50-2.00", "> 2.00"),
                right = FALSE)
  RPTEC_IPA_Upload_effect <-  cbind(RPTEC_IPA_Upload_effect, levels)
  names(RPTEC_IPA_Upload_effect)[ncol(RPTEC_IPA_Upload_effect)] <- paste0(names(RPTEC_IPA_Upload_effect)[i], "_levels")
}

#change adj_pvalue to adjpvalue
names(RPTEC_IPA_Upload_effect) <- gsub("adj_p", "adjp", names(RPTEC_IPA_Upload_effect))

#include asterisks for pvaalues below 0.05
for(i in grep("_pvalue", names(RPTEC_IPA_Upload_effect), ignore.case = T)){
  asterisks <- cut(RPTEC_IPA_Upload_effect[,i], 
                   breaks = c(0, 0.05, Inf), 
                   labels = c("*", ""),
                   right = FALSE)
  RPTEC_IPA_Upload_effect <-  cbind(RPTEC_IPA_Upload_effect, asterisks)
  names(RPTEC_IPA_Upload_effect)[ncol(RPTEC_IPA_Upload_effect)] <- paste0(names(RPTEC_IPA_Upload_effect)[i], "_asterisks")
}

#include circles for adj pvalues below 0.05
for(i in grep("adjpvalue", names(RPTEC_IPA_Upload_effect), ignore.case = T)){
  asterisks <- cut(RPTEC_IPA_Upload_effect[,i], 
                   breaks = c(0, 0.05, Inf), 
                   labels = c("°", ""),
                   right = FALSE)
  RPTEC_IPA_Upload_effect <-  cbind(RPTEC_IPA_Upload_effect, asterisks)
  names(RPTEC_IPA_Upload_effect)[ncol(RPTEC_IPA_Upload_effect)] <- paste0(names(RPTEC_IPA_Upload_effect)[i], "_asterisks")
}

#combine asterisks and circles to be displayed in heatmap 
for (i in grep("_pvalue_asterisks", names(RPTEC_IPA_Upload_effect), ignore.case = T)) {
  col_display <- paste0(RPTEC_IPA_Upload_effect[,i], RPTEC_IPA_Upload_effect[,i+5])
  RPTEC_IPA_Upload_effect <-  cbind(RPTEC_IPA_Upload_effect, col_display)
  names(RPTEC_IPA_Upload_effect)[ncol(RPTEC_IPA_Upload_effect)] <- paste0(names(RPTEC_IPA_Upload_effect)[i], "_col_display")
}

# set colors for heatmap
colors <- structure(
  c("dodgerblue2", "lightskyblue1", "white", "lightpink", "red2"), 
  names = c("< 0.50", "0.50-0.67", "0.67-1.50", "1.50-2.00", "> 2.00"))

#adapt row and col titles 
RPTEC_IPA_Upload_effect$effect <-  gsub("%", "\n", RPTEC_IPA_Upload_effect$effect)
RPTEC_IPA_Upload_effect <- RPTEC_IPA_Upload_effect[c(1:21, 23,22, 24:26, 28, 27, 29:31, 33, 32, 34:36, 38, 37)]
col_title <- c("PFOS 100 µM", "PFOA 480 µM", expression(AB[1]*" 10 µM"), "Las 10 µM", expression(CdCl[2]*" 10 µM"))

#Heatmap
png(filename = "RPTEC_Array.png", 
    width = 70,
    height = 280,
    units = "mm",
    res = 900)

Heatmap(RPTEC_IPA_Upload_effect[, c(grep("_levels", names(RPTEC_IPA_Upload_effect)))],
        name = "fold change",
        width = unit(20, "mm"),
        col = colors,
        column_names_rot = 90,
        row_labels = RPTEC_IPA_Upload_effect$Gene,
        row_names_gp = gpar(fontsize = 6),
        row_title = gsub(" ", "\n", unique(RPTEC_IPA_Upload_effect$effect)),
        row_title_rot = 360,
        row_title_gp = gpar(fontsize = 7, just = "center"),
        column_labels = col_title,
        column_names_centered = F,
        column_names_gp = gpar(fontsize = 7),
        cluster_columns = F, 
        cluster_rows = F, 
        row_names_side = "right",
        row_title_side = "right",
        split = RPTEC_IPA_Upload_effect$effect,
        border = T,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize =7)),
        cell_fun = function(j, i, x, y, width, height, fill){
          grid.text(RPTEC_IPA_Upload_effect[i,j+33], x, y - height*0.3, gp= gpar(fontsize = 10))
        })

dev.off()

#percentage of DEGs per effect
df <- RPTEC_IPA_Upload_effect[, order(names(RPTEC_IPA_Upload_effect))]
effects <- unique(RPTEC_IPA_Upload_effect$effect)
new_df <- data.frame(0)

#count times per effect when FR>1.5 for up or FC<-1.5 for down + significant
for(k in grep("_Fold_regulation$", names(df), ignore.case = T)){
  affected <- numeric(0)
  for(m in 1:nrow(df)){
    if(df[m,k] > 1.5 & df[m, k+3] == "*"){
      affected <- c(affected, df$effect[m])
    }
    if(df[m,k] < -1.5 & df[m, k+3] == "*"){
      affected <- c(affected, df$effect[m])
    }
  }
  percent_vector <- numeric(0)
  for (i in effects) {
    Summe <- numeric(0)
    Summe <- length(grep(i, affected))
    total <- nrow(df[c(grep(i, df$effect)),])
    Percent <- round(Summe*100/total)
    percent_vector <- c(percent_vector, Percent)
  }
  new_df <- cbind(new_df, percent_vector)
}

names(new_df)[2:6] <- names(df[grep("_Fold_regulation$", names(df), ignore.case = T)])
row.names(new_df) <- effects
new_df$effects <- effects

new_df[,1] <- c(rep("RPTEC", nrow(new_df)))
names(new_df)[1] <- "Cell Line"
percentage_RPTEC <- new_df[,c(7,1,6,5,2,4,3)]

#combine HepaRG and RPTEC
percentage_HepaRG_RPTEC <- rbind(percentage_HepaRG, percentage_RPTEC)
names(percentage_HepaRG_RPTEC) <- gsub("_Fold_Regulation", "", names(percentage_HepaRG_RPTEC))

#write_csv2(x = percentage_HepaRG_RPTEC, file = "percentage_HepaRG_RPTEC_IPA_1_5.csv")

#Supplementary
#import Ct values per effect

RPTEC_Array_rep1 <- as.data.frame(read_delim("RPTEC_Array_rep1_data.txt", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE))

RPTEC_Array_rep2 <- as.data.frame(read_delim("RPTEC_Array_rep2_data.txt", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE))

RPTEC_Array_rep3 <- as.data.frame(read_delim("RPTEC_Array_rep3_data.txt", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE))

#cutoff every Ct above 35
cutoff <- function(dataframe_column){
  for(i in 1:length(dataframe_column)){
    if(is.na(dataframe_column[i]) == TRUE){
      dataframe_column[i] <- 35
    }
    if(dataframe_column[i] > 35){
      dataframe_column[i] <- 35
    }
  }
  return(dataframe_column[c(1:length(dataframe_column))])
}

RPTEC_Array_rep1 <- as.data.frame(apply(RPTEC_Array_rep1, 2, cutoff))
RPTEC_Array_rep2 <- as.data.frame(apply(RPTEC_Array_rep2, 2, cutoff))
RPTEC_Array_rep3 <- as.data.frame(apply(RPTEC_Array_rep3, 2, cutoff))

#assign gene names as rownames
row.names(RPTEC_Array_rep1) <- make.names(RPTEC_Array_GeneList$`Gene Symbol`, unique = T)
row.names(RPTEC_Array_rep2) <- make.names(RPTEC_Array_GeneList$`Gene Symbol`, unique = T)
row.names(RPTEC_Array_rep3) <- make.names(RPTEC_Array_GeneList$`Gene Symbol`, unique = T)

#explanation see HepaRG
calculate_ddCt <- function(df_Ct){
  HKG_mean_NC <- mean(df_Ct[c((nrow(df_Ct) - 11):(nrow(df_Ct) - 7)),1])
  dCt_NC <- numeric(0)
  for(k in 1:nrow(df_Ct)){
    dCt_NC_gene <- df_Ct[k,1] - HKG_mean_NC
    dCt_NC <- c(dCt_NC, dCt_NC_gene)
  }
  df_ddCt <- numeric(0)
  for(i in 2:ncol(df_Ct)){
    HKG_mean_substance <- mean(df_Ct[c((nrow(df_Ct) - 11):(nrow(df_Ct) - 7)),i])
    ddCt_substance <- numeric(0)
    for(m in 1:(nrow(df_Ct) - 12)){
      dCt_substance_gene <- df_Ct[m,i] - HKG_mean_substance
      ddCt_substance_gene <- dCt_substance_gene - dCt_NC[m]
      ddCt_substance <- c(ddCt_substance, ddCt_substance_gene)
    }
    df_ddCt <- as.data.frame(cbind(df_ddCt, ddCt_substance))
    colnames(df_ddCt)[i-1] <- colnames(df_Ct)[i]
  }
  row.names(df_ddCt) <- row.names(df_Ct)[1:(nrow(df_Ct) - 12)]
  return(df_ddCt)
}

RPTEC_Array_ddCt_rep1 <- calculate_ddCt(RPTEC_Array_rep1)
RPTEC_Array_ddCt_rep2 <- calculate_ddCt(RPTEC_Array_rep2)
RPTEC_Array_ddCt_rep3 <- calculate_ddCt(RPTEC_Array_rep3)

colnames(RPTEC_Array_ddCt_rep1) <- paste(colnames(RPTEC_Array_ddCt_rep1), "rep1", sep = "_")
colnames(RPTEC_Array_ddCt_rep2) <- paste(colnames(RPTEC_Array_ddCt_rep2), "rep2", sep = "_")
colnames(RPTEC_Array_ddCt_rep3) <- paste(colnames(RPTEC_Array_ddCt_rep3), "rep3", sep = "_")

RPTEC_Array_ddCt_allrep <- cbind(RPTEC_Array_ddCt_rep1, RPTEC_Array_ddCt_rep2, RPTEC_Array_ddCt_rep3)

#function for boxplot with ggplot2; define 2 substances (plot is for 2 substances, can be NA for only one substance), if you want to save the plot ("YES"/"NO") and the height (plot with only one substance can be shorter)

library(ggplot2)

boxplot_Array_RPTEC <- function(substance1, substance2, save, plot_height){
  names(RPTEC_Array_ddCt_allrep) <- gsub("\\.", " ", names(RPTEC_Array_ddCt_allrep))
  substances1_2 <- c(grep(substance1, names(RPTEC_Array_ddCt_allrep)), grep(substance2, names(RPTEC_Array_ddCt_allrep)))
  datamelt_boxplot <- melt(cbind(Genes = row.names(RPTEC_Array_ddCt_allrep), RPTEC_Array_ddCt_allrep[, substances1_2]))
  datamelt_boxplot$variable <- gsub("_rep[0-9]+", "", datamelt_boxplot$variable)
  print(unique(datamelt_boxplot$variable))
  datamelt_boxplot$value <- (-1)*datamelt_boxplot$value
  col_title <- c("PFOS 100 µM", "PFOA 480 µM", "Aflatoxin B1 10 µM", "Cadmium chloride 10 µM", "Lasiocarpine 10 µM")
  datamelt_boxplot$variable <- ifelse(datamelt_boxplot$variable == substance1, col_title[grep(substance1, col_title)], col_title[grep(substance2, col_title)])
  print(unique(datamelt_boxplot$variable))
  datamelt_boxplot$variable <- factor(datamelt_boxplot$variable, levels = c(col_title[grep(substance1, col_title)], col_title[grep(substance2, col_title)]))
  Plot <- ggplot(datamelt_boxplot, aes(Genes, value))+
    geom_point(size = 1)+
    geom_boxplot(width = 0.5, linewidth = 0.3)+
    stat_summary(
      fun.data = "mean_sdl",
      fun.args = list(mult = 1),
      geom = "errorbar",
      width = 0.2,
      linewidth = 0.3
    )+
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size = 1
    )+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)+
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.5)+
    geom_hline(yintercept = -1, linetype = "dashed", color = "red", linewidth = 0.5)+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_blank(),
      strip.text = element_text(size = 8))+
    ylab(expression("relative gene transcription "*log[2]*"FC"))+
    facet_wrap(~variable, 
               ncol = 1)+
    ggtitle("RPTEC")
  if(save == "YES"){
    print("save")
    print(getwd())
    ggsave(
      filename = paste0(substance1, "_", substance2, ".png"),
      plot = Plot,
      device = "png",
      width = 297,
      height = plot_height,
      units = "mm",
      dpi = 900)
  }
  return(Plot)
}

boxplot_Array_RPTEC("PFOS", "PFOA", "NO", 210)

boxplot_Array_RPTEC("Aflatoxin B1", "Cadmium chloride", "NO", 210)

boxplot_Array_RPTEC("Lasiocarpine", "NA", "NO", 105)
