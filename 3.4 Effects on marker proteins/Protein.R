# Skript to evaluate Multiplexed microsphere-based sandwich immunoassays
##Import data

library(readr)
setwd("3.4 Effects on marker proteins/")

HepaRG_Protein_Result <- read_delim("HepaRG_Protein.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE, col_select = c(1:15))


RPTEC_Protein_Result <- read_delim("RPTEC_Protein.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)


# create dataframe with rep1-3 side by side
##concentration as character for row.name
HepaRG_Protein_Result$concentration_uM <- as.character(HepaRG_Protein_Result$concentration_uM)

RPTEC_Protein_Result$concentration_uM <- as.character(RPTEC_Protein_Result$concentration_uM)

## one df per rep

HepaRG_Protein_Result_rep1 <-  HepaRG_Protein_Result[which(HepaRG_Protein_Result$biol_rep == 1),]
HepaRG_Protein_Result_rep2 <-  HepaRG_Protein_Result[which(HepaRG_Protein_Result$biol_rep == 2),]
HepaRG_Protein_Result_rep3 <-  HepaRG_Protein_Result[which(HepaRG_Protein_Result$biol_rep == 3),]
HepaRG_Protein_Result_rep4 <-  HepaRG_Protein_Result[which(HepaRG_Protein_Result$biol_rep == 4),]

RPTEC_Protein_Result_rep1 <-  RPTEC_Protein_Result[which(RPTEC_Protein_Result$biol_rep == 1),]
RPTEC_Protein_Result_rep2 <-  RPTEC_Protein_Result[which(RPTEC_Protein_Result$biol_rep == 2),]
RPTEC_Protein_Result_rep3 <-  RPTEC_Protein_Result[which(RPTEC_Protein_Result$biol_rep == 3),]
RPTEC_Protein_Result_rep4 <-  RPTEC_Protein_Result[which(RPTEC_Protein_Result$biol_rep == 4),]

## change rownames to apply merge
### function: enter dataframe name
### creates empty vector from rownames
### creates names for column 1-3
### apply make.names to be added to df
### add to df

change_rownames <- function(mydataframe){
  rowname_vector <- numeric(0)
  for(i in 1:nrow(mydataframe)){
    rowname <- paste(mydataframe[i,5], mydataframe[i,7], mydataframe[i,2], sep = "_")
    rowname_vector <- c(rowname_vector, rowname)
  }
  mydataframe <- cbind(rowname_vector, mydataframe)
  make.names(rowname_vector, unique = T)
  .rowNamesDF(x = mydataframe, make.names = T) <- rowname_vector
    return(mydataframe)
}

HepaRG_Protein_Result_rep1 <- change_rownames(HepaRG_Protein_Result_rep1)
HepaRG_Protein_Result_rep2 <- change_rownames(HepaRG_Protein_Result_rep2)
HepaRG_Protein_Result_rep3 <- change_rownames(HepaRG_Protein_Result_rep3)
HepaRG_Protein_Result_rep4 <- change_rownames(HepaRG_Protein_Result_rep4)

RPTEC_Protein_Result_rep1 <- change_rownames(RPTEC_Protein_Result_rep1)
RPTEC_Protein_Result_rep2 <- change_rownames(RPTEC_Protein_Result_rep2)
RPTEC_Protein_Result_rep3 <- change_rownames(RPTEC_Protein_Result_rep3)
RPTEC_Protein_Result_rep4 <- change_rownames(RPTEC_Protein_Result_rep4)

# calulate T/C per rep
##create function with input dataframe of rep
### attach df, create empty vectors
### go through whole df and look for substance_number which is not 0
### if so, go through whole df and look if a row with substance_number 0 matches regarding timepoint and concentration number
### if so, calculate T/C (100%), add to df
### df is list, has to be converted to list and colums converted to numeric
### return new df

calculate_T_C <- function(replicate){
  attach(replicate)
  newtable <- numeric(0)
  newthing <- numeric(0)
  for(i in 1:nrow(replicate)){
    if(substance_number[i] != 0){
      for(k in 1:nrow(replicate)){
        if(time_point[k] == time_point[i] & concentration_number[k] == concentration_number[i] & substance_number[k] == 0){
          newthing <- replicate[i, c(9:16)]/replicate[k,c(9:16)]*100
          rowofnewthing <- c(replicate[i,1:8], newthing)
          newtable <- rbind(newtable, rowofnewthing)
        }
      }
    }
  }
  detach(replicate)
  for(m in c(1,3,6)){
    newtable[,m] <- as.character(newtable[,m])
  }
  newname <- as.data.frame(newtable)
  rownames(newname) <- newname$rowname_vector
    for(l in c(2,4,5,7,8:16)){
    newname[,l] <- as.numeric(newname[,l])
  }
  return(newname)
}

HepaRG_Protein_Result_rep1_T_C <- calculate_T_C(HepaRG_Protein_Result_rep1)
HepaRG_Protein_Result_rep2_T_C <- calculate_T_C(HepaRG_Protein_Result_rep2)
HepaRG_Protein_Result_rep3_T_C <- calculate_T_C(HepaRG_Protein_Result_rep3)
HepaRG_Protein_Result_rep4_T_C <- calculate_T_C(HepaRG_Protein_Result_rep4)

RPTEC_Protein_Result_rep1_T_C <- calculate_T_C(RPTEC_Protein_Result_rep1)
RPTEC_Protein_Result_rep2_T_C <- calculate_T_C(RPTEC_Protein_Result_rep2)
RPTEC_Protein_Result_rep3_T_C <- calculate_T_C(RPTEC_Protein_Result_rep3)
RPTEC_Protein_Result_rep4_T_C <- calculate_T_C(RPTEC_Protein_Result_rep4)

# merge 4 dataframes
## first merge 2 twice, merge these
## change rownames back, as they get set as new colomns
## get 3 rows with indications of the rows (as originally) from the new row.names
## create one dataframe and sort by substance name and concentration
## set rownames again
### input: dataframes, suffixes as "_suffix"

merge_4_data.frames <- function(dataframe1, suffix1, dataframe2, suffix2, dataframe3, suffix3, dataframe4, suffix4){
  firstmerge <- merge(dataframe1[,c(9:16)],dataframe2[,c(9:16)], by = "row.names", all.x = T, all.y = T, suffixes = c(suffix1, suffix2))
  rownames(firstmerge) <- firstmerge$Row.names
  secondmerge <- merge(dataframe3[,c(9:16)],dataframe4[,c(9:16)], by = "row.names", all.x = T, all.y = T, suffixes = c(suffix3, suffix4))
  rownames(secondmerge) <- secondmerge$Row.names
  thirdmerge <- merge(firstmerge[,c(2:ncol(firstmerge))],secondmerge[,c(2:ncol(secondmerge))], by = "row.names", all.x = T, all.y = T)
  rownames(thirdmerge) <- thirdmerge$Row.names
  Beschriftung <- data.frame(t(data.frame(strsplit(as.character(thirdmerge[c(1:length(thirdmerge$Row.names)),1]), split = "_"))))
  dataframename <- cbind(Beschriftung, thirdmerge[,c(2:ncol(thirdmerge))])
  row.names(dataframename) <- thirdmerge$Row.names
  names(dataframename)[1] <- "substance"
  names(dataframename)[2] <- "concentration_µM"
  names(dataframename)[3] <- "timepoint"
  return(dataframename)
}

HepaRG_Protein_Result_merge <- merge_4_data.frames(HepaRG_Protein_Result_rep1_T_C, "_rep1", HepaRG_Protein_Result_rep2_T_C, "_rep2", HepaRG_Protein_Result_rep3_T_C, "_rep3", HepaRG_Protein_Result_rep4_T_C, "_rep4")

RPTEC_Protein_Result_merge <- merge_4_data.frames(RPTEC_Protein_Result_rep1_T_C, "_rep1", RPTEC_Protein_Result_rep2_T_C, "_rep2", RPTEC_Protein_Result_rep3_T_C, "_rep3", RPTEC_Protein_Result_rep4_T_C, "_rep4")

#calculate mean
## create vector for new colnames
## create df with row information
### calculate rowmean for each analyte and attach it to df
### attach colnames
### same process for sd
## apply to merge df 
#statistics
## same process for stat
### extract p-value into a separate vector ---> add to df

library(boot)

calculate_mean_sd <- function(dataframe_for_mean_sd){
  coloumnames <- gsub("_rep[0-9]", "", names(dataframe_for_mean_sd))
  mean_SD_df <- data.frame(dataframe_for_mean_sd[,c(1:3)])
  for (i in 4:11){
    analyte_mean <- rowMeans(dataframe_for_mean_sd[, c(seq(from = i, by = 8, length.out = 4))], na.rm = T)
    mean_SD_df <- cbind(mean_SD_df, analyte_mean)
    names(mean_SD_df)[i] <- paste(coloumnames[i], "mean", sep = "_")
  }
  for (i in 4:11){
    analyte_sd <- apply(dataframe_for_mean_sd[, c(seq(from = i, by = 8, length.out = 4))], MARGIN = 1, FUN = sd, na.rm = T)
    mean_SD_df <- cbind(mean_SD_df, analyte_sd)
    names(mean_SD_df)[i+8] <- paste(coloumnames[i], "sd", sep = "_")
  }
  statistic_function <- function(data, indices){
    sample_data <- data[indices]
    return(mean(sample_data))
  }
  for (i in 4:11) {
    t_test_substance <- numeric(0)
    p_value_as_asterisks <- numeric(0)
    bootstrap_as_asterisks <- numeric(0)
    for(k in 1:nrow(dataframe_for_mean_sd)){
      t_test <- t.test(as.numeric(dataframe_for_mean_sd[k, c(seq(from = i, by = 8, length.out = 4))]), mu = 100)
      t_test_substance <- c(t_test_substance, t_test$p.value)
      if(t_test$p.value < 0.05){
        as_asterisks <- "*"
      }else{
        as_asterisks <- ""
      }
      p_value_as_asterisks <- c(p_value_as_asterisks, as_asterisks)
      boot_result <- boot(na.omit(as.numeric(dataframe_for_mean_sd[k, c(seq(from = i, by = 8, length.out = 4))])), statistic_function, R = 1000)
      boostrap_means <- boot_result$t
      boot.ci(boot_result, type = "basic")
      confidence_interval <- boot.ci(boot_result, type = "basic")$basic
      if(confidence_interval[1,4] < 100 & 100 < confidence_interval[1,5]){
        bootstrap <- ""
      }else{
        bootstrap <- "*"
      }
      bootstrap_as_asterisks <- c(bootstrap_as_asterisks, bootstrap)
    }
    mean_SD_df <- cbind(mean_SD_df, t_test_substance, p_value_as_asterisks, bootstrap_as_asterisks)
    names(mean_SD_df) <- gsub("t_test_substance", paste(names(mean_SD_df[i]), "p-value", sep = "_"), names(mean_SD_df))
    names(mean_SD_df) <- gsub("p_value_as_asterisks", paste(names(mean_SD_df[i]), "asterisks", sep = "_"), names(mean_SD_df))
    names(mean_SD_df) <- gsub("bootstrap_as_asterisks", paste(names(mean_SD_df[i]), "boot", sep = "_"), names(mean_SD_df))
  }
  return(mean_SD_df)
}

HepaRG_Protein_Result_mean_SD <- calculate_mean_sd(HepaRG_Protein_Result_merge)

RPTEC_Protein_Result_mean_SD <- calculate_mean_sd(RPTEC_Protein_Result_merge)

HepaRG_Protein_Result_mean_SD <- cbind(row.names(HepaRG_Protein_Result_mean_SD), HepaRG_Protein_Result_mean_SD)
names(HepaRG_Protein_Result_mean_SD)[1] <- "longname"

RPTEC_Protein_Result_mean_SD <- cbind(row.names(RPTEC_Protein_Result_mean_SD), RPTEC_Protein_Result_mean_SD)
names(RPTEC_Protein_Result_mean_SD)[1] <- "longname"

###Heatmap with Complex Heatmaps

library(ComplexHeatmap)

HepaRG_Protein_Result_mean_SD_CH <- cbind(CellLine = rep("HepaRG", times = nrow(HepaRG_Protein_Result_mean_SD)), HepaRG_Protein_Result_mean_SD[, c(2,3,4,5:12, seq(from = 23, by = 3, to = 44))])
RPTEC_Protein_Result_mean_SD_CH <- cbind(CellLine = rep("RPTEC", times = nrow(RPTEC_Protein_Result_mean_SD)), RPTEC_Protein_Result_mean_SD[, c(2,3,4,5:12, seq(from = 23, by = 3, to = 44))])

make_levels <- function(mean_SD) {
  library(dplyr)
  for(i in grep("_mean$", names(mean_SD), ignore.case = T)){
    levels <- cut(mean_SD[,i], 
                  breaks = c(-Inf, 75, 150, 200, Inf), 
                  labels = c("< 0.75", "0.75-1.50", "1.50-2.00", "> 2.00"),
                  right = FALSE)
    mean_SD <-  cbind(mean_SD, levels)
    names(mean_SD)[ncol(mean_SD)] <- paste0(gsub("mean", "", names(mean_SD)[i]), "_levels")
  }
  mean_SD[, "concentration_µM"] <- as.numeric(mean_SD[, "concentration_µM"])
  ordered_df <- mean_SD %>%
    arrange(substance, concentration_µM)
  return(ordered_df)
}

HepaRG_Protein_Result_mean_SD_CH <- make_levels(HepaRG_Protein_Result_mean_SD_CH)
RPTEC_Protein_Result_mean_SD_CH <-  make_levels(RPTEC_Protein_Result_mean_SD_CH)

###Heatmap parameters

ht_opt("legend_border" = "black")
ht_opt$COLUMN_ANNO_PADDING = unit(0.5, "cm")
ht_opt$ROW_ANNO_PADDING = unit(0.2, "cm")

colors <- structure(
  c("lightskyblue1", "white", "lightpink", "red2"), 
  names = c("< 0.75", "0.75-1.50", "1.50-2.00", "> 2.00"))

row_title <- c("PFOS", "PFOA", expression(AB[1]), "Las", expression(CdCl[2]))

ha_HepaRG = rowAnnotation(
  conc = anno_text(
    paste0(" ", HepaRG_Protein_Result_mean_SD_CH$concentration_µM, " µM"),
    just = "right",
    location = 1,
    gp = gpar(fontsize = 8))
  )

title_HepaRG <-  HeatmapAnnotation(
  title = anno_text(c(rep("", 5), "A", rep("", 2)), rot = 0, 
                    gp = gpar(fontsize = 15))
)

lgd = list(labels = c("< 0.75", "0.75-1.50", "1.50-2.00", "> 2.00"),
           legend_gp = gpar(fill = c("lightskyblue1", "white", "lightpink", "red2")), 
           labels_gp = gpar(fontsize = 8),
           title_gp = gpar(fontsize = 10),
           nrow = 1,
           title_position = "leftcenter")


HT_HepaRG <- Heatmap(HepaRG_Protein_Result_mean_SD_CH[, c(21:28)],
                     row_split = factor(HepaRG_Protein_Result_mean_SD_CH[, "substance"],
                                        levels = c("PFOS", "PFOA", "Aflatoxin B1", "Lasiocarpine", "Cadmium chloride")),
                     row_gap = unit(0.2, "cm"),
                     col = colors,
                     name = "fold change",
                     column_title_gp = gpar(just = "left"),
                     cluster_row_slices = F,
                     cluster_rows = F,
                     cluster_columns = F,
                     row_title_rot = F,
                     row_title = row_title,
                     row_title_gp = gpar(fontsize = 10),
                     row_names_side = "left",
                     row_names_gp = gpar(fontsize = 8),
                     column_labels = gsub("_mean", "", names(HepaRG_Protein_Result_mean_SD_CH[, c(5:12)])),
                     column_names_rot = 45,
                     column_names_gp = gpar(fontsize = 8),
                     column_order = order(gsub("_mean", "", names(HepaRG_Protein_Result_mean_SD_CH[, c(5:12)]))),
                     row_labels = gsub("h", " h", HepaRG_Protein_Result_mean_SD_CH$timepoint),
                     border = T,
                     top_annotation = title_HepaRG,
                     left_annotation = ha_HepaRG,
                     cell_fun = function(j, i, x, y, width, height, fill){
                       grid.text(sprintf("%.1f", (HepaRG_Protein_Result_mean_SD_CH[i, j+4]/100)), x, y, 
                                 gp= gpar(fontsize = 6))
                     },
                     layer_fun = function(j, i, x, y, w, h, f){
                       grid.text(pindex(as.matrix(HepaRG_Protein_Result_mean_SD_CH), i, j+12), x+w*0.35, y+h*0.25, 
                                 gp = gpar(fontsize = 10, fontface = "bold"))
                     },
                     heatmap_legend_param = lgd
)

HT_HepaRG

title_RPTEC = HeatmapAnnotation(
  title = anno_text(c(rep("", 5), "B", rep("", 2)), rot = 0,
                    gp = gpar(fontsize = 15)
                    )
  )

ha_RPTEC = rowAnnotation(
  conc = anno_text(
    paste0(RPTEC_Protein_Result_mean_SD_CH$concentration_µM, " µM"),
    just = "right",
    location = 1,
    gp = gpar(fontsize = 8)
    )
  )

HT_RPTEC <- Heatmap(RPTEC_Protein_Result_mean_SD_CH[, c(21:28)],
                     row_split = factor(HepaRG_Protein_Result_mean_SD_CH[, "substance"],
                                        levels = c("PFOS", "PFOA", "Aflatoxin B1", "Lasiocarpine", "Cadmium chloride")),
                     col = colors,
                     name = "fold change",
                     column_title_gp = gpar(just = "left"),
                     cluster_row_slices = F,
                     cluster_rows = F,
                     cluster_columns = F,
                     row_title_rot = F,
                     row_title = row_title,
                     row_names_side = "left",
                     column_labels = gsub("_mean", "", names(RPTEC_Protein_Result_mean_SD_CH[, c(5:12)])),
                     column_names_rot = 45,
                     column_order = order(gsub("_mean", "", names(RPTEC_Protein_Result_mean_SD_CH[, c(5:12)]))),
                     column_names_gp = gpar(fontsize = 8),
                     row_labels = gsub("h", " h", RPTEC_Protein_Result_mean_SD_CH$timepoint),
                     border = T,
                     left_annotation = ha_RPTEC,
                     top_annotation = title_RPTEC,
                     cell_fun = function(j, i, x, y, width, height, fill){
                       grid.text(sprintf("%.1f", (RPTEC_Protein_Result_mean_SD_CH[i, j+4]/100)), x, y, 
                                 gp= gpar(fontsize = 6))
                     },
                     layer_fun = function(j, i, x, y, w, h, f){
                       grid.text(pindex(as.matrix(RPTEC_Protein_Result_mean_SD_CH), i, j+12), x+w*0.35, y+h*0.25, 
                                 gp = gpar(fontsize = 10, fontface = "bold"))
                     },
                    show_heatmap_legend = F
)

HT_RPTEC

HT_combined <- HT_HepaRG + HT_RPTEC

HT_combined

png("Protein.png", width = 21, height = 15, units = "cm", res = 900)

draw(HT_combined, 
     heatmap_legend_side = "bottom",
     ht_gap = unit(1, "cm"),)

dev.off()

#spupplementary material
##boxplot
###HepaRG

library(reshape)
library(ggplot2)

HepaRG_Protein_datamelt_boxplot <- melt(cbind(longnames = row.names(HepaRG_Protein_Result_merge), HepaRG_Protein_Result_merge))

HepaRG_Protein_datamelt_boxplot$variable <- gsub("_rep[0-9]+", "", HepaRG_Protein_datamelt_boxplot$variable)

Sup_Fig_HepaRG <- ggplot(HepaRG_Protein_datamelt_boxplot[which(is.na(HepaRG_Protein_datamelt_boxplot$value) == F),], aes(longnames, value))+
  geom_point(size = 1)+
  geom_boxplot(width = 0.5, size = 0.3)+
  stat_summary(
    fun.data = "mean_sdl",
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.2,
    size = 0.3
  )+
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 1
  )+
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", linewidth = 0.5)+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_blank())+
  ylab("relative protein expression T/C (%)")+
  facet_wrap(~variable, scales = "free_y", ncol = 4)+
  ggtitle("HepaRG")

ggsave("Sup_Fig_HepaRG_Protein.png",width = 9.21, height = 6.85, dpi = 1200)
 
###RPTEC

RPTEC_Protein_datamelt_boxplot <- melt(cbind(longnames = row.names(RPTEC_Protein_Result_merge), RPTEC_Protein_Result_merge))

RPTEC_Protein_datamelt_boxplot$variable <- gsub("_rep[0-9]+", "", RPTEC_Protein_datamelt_boxplot$variable)

ggplot(RPTEC_Protein_datamelt_boxplot[which(is.na(RPTEC_Protein_datamelt_boxplot$value) == F),], aes(longnames, value))+
  geom_point(size = 1)+
  geom_boxplot(width = 0.5, size = 0.3)+
  stat_summary(
    fun.data = "mean_sdl",
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.2,
    size = 0.3
  )+
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 1
  )+
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", size = 0.5)+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_blank())+
  ylab("relative protein expression T/C (%)")+
  facet_wrap(~variable, scales = "free_y", ncol = 4)+
  ggtitle("RPTEC")

ggsave("Sup_Fig_RPTEC_Protein.png",width = 9.21, height = 6.85, dpi = 1200)


