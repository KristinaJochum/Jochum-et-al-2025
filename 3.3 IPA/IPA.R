#IPA Graph for Paper
#Y-Axis: Diseases (top 10 overall), X-Axis: p-value, Bars: HepaRG, RPTEC, combined

#1.upload all results from HepaRG, RPTEC, combi
#2.combine all of theme in one df by effect
#3.determine the effects coming up the most based on p-value as well
#4.make bar graph as described above

## get library
library(readr)
library(tidyverse)

## data files as txt
### generate list of all files .txt (=datafiles)
### read all files from list, 
### skip first 2 rows, only cat and p-value cols
### remove doubled categories and only keep the one with the lowest p-value
#### Group by category and select the row with the minimum number,  ungroup
### change colname to include cell line, generate name from temp as substance
### assign file with corresponding filename, and add to result_list

setwd("3.3 IPA")

Import <- function(CellLine){
  setwd(paste0(CellLine, "/"))
  temp <- list.files(pattern = "*.txt")
  result_list <- list()
  for (i in 1:length(temp)) {
    file <- read_tsv(temp[i], col_names = T, trim_ws = TRUE, skip = 2, locale = locale(decimal_mark = ","))
    file <- file[,c(1,3)]
    result_df <- file %>%
      group_by(Categories) %>%
      summarise(`p-value` = min(`p-value`))
    result_df <- result_df %>% ungroup()
    colnames(result_df) <- gsub("p-value", paste(colnames(result_df[2]), CellLine, gsub(".txt", "", temp[i]), sep = "_"), colnames(result_df), )
    filename <- paste(gsub(".txt", "", temp[i]))
    result_list[[filename]] <- result_df
  }
  setwd(("../"))
  return(result_list)
}

HepaRG_IPA <- Import("HepaRG")
RPTEC_IPA <- Import("RPTEC")
combined_IPA <- Import("combined")

### combine more than one df
#### generate list of all df

df_list <- list(as.data.frame(HepaRG_IPA[["PFOS"]]),
                as.data.frame(HepaRG_IPA[["PFOA"]]), 
                as.data.frame(HepaRG_IPA[["Aflatoxin B1"]]),
                as.data.frame(HepaRG_IPA[["Cadmium chloride"]]),
                as.data.frame(HepaRG_IPA[["Lasiocarpine"]]),
                as.data.frame(RPTEC_IPA[["PFOS"]]),
                as.data.frame(RPTEC_IPA[["PFOA"]]), 
                as.data.frame(RPTEC_IPA[["Aflatoxin B1"]]),
                as.data.frame(RPTEC_IPA[["Cadmium chloride"]]),
                as.data.frame(combined_IPA[["PFOS"]]),
                as.data.frame(combined_IPA[["PFOA"]]), 
                as.data.frame(combined_IPA[["Aflatoxin B1"]]),
                as.data.frame(combined_IPA[["Cadmium chloride"]]),
                as.data.frame(combined_IPA[["Lasiocarpine"]])
                )

#### merged by Categories

allCat_pvalue <- Reduce(function(x, y) merge(x, y, by = "Categories", all=TRUE), df_list)

###remove non-liver/-kidney effects

Cardiac_effects <- c(
  grep("Heart", allCat_pvalue$Categories),
  grep("Cardiac", allCat_pvalue$Categories),
  grep("Pulmonary", allCat_pvalue$Categories)
)

Cat_pvalue <- allCat_pvalue[-Cardiac_effects, ]
row.names(Cat_pvalue) <- Cat_pvalue$Categories

##applying evaluation criteria

for (i in c(grep("p-value", names(Cat_pvalue)))){
  evaluation <- cut(Cat_pvalue[[i]], 
                    breaks = c(0, 0.0005, 0.005, 0.05, 1), 
                    labels = c("+++", "++", "+", NA),
                    right = FALSE)
  Cat_pvalue <- cbind(Cat_pvalue, evaluation)
  names(Cat_pvalue)[ncol(Cat_pvalue)] <- paste0(names(Cat_pvalue)[i], "_asterix")
}


Cat_pvalue$rowSums <- rowSums(!is.na(Cat_pvalue))

#write.csv2(Cat_pvalue, "Cat_pvalue_evaluation.csv", row.names = F)

top10 <- tail(Cat_pvalue[order(rowSums(!is.na(Cat_pvalue))),-c(grep("_asterix", names(Cat_pvalue)))], 10)

top10$`p-value_RPTEC_Lasiocarpine` <- rep(NA, 10)

#transform for plot
library(reshape)
top10_melt <- melt(top10[, c(1, grep("p-value", names(top10)))])
Beschriftung <- data.frame(t(data.frame(strsplit(as.character(top10_melt[,2]), split = "_"))))
top10_melt <- cbind(Beschriftung[, c(2,3)], top10_melt)

##graph
library(ggplot2)
##use ggpattern to get rid of black color in graph
library(ggpattern)

##transform data
transformed_data <- transform(
  top10_melt,
  X2 = factor(X2, levels = c("HepaRG", "RPTEC", "combined")),
  Categories = factor(Categories, levels = row.names(top10)),
  pattern = ifelse(
    factor(top10_melt$X3, levels = rev(c("PFOS", "PFOA", "Aflatoxin B1", "Lasiocarpine", "Cadmium chloride"))) %in% 
      rev(c("Lasiocarpine", "Cadmium chloride")), "stripe", "none"
  ))

## Create the plot
gg <- ggplot(
  transformed_data,
  aes(
    y = -log10(value),
    x = Categories,
    fill = factor(X3, levels = rev(c("PFOS", "PFOA", "Aflatoxin B1", "Lasiocarpine", "Cadmium chloride"))),
    pattern = pattern
  )
)

## Add geom_col_pattern
gg <- gg + 
  geom_col_pattern(
    width = 0.6,
    position = position_dodge(width = 0.7),
    color = "black",
    lwd = 0.8,
    pattern_colour = "grey34",
    pattern_fill = "grey34",
    pattern_density = 0.4,  # Adjust stripe density
    pattern_spacing = 0.02  # Adjust stripe spacing
  )+
  scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe"))

## Add coordinates and facets
gg <- gg +
  coord_flip() +
  scale_x_discrete(labels = gsub(",", "\n", row.names(top10))) +
  facet_wrap(~X2, ncol = 3, labeller = as_labeller(c("HepaRG" = "HepaRG", "RPTEC" = "RPTEC/TERT1", "combined" = "combined")))

## Add guides and scales
gg <- gg +
  scale_fill_manual(
    values = c("white", "grey74", "grey24", "white", "grey74"),
    breaks = c("PFOS", "PFOA", "Aflatoxin B1", "Lasiocarpine", "Cadmium chloride"),
    labels = c("PFOS", "PFOA", expression(AB[1]), "Las", expression(CdCl[2]))
  )+ theme_bw() 

## Add labels and theme
gg <- gg +
  labs(x = "", y = expression("-"*log[10]*"(p-value)")) +
  guides(
    pattern = "none",
    fill = guide_legend(
      title = "",
      override.aes = list(
        pattern = c("none", "none", "none", "stripe", "stripe"),
        fill = c("white", "grey74", "grey24", "white", "grey74")
      )
    )
  )+
  theme(strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10), 
        axis.title.x = element_text(size = 12), 
        legend.text = element_text(size = 12))

gg

ggsave(
  filename = "IPA.png", device = "png",
  width = 297,
  height = 210,
  units = "mm",
  dpi = 900)

