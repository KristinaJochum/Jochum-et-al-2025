#Cell viability testing with WST-1 and Neutral Red uptake assay
#import test over control values per replicate
#WST1
RPTEC_WST1_72h <- as.data.frame(read_csv2("RPTEC_WST1_72h.csv", trim_ws = TRUE))
rownames(RPTEC_WST1_72h) <- RPTEC_WST1_72h[, 1]  ## set rownames
RPTEC_WST1_72h <- RPTEC_WST1_72h[, -1] 

#NRU
RPTEC_NRU_72h <- as.data.frame(read_csv2("RPTEC_NRU_72h.csv", trim_ws = TRUE))
rownames(RPTEC_NRU_72h) <- RPTEC_NRU_72h[, 1]  ## set rownames
RPTEC_NRU_72h <- RPTEC_NRU_72h[, -1] 


#create dataframe with Means and SD, substances only
##substances only (rownumber)
library(ggplot2)

Cytotox_mean_SD <- function(dataframe){
  needed_rows <- numeric()
  substances <- c("PFOS", "PFOA", "AB1", "Cadmium", "Lasiocarpine")
  n_above_3 <- which(rowSums(!is.na(dataframe[, c(seq(4, ncol(dataframe), 2))])) >= 3)
  print(names(n_above_3))
  for(i in substances){
    row_addition <- n_above_3[grep(pattern = i, names(n_above_3))]
    needed_rows <- c(needed_rows, row_addition)
  }
  substance_means <- rowMeans(dataframe[needed_rows,c(seq(4, ncol(dataframe), 2))], na.rm = T)
  substance_SD <- apply(X = dataframe[needed_rows,c(seq(4, ncol(dataframe), 2))], MARGIN = 1, FUN = sd, na.rm = T)
  means_SD_allrep_dataframe <- data.frame(
    substance = factor(dataframe$substance[needed_rows], 
                       levels = c("PFOS", "PFOA", "AB1", "Cadmium", "Lasiocarpine"),
                       labels = c("PFOS", "PFOA", "Aflatoxin B1", "Cadmium chloride", "Lasiocarpine")), 
    concentration = dataframe$concentration[needed_rows], 
    unit = dataframe$unit[needed_rows], 
    substance_means,
    substance_SD,
    row.names = row.names(dataframe[needed_rows,])
  )
  return(means_SD_allrep_dataframe)
}

means_SD_WST1 <- Cytotox_mean_SD(RPTEC_WST1_72h)
means_SD_NRU <- Cytotox_mean_SD(RPTEC_NRU_72h)

#create basic bar plot
##define data, concentration as factor
df_base_WST1 <- ggplot(
  data = means_SD_WST1,
  aes(x = factor(concentration), y = ifelse(means_SD_WST1$substance_means < 0, 0, means_SD_WST1$substance_means))
)

##specify tings
### errrorbar width and everything, horizontal line, bar width, facetwrap, theme
RPTEC_Cytotox_WST1 <- df_base_WST1 + 
  geom_errorbar(aes(
    ymin = ifelse(means_SD_WST1$substance_means -
      means_SD_WST1$substance_SD < 0, 0, means_SD_WST1$substance_means -
        means_SD_WST1$substance_SD) ,
    ymax = ifelse(means_SD_WST1$substance_means +
      means_SD_WST1$substance_SD < 0, 0, means_SD_WST1$substance_means +
        means_SD_WST1$substance_SD)
  ), width = 0.2, linewidth = 0.3)+
  geom_hline(aes(yintercept = 100), linetype = 2, size = 0.5, linewidth = 0.3)+
  geom_bar(stat = "identity", color = "black", width = 0.5, linewidth = 0.3)+
  facet_wrap(~substance, scales = "free_x", ncol = 6)+
  labs(title = "a", x = "concentration (µM)", y = "cell viability T/C (%)")+
  theme_bw()+
  theme(axis.text = element_text(size = 6),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        plot.title.position = "plot",
        strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 215), breaks = c(seq(50, 200, 50)))

RPTEC_Cytotox_WST1

#create basic bar plot
##define data, concentration as factor
df_base_NRU <- ggplot(
  data = means_SD_NRU,
  aes(x = factor(concentration), y = ifelse(means_SD_NRU$substance_means < 0, 0, means_SD_NRU$substance_means))
)

##specify tings
### errrorbar width and everything, horizontal line, bar width, facetwrap, theme
RPTEC_Cytotox_NRU <- df_base_NRU + 
  geom_errorbar(aes(
    ymin = ifelse(means_SD_NRU$substance_means -
                    means_SD_NRU$substance_SD < 0, 0, means_SD_NRU$substance_means -
                    means_SD_NRU$substance_SD) ,
    ymax = ifelse(means_SD_NRU$substance_means +
                    means_SD_NRU$substance_SD < 0, 0, means_SD_NRU$substance_means +
                    means_SD_NRU$substance_SD)
  ), width = 0.2, linewidth = 0.3)+
  geom_hline(aes(yintercept = 100), linetype = 2, linewidth = 0.3)+
  geom_bar(stat = "identity", color = "black", width = 0.5, linewidth = 0.3)+
  facet_wrap(~substance, scales = "free_x", ncol = 6)+
  labs(title = "b", x = "concentration (µM)", y = "cell viability T/C (%)")+
  theme_bw()+
  theme(axis.text = element_text(size = 6),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        plot.title.position = "plot",
        strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 10))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 215), breaks = c(seq(50, 200, 50)))

RPTEC_Cytotox_NRU

library(gridExtra)

RPTEC_Cytotox <- grid.arrange(RPTEC_Cytotox_WST1, RPTEC_Cytotox_NRU, ncol = 1)

ggsave(
  filename = paste0("Cytotox_Sup", ".png"),
  plot = RPTEC_Cytotox,
  device = "png",
  width = 220,
  height = 130,
  units = "mm",
  dpi = 900)



