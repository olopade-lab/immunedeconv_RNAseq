library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(scales)

deconv_Nigerian <- read.csv("output/Nigerian_cibersortx_results_reformatted.csv")
deconv_TCGA <- read.csv("output/TCGA_cibersortx_results_reformatted.csv", check.names=FALSE)

NigerianPAM50 <- read.csv("Nigerian_PAM50.csv")
TCGAPAM50 <- read.csv("TCGA_PAM50.csv")

rownames(deconv_Nigerian) <- gsub(" ", ".", gsub("\\+", "", gsub(" \\(Tregs\\)", "", deconv_Nigerian$cell_type)))
deconv_Nigerian <- subset(deconv_Nigerian, select = -c(cell_type))
deconv_Nigerian <- as.data.frame(t(deconv_Nigerian))
deconv_Nigerian$sample_name <- rownames(deconv_Nigerian)
rownames(deconv_Nigerian) <- 1:nrow(deconv_Nigerian)
deconv_Nigerian$sample_name <- gsub("\\.", "-", deconv_Nigerian$sample_name)
deconv_Nigerian <- merge(deconv_Nigerian, NigerianPAM50)
deconv_Nigerian <- deconv_Nigerian[deconv_Nigerian["PAM50_subtype"] == "Basal", ]

rownames(deconv_TCGA) <- gsub(" ", ".", gsub("\\+", "", gsub(" \\(Tregs\\)", "", deconv_TCGA$cell_type)))
deconv_TCGA <- subset(deconv_TCGA, select = -c(cell_type))
deconv_TCGA <- as.data.frame(t(deconv_TCGA))
deconv_TCGA$sample_name <- rownames(deconv_TCGA)
rownames(deconv_TCGA) <- 1:nrow(deconv_TCGA)
# deconv_TCGA$sample_name <- gsub("\\.", "-", deconv_TCGA$sample_name)
# deconv_TCGA$sample_name <- gsub("X", "", rownames(deconv_TCGA))
deconv_TCGA <- merge(deconv_TCGA, TCGAPAM50)
deconv_TCGA <- deconv_TCGA[deconv_TCGA["PAM50_subtype"] == "Basal", ]

deconv <- rbind(deconv_Nigerian, deconv_TCGA)

total_dataframe = data.frame(sample_name = character(),
                             PAM50_subtype = character(),
                             score = double(),
                             cell_type = character())

for (cell_type in list("B.cell.naive", "T.cell.CD8", "T.cell.follicular.helper",
                       "T.cell.regulatory", "NK.cell.activated", "Mast.cell.activated")) {
  specific_dataframe <- deconv[, c("sample_name", "group", cell_type)]
  colnames(specific_dataframe)[3] <- "score"
  specific_dataframe$cell_type <- cell_type
  total_dataframe <- rbind(total_dataframe, specific_dataframe)
}

total_dataframe$group <- revalue(total_dataframe$group, c("TCGA_black"="TCGA-AA",
                                                    "TCGA_white"="TCGA-EA"))

get_significant_pairs <- function(deconv_dataframe, cell_type, group_pairs) {
  significant_pairs <- list()
  i <- 1
  for (row in 1:nrow(group_pairs)) {
    list_x <- deconv_dataframe[
      deconv_dataframe["group"] == as.character(group_pairs[row, 1]), ][[cell_type]]
    list_y <- deconv_dataframe[
      deconv_dataframe["group"] == as.character(group_pairs[row, 2]), ][[cell_type]]
    if (t.test(list_x, list_y)[["p.value"]] < 0.05) {
      significant_pairs[[i]] <- c(as.character(group_pairs[row, 1]), as.character(group_pairs[row, 2]))
      i <- i + 1
    }
  }
  return(significant_pairs)
}

group_pairs <- data.frame(c("Nigerian", "Nigerian", "TCGA_black"),
                          c("TCGA_black", "TCGA_white", "TCGA_white"))

xlabs <- paste(gsub("TCGA_white", "TCGA-EA", gsub("TCGA_black", "TCGA-AA", levels(deconv$group))),
               "\n(N=",table(deconv$group),")",sep="")

cibersortx_stacked_barplot <- ggplot(total_dataframe, aes(fill=cell_type, y=score, x=sample_name)) + geom_col() +
  geom_bar(position="stack", stat="identity") + facet_grid(~ group, scales = "free", space = "free") + 
  scale_y_continuous(name = "Score", expand = c(0, 0)) + theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  scale_fill_discrete(name="Cell type",
                      breaks=c("B.cell.naive", "T.cell.CD8", "T.cell.follicular.helper",
                              "T.cell.regulatory", "NK.cell.activated", "Mast.cell.activated"),
                      labels=c("B cell naive", "T cell CD8+", "T cell follicular helper",
                               "T cell regulatory", "NK cell activated", "mast cell activated")) +
  labs(x="Samples")

deconv$score <- rowSums(deconv[c("B.cell.naive", "T.cell.CD8", "T.cell.follicular.helper",
                                 "T.cell.regulatory", "NK.cell.activated", "Mast.cell.activated")])

significant_pairs <- get_significant_pairs(deconv, "score", group_pairs)
cibersortx_boxplot <- ggplot(deconv, aes(x = group, y = score, fill = group)) +
  geom_boxplot(show.legend = FALSE) + 
  scale_x_discrete(labels = xlabs, name = "Group") +
  scale_y_continuous(name = "Score") +
  theme_bw() + stat_compare_means(comparisons = significant_pairs, method = "t.test")

both_plots <- ggarrange(cibersortx_stacked_barplot, cibersortx_boxplot, ncol = 1)
ggsave("CIBERSORTx_barplot_boxplot.jpg", both_plots, dpi = 500, width = 10, height = 10)


horizontal_plot <- ggplot(total_dataframe, aes(fill=cell_type, y=score, x=sample_name)) + geom_col() +
  geom_bar(position="stack", stat="identity")  + facet_grid(subtype_SheilaDoc ~ ., scales = "free", space = "free") + coord_flip() +
  theme(axis.title.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), strip.text.y = element_text(angle = 0))

ggplot(total_dataframe, aes(fill=cell_type, y=score, x=sample_name)) + geom_col() +
  geom_bar(position="stack", stat="identity") + facet_grid(subtype_SheilaDoc ~ ., scales = "free", space = "free") + coord_flip() +
  theme(axis.title.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), strip.text.y = element_text(angle = 0))

ggsave("stacked_barplot_horizontal.jpg", horizontal_plot, dpi = 500)