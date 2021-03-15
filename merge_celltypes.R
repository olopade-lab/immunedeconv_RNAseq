library(immunedeconv)
library(dplyr)
library(textshape)

args <- commandArgs(trailingOnly = TRUE)
method <- args[1]
if (method == "mcpcounter") {
  method <- "mcp_counter"
} else if (method == "cibersortx") {
  method <- "cibersort"
}
inputFile <- read.csv(args[2])
outputFilePrefix <- args[3]
cellTypes = c("B cell", "T cell CD4+", "T cell CD8+", "T cell",
              "NK cell", "Macrophage", "Neutrophil", "Myeloid dendritic cell")
mappedResults <- inputFile %>% map_result_to_celltypes(cellTypes, method)
mappedResults["cell_type"] <- rownames(mappedResults)
mappedResults <- mappedResults %>%
  select(cell_type, everything())
write.csv(mappedResults,
          paste(outputFilePrefix, "_merged.csv", sep = ""),
          row.names = FALSE)