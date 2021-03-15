library(immunedeconv)

args <- commandArgs(trailingOnly = TRUE)
method <- args[1]
inputFile <- read.csv(args[2])
outputFilePrefix <- args[3]
outputFilename <- paste(outputFilePrefix, "_", method, "_results.csv", sep = "")
rowInputFile <- as.character(inputFile[,1])
inputFile <- inputFile[,2:dim(inputFile)[2]]
row.names(inputFile) <- rowInputFile

if (method == "timer") {
  indicVec <- rep('brca', ncol(inputFile))
  deconvolutionResults <- immunedeconv::deconvolute(inputFile, "timer", indicVec)
} else if (method == "xcell") {
  deconvolutionResults <- immunedeconv::deconvolute(inputFile, "xcell")
} else if (method == "epic") {
  deconvolutionResults <- immunedeconv::deconvolute(inputFile, "epic", tumor = TRUE)
} else if (method == "mcpcounter") {
  deconvolutionResults <- immunedeconv::deconvolute(inputFile, "mcp_counter")
} else if (method == "quantiseq") {
  deconvolutionResults <- immunedeconv::deconvolute(inputFile, "quantiseq",
                                                    tumor = TRUE, scale_mrna = TRUE)
} else {
  print(paste("Method name not recognized. ",
              "Please use one of the following: ",
              "epic, mcpcounter, quantiseq, timer, xcell",
              sep = ""))
  quit()
}
write.csv(deconvolutionResults, outputFilename, row.names = FALSE)

