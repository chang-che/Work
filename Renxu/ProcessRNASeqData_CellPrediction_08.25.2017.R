

dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_533Samples/'
load(file=paste(dataPathRNASeq,'FPKM533.RData',sep=''))
write.table(fpkm, file = 'FPKM533.txt', sep = "\t", quote = FALSE)
## note: I downloaded the fpkm data from server, and run xCell on the webpage
proportion = read.table("xCell_FINAL_WTCProjectFPKM533_xCell_1830082517.txt",
                        sep = "\t", header = TRUE, row.names = 1)
dim(proportion)
raw_score = read.table("xCell_RAW_WTCProjectFPKM533_xCell_1830082517.txt",
                       sep = "\t", header = TRUE, row.names = 1)
save(proportion, file = "WTCProjectCellProportion.RData")
save(raw_score, file = "WTCProjectRawScore.RData")






















