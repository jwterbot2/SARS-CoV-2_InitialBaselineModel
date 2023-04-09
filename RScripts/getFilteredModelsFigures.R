#Loads data from a file listing models to analyze by sourcing loadDataList.R
#Then draws and saves the diversity figure for each model contained within $Data_Storage

for (i in c(1:length(Data_Storage)))
{
  fig = getBranchDiversityFigure_v1(Data_Storage[i], setNames("100", "_partial_100"), setNames("1000", "1000"), 0.02, snp_y_max = 25, tpi_y_max = 0.003);
  ggsave(paste("Figure", names(Data_Storage[i]),"_100.png", sep=""), fig, width = 5, height = 3, dpi = 300)
  
  fig = getBranchDiversityFigure_v1(Data_Storage[i], setNames("1000", "_partial_1000"), setNames("1000", "1000"), 0.02, snp_y_max = 25, tpi_y_max = 0.003);
  ggsave(paste("Figure", names(Data_Storage[i]),"_1000.png", sep=""), fig, width = 5, height = 3, dpi = 300)
}