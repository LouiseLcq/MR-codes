######F statistics######
setwd("D:/postgraduate.research/my subject/Mendelian randomization/outcome data Finngen to CIF/TNK/10scatter/harmonise")
file_list <- list.files()
for (file in file_list) {
  data <- read.csv(file, header = TRUE, check.names = FALSE)
  
  data$R2<-data$beta.exposure*data$beta.exposure*2*(data$eaf.exposure)*(1-data$eaf.exposure)
  data$F<-(data$samplesize.exposure-2)*data$R2/(1-data$R2)
  
  base_name <- tools::file_path_sans_ext(basename(file))
  output_file <- paste0("data_R2_F_", base_name, ".csv")
  write.csv(data[, c("SNP", "R2", "F")], output_file, quote = FALSE, row.names = FALSE)
}
                                                                                
