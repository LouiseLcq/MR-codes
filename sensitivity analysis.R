library(TwoSampleMR)
library(data.table)
setwd("D:/postgraduate.research/my subject/Mendelian randomization/outcome data Finngen to CIF/CLL/harmonize")
file_list <- list.files()

#多效性检验
for (file in file_list) {
  data<- fread(file)
  pleiotropy <- mr_pleiotropy_test(data)
  write.table(pleiotropy, paste0("D:/postgraduate.research/my subject/Mendelian randomization/outcome data Finngen to CIF/CLL/pleiotropy/pleiotropy_",file,".txt"),sep = "\t", col.names = T, row.names = FALSE)
}

setwd("D:/postgraduate.research/my subject/Mendelian randomization/outcome data Finngen to CIF/CLL/pleiotropy")
file_list <- list.files()
for (file in file_list) {
  pleiotropy<- fread(file)
pleiotropy$inter1 <- round(pleiotropy$egger_intercept - 1.96*pleiotropy$se,3)
pleiotropy$inter2 <- round(pleiotropy$egger_intercept + 1.96*pleiotropy$se,3)
pleiotropy$egger_intercept_CI<-paste0(pleiotropy$inter1,',',pleiotropy$inter2) 
write.table(pleiotropy, paste0("D:/postgraduate.research/my subject/Mendelian randomization/outcome data Finngen to CIF/CLL/egger/egger_",file),sep = "\t", col.names = T, row.names = FALSE)}
#这个pleiotropy需要p＞0。05，且egger足够小，最好置信区间跨越0

#异质性检验
setwd("D:/postgraduate.research/my subject/Mendelian randomization/outcome data Finngen to CIF/CLL/harmonize")
file_list <- list.files()
for (file in file_list) {
  data<- fread(file)
het <- mr_heterogeneity(data)
write.table(het, paste0("D:/postgraduate.research/my subject/Mendelian randomization/outcome data Finngen to CIF/CLL/heterogeneity/heterogeneity_",file),sep = "\t", col.names = T, row.names = FALSE)}
#这个het不需要算置信区间，p＞0.05就行

#result_once <- mr(dat,method_list=subset(mr_method_list(),name=='Inverse variance weighted'|name=='MR Egger')$obj)
setwd("D:/postgraduate.research/my subject/Mendelian randomization/outcome data Finngen to CIF/CLL/mr")
file_list <- list.files()
for (file in file_list) {
  data<- fread(file)
  data$b_l<-data$b-1.96*data$se
  data$b_h<-data$b+1.96*data$se
  data$OR<-exp(data$b)
  data$OR_l<-exp(data$b_l)
  data$OR_h<-exp(data$b_h)
  write.table(data, paste0("D:/postgraduate.research/my subject/Mendelian randomization/outcome data Finngen to CIF/CLL/OR/OR_",file),sep = "\t", col.names = T, row.names = FALSE)}

