####p adjust#######
setwd('D:/postgraduate.research/my subject/Mendelian randomization/outcome data Finngen to CIF/AML/4mr')
file_list <- list.files()
tab<-NULL
for (file in file_list) {
  data <- read.csv(file, header = TRUE, check.names = FALSE)
  x<-data.frame(data)
  tab<-rbind(tab,x)
}
tab1<-tab[tab$method=="Inverse variance weighted"|tab$method=="Wald ratio",]

is.numeric(tab1$pval)
tab1$pval_adjust<-p.adjust(tab1$pval,'fdr')
file_path <- file.path("D:/postgraduate.research/my subject/Mendelian randomization/outcome data Finngen to CIF/mrsummary91", "tab1_91_AML.csv")
write.csv(tab1, file_path)