##### MR Merge according rsids- example #####
setwd("D:/postgraduate.research/my subject/Mendelian randomization/exposed data/output cytokines/F 10 r^2")
list.files()
file_list <- list.files()
BNGF<- read.table("F 10 r^2LD_B_NGF_F.txt",header = T)
View(BNGF)
library(data.table)
setwd("D:/postgraduate.research/my subject/Mendelian randomization/outcome data/finngen_R10_C3_CML_EXALLC")
out<-fread('finngen_R10_C3_CML_EXALLC')
out$SNP <- out$rsids
out$'chr:pos'<-paste(out$`#chrom`, out$pos, sep = ":")
out$phenotype<-'CML'
View(out)
BNGFout<-merge(BNGF,out,by = "SNP")
View(BNGFout)


##### MR Merge according rsids- cycle #####
for (file in file_list) {
  file=file_list[1]
  BNGF<- read.table(file,header = T)
  BNGFout<-merge(BNGF,out,by = "SNP")
  write.csv(BNGFout,paste0('D:/postgraduate.research/my subject/Mendelian randomization/CML/MR outcome/Merge/Merge_',file,'.csv'))
}

sum(out$SNP== "")
outblank <- out[out$SNP == "", ]

###### Merge outblank and snp150database defeat#####
######setwd("D:/postgraduate.research/my subject/Mendelian randomization/MR process learning resource/snp150_hg19.txt")
######snp150 <- fread("snp150_hg19.txt")
######outblank$`chromosome:start` <- paste(outblank$`#chrom`,outblank$pos,sep = ":")
#fillblank <- merge(outblank, snp150, by ="pos")
#out<-merge(fillblank,out,by="pos)
# rerun the MR Merge according rsids- cycle

###### use biomaRt to see snpsblank######
library(biomaRt)
outblanksmall <- head(outblank,10)
snpmart <- useMart(biomart = "ENSEMBL_MART_SNP",dataset = "hsapiens_snp")
snps<-getBM(attributes = c('refsnp_id','chr_name','chrom_start'),
            filters = c('chr_name','start'),
            values = list(outblanksmall$`#chrom`,outblanksmall$pos),
            mart = snpmart
            
            )
# see the first few lines in the snps（fillblank) data

##### use biomaRt to see location-example######
setwd("D:/postgraduate.research/my subject/Mendelian randomization/exposed data/output cytokines/F 10 r^2")
BNGF<- read.table("F 10 r^2LD_B_NGF_F.txt",header = T)
View(BNGF)
ids <- BNGF[, 1]
snp_attributes <- c("refsnp_id","chr_name", "chrom_start")
snp_locations <- getBM(attributes = snp_attributes, filters = "snp_filter", 
                       values = ids, mart = snpmart)
merge_BNGF<-merge(BNGF,snp_locations,by.x="SNP",by.y = "refsnp_id")
write.table(merge_BNGF,col.names = T,row.names = F,sep = "\t")

###### use biomaRt to see location-cycle NO REPEAT#####
setwd("D:/postgraduate.research/my subject/Mendelian randomization/exposed data/output cytokines/F 10 r^2")
list.files()
file_list <- list.files()
snp_attributes <- c("refsnp_id","chr_name", "chrom_start")
for (file in file_list){
  data<-read.table(file,header = T)
  ids <- data[, 1]
  snp_locations <- getBM(attributes = snp_attributes, filters = "snp_filter", 
                         values = ids, mart = snpmart)
  merge_data<-merge(data,snp_locations,by.x="SNP",by.y = "refsnp_id")
  write.table(merge_data,paste0("D:/postgraduate.research/my subject/Mendelian randomization/exposed data/output cytokines/exp_chr/chr_",file),col.names = T,row.names = F,sep = "\t")
}

##### merge the exposure(loc) and CML outcome with locations-example#####
setwd("D:/postgraduate.research/my subject/Mendelian randomization/exposed data/output cytokines/exp_chr")
BNGF<-read.table("chr_F 10 r^2LD_B_NGF_F.txt",header = T)
View(BNGF)
BNGF$'chr:pos' <- paste(BNGF$chr_name, BNGF$chrom_start, sep = ":")
out$'chr:pos'<-paste(out$`#chrom`, out$pos, sep = ":")
merge_BNGF<-merge(BNGF, out, by = "chr:pos")

##### merge the exposure(loc) and CML outcome with locations-cycle START there#####
setwd("D:/postgraduate.research/my subject/Mendelian randomization/exposed data/output cytokines/exp_chr")
list.files()
file_list <- list.files()
out$'chr:pos'<-paste(out$`#chrom`, out$pos, sep = ":")
for (file in file_list){
  data<-read.table(file,header = T) 
  data$'chr:pos' <- paste(data$chr_name, data$chrom_start, sep = ":")
  merge_data<-merge(data, out, by = "chr:pos")
  write.table(merge_data,paste0("D:/postgraduate.research/my subject/Mendelian randomization/CML/MR outcome/Merge/Merge_",file),col.names = T,row.names = F,sep = "\t")
}

#####format merged outcome#####
setwd("D:/postgraduate.research/my subject/Mendelian randomization/outcome data Finngen to cytokines/CML/MR outcome/Merge")
list.files()
file_list <- list.files()
for (file in file_list) {
  out_data <- fread(file)
  outcome <- format_data(
    dat= out_data,
    type = "outcome",
    snps = out_data$SNP,
    header = TRUE,
    phenotype_col = "phenotype",
    snp_col = "SNP.y",
    beta_col = "beta",
    se_col = "sebeta",
    eaf_col="af_alt",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    pval_col = "pval"
  )
  output_file <- paste0("D:/postgraduate.research/my subject/Mendelian randomization/CML/MR outcome/format/format_", file, ".txt")
  
  # 将结果写入输出文件
  #write.table(outcome, output_file, sep = "\t", col.names = T, row.names = FALSE)
  
}
###Warning message:In .fun(piece, ...) :Duplicated SNPs present in exposure data for phenotype 'CML. Just keeping the first instance:rs6764884###

#####harmonise-cycle#####
library(TwoSampleMR)
setwd("D:/postgraduate.research/my subject/Mendelian randomization/exposed data/output cytokines/F 10 r^2")
list.files()
file_list <- list.files()
for (file in file_list) {
  exp_data <- read.table(file, header = TRUE, check.names = FALSE, sep = "\t")
  changedname <- sub("^F 10 r\\^2LD_(.*?)_F\\.txt$", "\\1", file)
  out_file <- paste0("D:/postgraduate.research/my subject/Mendelian randomization/CML/MR outcome/format/format_Merge_chr_F 10 r^2LD_", changedname, "_F.txt.txt")
  out_data <- read.table(out_file, header = TRUE, sep = "\t")
  harmonised_data <- harmonise_data(exp_data, out_data, action = 2)
  output_file <- paste0("D:/postgraduate.research/my subject/Mendelian randomization/CML/MR outcome/harmonize/harmonize_", changedname, ".txt")
  write.table(harmonised_data, file = output_file, sep = "\t", col.names = T, row.names = FALSE)
}

####mr#######
setwd("D:/postgraduate.research/my subject/Mendelian randomization/CML/MR outcome/harmonize")
file_list <- list.files()
for (file in file_list) {
  data <- read.table(file, header = TRUE, check.names = FALSE, sep = "\t")
  changedname <- sub("^harmonised_(.*)$", "\\1", file)
  mrdata <- mr(data)
  output_file <- paste0("D:/postgraduate.research/my subject/Mendelian randomization/CML/MR outcome/mr/mr_", changedname, ".txt")
  write.table(mrdata, file = output_file, sep = "\t", col.names = T, row.names = FALSE)
}

####p adjust#######
library(stringr)
setwd('D:/postgraduate.research/my subject/Mendelian randomization/outcome data Finngen to CIF/CLL/4mr')
file_list <- list.files()
tab<-NULL
for (file in file_list) {
  data <- read.table(file, header = TRUE, check.names = FALSE, sep = "\t")
  x<-data.frame(data)
  tab<-rbind(tab,x)
}
tab1<-tab[tab$method=="Inverse variance weighted"|tab$method=="Wald ratio",]

is.numeric(tab1$pval)
tab1$pval_adjust<-p.adjust(tab1$pval,'fdr')
#tab$probeID<-str_extract(tab$V1,"OID\\d+")

# for (i in 1:length(file_list)) {
#     data <- read.table(file_list[i], header = TRUE, check.names = FALSE, sep = "\t")
#     print(i)
#     x<-data.frame(data)
#     tab<-rbind(tab,x)
# }
# 
# tab1<-subset(tab,method=='Inverse variance weighted')

#####sensitivity analysis######


