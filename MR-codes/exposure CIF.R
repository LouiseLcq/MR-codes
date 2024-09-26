######pick exposure data#####
#####install packaged#####
library(TwoSampleMR)
library(biomaRt)
library(data.table)

#####input files######
setwd("D:/postgraduate.research/my subject/Mendelian randomization/exposed data CIF/GWAS CIF")
file_list<-list.files()


#####pick suitable p value######
pFilter = 5e-8
for (file in file_list){
  exp_data <- read.table(file, header = TRUE, check.names = FALSE, sep = "\t")
  exp_data_p <- subset(exp_data, exp_data$`p_value` < pFilter)
  changedname <- sub("^.*([0-9]{3})\\.csv$", "\\1", file)
  output_file <- paste0("D:/postgraduate.research/my subject/Mendelian randomization/exposed data CIF/output CIF/p value/5e-8_", changedname)
  write.table(exp_data_p, file = output_file, sep = "\t", col.names = T, row.names = FALSE)
}


######compensate NA rsids(didn't finish)#######
setwd("D:/postgraduate.research/my subject/Mendelian randomization/exposed data CIF/output CIF/p value")
file_list<-list.files()
for (file in file_list){
  exp_data <- read.table(file, header = TRUE, check.names = FALSE, sep = "\t")
  na_rows <- is.na(file$rsid)
  for (i in which(na_rows)) {
    if (is.na(file$rsid[i])) {
      chromosome <- file$chromosome[i]
      base_pair_location <- file$base_pair_location[i]
      mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                      dataset = "hsapiens_gene_ensembl")
      filter_query <- paste0(chromosome, ":", base_pair_location, "-", base_pair_location)
      rsid_result <- getBM(attributes = c("snp"),
                           filters = c("region"),
                           values = filter_query,
                           mart = mart)
      if (nrow(rsid_result) > 0) {
        fileA$rsid[i] <- rsid_result$snp[1]
      }
    }
  }
}
######de LD######
setwd("D:/postgraduate.research/my subject/Mendelian randomization/exposed data CIF/output CIF/p value")
file_list<-list.files()
for (file in file_list){
  exp_data <- read.table(file, header = TRUE, check.names = FALSE, sep = "\t")
  changedname <- sub("^.*([0-9]{3})\\.csv$", "\\1", file)
  exp_data$PHENO <- changedname
  formatted_data <- format_data(exp_data,
                                type = "exposure",
                                phenotype_col = "PHENO",
                                snp_col = "rsid",
                                beta_col = "beta",
                                se_col = "standard_error",
                                pval_col = "p_value",
                                eaf_col = "effect_allele_frequency",
                                samplesize_col = "n",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele")
  clumped_data <- clump_data(
    formatted_data,
    clump_kb = 10000,
    clump_r2 = 0.1,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EUR"
  )
  output_file <- paste0("D:/postgraduate.research/my subject/Mendelian randomization/exposed data CIF/output CIF/Linkage disequilibrium/LD_", changedname)
  write.table(clumped_data, output_file, sep = "\t", col.names = T, row.names = FALSE)
}

###### add code to deal with clump 502 ########
setwd("D:/postgraduate.research/my subject/Mendelian randomization/exposed data CIF/output CIF/p value")
file_list<-list.files()
max_retry_per_file <- 10 
for (file in file_list){
  success <- FALSE
  retry_count <- 0
  last_error_message <- NULL
  while(!success && retry_count < max_retry_per_file) {
    tryCatch({
      exp_data <- read.table(file, header = TRUE, check.names = FALSE, sep = "\t")
      changedname <- sub("^.*([0-9]{3})\\.csv$", "\\1", file)
      exp_data$PHENO <- changedname
      formatted_data <- format_data(exp_data,
                                    type = "exposure",
                                    phenotype_col = "PHENO",
                                    snp_col = "rsid",
                                    beta_col = "beta",
                                    se_col = "standard_error",
                                    pval_col = "p_value",
                                    eaf_col = "effect_allele_frequency",
                                    samplesize_col = "n",
                                    effect_allele_col = "effect_allele",
                                    other_allele_col = "other_allele")
      clumped_data <- clump_data(
        formatted_data,
        clump_kb = 10000,
        clump_r2 = 0.1,
        clump_p1 = 1,
        clump_p2 = 1,
        pop = "EUR"
      )
      output_file <- paste0("D:/postgraduate.research/my subject/Mendelian randomization/exposed data CIF/output CIF/Linkage disequilibrium/LD_", changedname)
      write.table(clumped_data, output_file, sep = "\t", col.names = T, row.names = FALSE)
      cat("Processed file:", file, "\n")
      success <- TRUE 
    }, error = function(e) {
      last_error_message <- conditionMessage(e)
      cat("Error processing file:", file, "\n")
      cat("Retrying...\n")
    })
    retry_count <- retry_count + 1
  }
  
  if (!success && !is.null(last_error_message)) {
    cat("Final error message for file:", file, "\n")
    cat(last_error_message, "\n")
  }
}
  if (!success && !is.null(last_error_message)) {
    cat("Final error message for file:", file, "\n")
    cat(last_error_message, "\n")
  }
}



######calculate total lines########
folder_path <- "D:/postgraduate.research/my subject/Mendelian randomization/exposed data CIF/output CIF/Linkage disequilibrium"
total_lines <- 0
file_names <- list.files(folder_path, pattern = "\\.tsv$")
for (filename in file_names) {
  file_path <- file.path(folder_path, filename)
  lines <- length(readLines(file_path)) - 1 
  total_lines <- total_lines + lines
}
print(paste("total line is:", total_lines))

