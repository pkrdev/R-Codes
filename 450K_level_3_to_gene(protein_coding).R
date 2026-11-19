#!/usr/bin/Rscript --slave

library(stringr)

library(data.table)

library(futile.logger)

library(MASS)


cr.dir <-list.files(path=".", full.names = TRUE,recursive = TRUE,pattern="\\gdc_hg38.txt$")

dir.create("output")

G_list <- read.csv("matched_ens_genes.csv",header=TRUE,stringsAsFactors=FALSE)

G_list <- G_list[-1]

for(i in 1:length(cr.dir)) {
 file.copy(cr.dir,"output")
 }

meta_data <- read.csv("meta_gbm_methy_pri.csv",stringsAsFactor=FALSE,header=TRUE)

setwd(output)


for(i in 1:length(cr.dir)) {

print("##################################   Reading File using parallel processing ########################")

require(data.table)

file_name <- meta_data$File.Name[i]

system.time(data <- fread(file_name))

data <- as.data.frame(data)

print("##################################   Reading File Successful ########################")

print("##################################   Extracting Protein Coding CpG probes   ########################")

protein_coding_index <-which(data$Gene_Type %like% "protein_coding")

protein_coding <- data[protein_coding_index,]

loop_count <- length(protein_coding_index) 

##################################   Removing Raw Data to Release Memory   #########################

print("##################################   Removing Raw Data to Release Memory   #########################")

remove(data)

##################################   Pre defined Function for NA handling   ########################

check <- function(value) {

if (length(value)==0) {

return(as.character(NA))

}

else {

return(value)

}

}

################################# Replacing GeneName with ensemble id in Protein Coding File ########

print("################################# Replacing GeneName with ensemble id in Protein Coding File ########")

start_time <- Sys.time()

symb <- NULL;

for(i in 1:loop_count) {

x <-strsplit(protein_coding$Gene_Symbol[i],split=";",fixed=TRUE)

y <- unlist(x, use.names=FALSE)

y <-unique(y)

l <- length(y)

for(i in 1:l) {

m <- G_list[G_list$hgnc_symbol %in% y[i],]$ensembl_gene_id

y[i] <- check(m)

}

z <-paste(y, collapse = ";")

symb <- c(symb,z)

}

protein_coding$Gene_Symbol <- symb

protein_coding <- as.data.frame(protein_coding)

end_time <- Sys.time()

print((start_time - end_time))

############################### Fetching Enesmble ids #################################################

############################### Using parallel processing #############################################

print("############################### Fetching Enesmble ids Using Parallel Processing ###################")

Ensemble_Id <- NULL;

library(doParallel)

cr <- detectCores()

registerDoParallel(cores=(cr-2))

Ensemble_Id <- foreach(i=1:loop_count, .combine = 'c') %dopar% {

unique(unlist(strsplit(protein_coding$Gene_Symbol[i],split=";",fixed=TRUE), use.names=FALSE))

}

Ensemble_Id <- unique(Ensemble_Id)

Ensemble_Id <- Ensemble_Id[-which(Ensemble_Id=="NA")]

################## Parallel Processing Avg Beta Values Fetching #######################################

print("################## Parallel Processing Avg Beta Values Fetching #######################################")

betas <- NULL;

library(doParallel)

cr <- detectCores()

registerDoParallel(cores=(cr-2))

betas <- foreach(i=1:length(Ensemble_Id), .export=c('like'), .packages='data.table') %dopar% {

mean(protein_coding[which(protein_coding$Gene_Symbol %like% Ensemble_Id[i]),]$Beta_value)

}

betas <- unlist(betas, use.names=FALSE)

################# Creating Data Frame of Gene Vs Beta Values #########################################

print("################# Creating Data Frame of Gene Vs Beta Values #########################################")

gene_wise <- data.frame(Ensemble_Id,betas)

out_file_name <-meta_data[which(meta_data$File.Name==file_name),]$Case.ID

out_file_name <- paste0(out_file_name,".csv")

out_file_name <-str_replace_all(out_file_name, " " ,"")


write.csv(gene_wise,out_file_name)

}





