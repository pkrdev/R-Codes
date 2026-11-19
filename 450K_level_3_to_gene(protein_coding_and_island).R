#!/usr/bin/Rscript --slave

library(stringr)

library(data.table)

library(futile.logger)

library(MASS)

files <-list.files(path=".", full.names = TRUE,recursive = TRUE,pattern="\\gdc_hg38.txt$")

dir.create("output_new")

G_list <- read.delim("gene_data.csv",sep="\t",header=TRUE,stringsAsFactor=FALSE)


for(i in 1:length(files)) {
 
file.copy(files,"output_new")

}

meta_data <- read.csv("meta_gbm_methy_pri.csv",stringsAsFactor=FALSE,header=TRUE)

setwd("output_new")


for(i in 1:length(files)) {

print("##################################   Reading File using parallel processing ########################")

require(data.table)

file_name <- meta_data$File.Name[i]

system.time(data <- fread(file_name))

data <- as.data.frame(data)

data <- na.omit(data)

print("##################################   Reading File Successful ########################")

print("##################################   Extracting Protein Coding CpG probes   ########################")

protein_coding_index <-which(data$Gene_Type %like% "protein_coding")

protein_coding <- data[protein_coding_index,]

protein_coding_island_index <-which(protein_coding$Feature_Type %like% "Island")

protein_coding_island <- protein_coding[protein_coding_island_index,]

loop_count <- length(protein_coding_island_index) 

protein_coding <-protein_coding_island


##################################   Removing Raw Data to Release Memory   #########################

print("##################################   Removing Raw Data to Release Memory   #########################")

remove(data)

remove(protein_coding_island)

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

m <- G_list[G_list$Gene_Name %in% y[i],]$ensembl_gene_id

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

protein_ens <-G_list$ensembl_gene_id

unmatch_ids <- protein_ens[which(!protein_ens %in% Ensemble_Id)]

unmatch_count <- length(protein_ens[(which(!protein_ens %in% Ensemble_Id))])

NA_container <- rep("NA",unmatch_count)

df <- data.frame(as.character(unmatch_ids),NA_container)

colnames(df) <-colnames(gene_wise)


df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE)

final <-rbind(gene_wise,df)

write.csv(final,"temp.csv")

final<- read.csv("temp.csv",stringsAsFactor=FALSE,header=TRUE)

final <- final[-1]

final<-final[order(final$Ensemble_Id),]

out_file_name <-meta_data[which(meta_data$File.Name==file_name),]$Case.ID

out_file_name <- paste0(out_file_name,".csv")

out_file_name <-str_replace_all(out_file_name, " " ,"")

write.csv(final,out_file_name)

}





