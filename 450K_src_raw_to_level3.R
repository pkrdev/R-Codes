#!/usr/bin/Rscript --slave

########################## Input Data from Command Line ################################################

location <- commandArgs(TRUE)

print(location)

########################## Loading Essential Libraries #################################################

library(stringr)

library(minfi)

library(data.table)

library(futile.logger)

######################### Removing space and newline character from input #############################

location <- str_replace_all(location, "[\r\n]" , "")

print(location)

cr.dir <- getwd()

tempdir <- tempdir()

######################### Fetching list of file from location input ###################################

file_names <-list.files(path = location, full.names = TRUE,recursive = TRUE,pattern="\\.idat$")

fname = basename(file_names)

common_fname <- str_sub(fname, 1, 17)

beta_df <- data.frame()

dir.create("output")


########################## Copying files in temp directory for further processing #####################

for(i in 1:length(file_names)) {

file.copy(file_names[grep(common_fname[i], file_names)], tempdir)

msg <- paste("Getting File.....",fname[i])

print(msg) 

}

setwd(tempdir)

######################### Using minfi library to read idat files #####################################

for(i in 1:length(file_names)) {

setwd(tempdir)

RGSet <- read.metharray(common_fname[i], extended = TRUE, verbose =TRUE, force = TRUE)


######################### Noob Normalization for single sample Normalization #########################

RGSet_noob_nor <- preprocessNoob(rgSet = RGSet,dyeCorr = TRUE,verbose = TRUE,dyeMethod = "single")

grRatio <- ratioConvert(RGSet_noob_nor, what = c("both"), keepCN = TRUE)


######################## Mapping Probes to genome ####################################################

GMset <- mapToGenome(grRatio)

betas <- getBeta(GMset)

rng <- granges(GMset)

rng_df <- data.frame(rng)

beta19 <- cbind(betas,rng_df)

beta19 <- setDT(beta19, keep.rownames = TRUE)[]

colnm <- colnames(beta19)

colnm[1] <- "Composite_Element_REF"

colnm[2] <- "Beta_value"

colnames(beta19) <- colnm

if(!("Composite_Element_REF" %in% names(beta19) & "Beta_value" %in% names(beta19))) {

flog.error("Necessary columns do not exist in input %s", beta19)

}

######################## Reading Annotation file to convert hg19 to hg38 ############################

setwd(cr.dir)

meta = read.table("hm450.manifest.hg38.gencode.v22.tsv",header = TRUE,sep = "\t")

setwd(tempdir)

setnames(meta, c("Chromosome", "Start", "End", "Composite_Element_REF", "Gene_Symbol", "Gene_Type", "Transcript_ID", "Position_to_TSS", "CGI_Coordinate", "Feature_Type"))

######################## Hasing between hg19 and hg38 ##############################################

m = match(beta19$"Composite_Element_REF", meta$Composite_Element_REF)

if(!setequal(m, 1:length(m))) {

flog.warn("Probeset IDs do not match")

}

######################### Writting complete table to a variable ####################################

beta38 = cbind(beta19[, c("Composite_Element_REF", "Beta_value")], meta[m, -4])

setwd(cr.dir)
setwd("output")

####################### Writting methylation beta values to folder ################################## 

write.csv(beta38,paste(common_fname[i],".csv"))

}

setwd(cr.dir)



