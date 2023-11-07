#!/usr/bin/env Rscript

data_folder = "!{data_folder}"
tissue      = "!{tissue}"

cores = 2
suppressPackageStartupMessages(library(data.table))

options(width=210)
options(scipen=999)

################################
### I load all matrices here ###
################################
complete_regions = fread(file = sprintf("%s/temp/tissues/%s.regions.data.gz",data_folder,tissue), sep='\t', header=T, stringsAsFactors=FALSE, blank.lines.skip=T, nThread = cores, tmpdir="./")
activity_data    = fread(sprintf("%s/temp/tissues/Activity.data.gz",data_folder),                 sep='\t', header=T, stringsAsFactors=FALSE, blank.lines.skip=T, nThread = cores, tmpdir="./")

cat(sprintf("Load data for %s\n",tissue))

##############################################
### I replace Activity with normalized one ###
##############################################
fields = c("regions",tissue)

activity_tissue = activity_data[,..fields]
colnames(activity_tissue) = c("name","Activity.norm")
rm(activity_data)

complete_regions = merge(complete_regions, activity_tissue, by = "name", sort = F)

setcolorder(complete_regions, c("#chrom","start","end","name","score","strand","type","H3K27ac","Open","Cofactor","TPMs","Symbol","Sample","Activity.norm","Activity"))
complete_regions[,Activity:=NULL]
colnames(complete_regions)[14] = "Activity"

cat(sprintf("Wrote activity for %s\n",tissue))

###############################
### Round values and report ###
###############################
complete_regions$H3K27ac  = sprintf("%.4f",complete_regions$H3K27ac)
complete_regions$Open     = sprintf("%.4f",complete_regions$Open)
complete_regions$Cofactor = sprintf("%.4f",complete_regions$Cofactor)
complete_regions$TPMs     = sprintf("%.4f",complete_regions$TPMs)
complete_regions$Activity = sprintf("%.7f",complete_regions$Activity)
complete_regions[Activity == "0.0000000", Activity:="0.0"]

fwrite(complete_regions, file = "Tissue.data", sep = '\t', quote=F, na="0.0")
system(sprintf("mv Tissue.data %s/data/%s.regions.data", data_folder, tissue))

proc.time()
warnings()
