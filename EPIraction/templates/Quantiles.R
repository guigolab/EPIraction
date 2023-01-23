#!/usr/bin/env Rscript

data_folder = "!{data_folder}"
chromosomes = "!{chromosomes}"
regions     = "!{regions}"

chromosomes = sub('\\[','',chromosomes)
chromosomes = sub('\\]','',chromosomes)

suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(Rfast))
suppressPackageStartupMessages(library(data.table))
options(width=210)

chromosomes = strsplit(chromosomes,", ")[[1]]

data_matrix = fread(file = sprintf("%s/temp/normalize/Consensus.%s.raw.data.gz",data_folder,chromosomes[1]), sep='\t', header=T, stringsAsFactors=FALSE, nThread = 6)

for(index in 2:length(chromosomes))
{
    update = fread(file = sprintf("%s/temp/normalize/Consensus.%s.raw.data.gz",data_folder,chromosomes[index]), sep='\t', header=T, stringsAsFactors=FALSE, nThread = 6)
    merge = list(data_matrix,update)
    data_matrix = rbindlist(merge,use.names = TRUE)
    cat("Merge with",chromosomes[index],"\n")
}

all_regions = fread(file = sprintf("%s/files/%s",data_folder,regions), sep='\t', header=F, stringsAsFactors=FALSE)
colnames(all_regions) = c("chrom","start","end","region")
all_regions[,start:=NULL]
all_regions[,end:=NULL]

data_matrix  = merge(all_regions,data_matrix, by = "region", sort = F)
index_fields = data_matrix[,c("chrom","region")]

data_matrix[,chrom:=NULL]
data_matrix[,region:=NULL]
tissues = colnames(data_matrix)

data_matrix = as.matrix(data_matrix)
colnames(data_matrix) = tissues

cat("\nMaximums Before:\n")
print(sprintf("%.3f",colMaxs(data_matrix, value = T, parallel = T)))

new_matrix = normalize.quantiles(data_matrix, copy=T, keep.names=T)
colnames(new_matrix) = tissues

cat("Done quantiles\n")
proc.time()

cat("\nMaximums after:\n")
print(sprintf("%.3f",colMaxs(new_matrix, value = T, parallel = T)))

complete_data = cbind(index_fields,new_matrix)
setDT(complete_data)
cat("\n")

for(curr_chrom in chromosomes)
{
    chrom_data = complete_data[ chrom == curr_chrom]
    chrom_data[,chrom:=NULL]
    fwrite(chrom_data, file=sprintf("Consensus.%s.norm.data.gz",curr_chrom), sep='\t', col.names=T, quote = F, nThread = 6, compress = 'gzip')
    system(sprintf("mv Consensus.%s.norm.data.gz %s/temp/normalize/Consensus.%s.norm.data.gz", curr_chrom, data_folder, curr_chrom))
    cat("Saved",curr_chrom,"\n")
}
proc.time()
