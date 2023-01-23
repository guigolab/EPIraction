#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
options(width=180)
cores = 4

data_folder  = "!{data_folder}"
tissues_list = "!{tissues}"
print(tissues_list)
tissues_list = sub('\\[','',tissues_list)
tissues_list = sub('\\]','',tissues_list)
tissues      = sort(strsplit(tissues_list,", ")[[1]])
chromosome   = "!{chrom}"

cat("Data: ",data_folder,"\n")
cat("Chrom:",chromosome,"\n")

data_matrix = fread(file = sprintf("%s/temp/signal/%s/Quantiles.%s.data.gz", data_folder, tissues[1], chromosome), sep='\t', header=T, stringsAsFactors=FALSE, nThread = cores)
setDT(data_matrix)
cat("Added",tissues[1],"\n")

for(index in 2:length(tissues))
{
    work_matrix = fread(file = sprintf("%s/temp/signal/%s/Quantiles.%s.data.gz", data_folder, tissues[index], chromosome), sep='\t', header=T, stringsAsFactors=FALSE, nThread = cores)
    setDT(work_matrix)
    data_matrix = merge(data_matrix,work_matrix, by = "region", sort = F)
    cat("Added",tissues[index],"\n")
}
fwrite(data_matrix,file=sprintf("%s/data/yyy.H3K27ac.%s.matrix.gz",data_folder,chromosome), sep='\t', col.names=T, quote = F, nThread = cores, compress = "gzip")

proc.time()
