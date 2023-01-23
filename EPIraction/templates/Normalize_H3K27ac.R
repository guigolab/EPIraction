#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Rfast))

options(width=180)
cores = 3

data_folder = "!{data_folder}"
curr_tissue = "!{tissue}"
chromosomes = "!{chromosomes}"

chromosomes = sub('\\[','',chromosomes)
chromosomes = sub('\\]','',chromosomes)
chromosomes = strsplit(chromosomes,", ")[[1]]
fields = c("region",curr_tissue)

current_chrom  = chromosomes[1]
H3K27ac_matrix = fread(file = sprintf("%s/temp/signal/%s/Raw.%s.data.gz", data_folder, curr_tissue, current_chrom), sep='\t', header=T, stringsAsFactors=FALSE, nThread = cores)
data_scalings  = fread(file = sprintf("%s/temp/normalize/Consensus.%s.scale.data.gz", data_folder, current_chrom), sep='\t', header=T, stringsAsFactors=FALSE, nThread = cores)
data_scalings  = data_scalings[,..fields]
data_scalings[,chromosome:=current_chrom]
colnames(data_scalings) = c("region","scale","chrom")
setcolorder(data_scalings, c("region","chrom","scale"))

data_matrix = merge(data_scalings,H3K27ac_matrix, by = "region", sort = F)
scale = as.numeric(data_matrix$scale)
for(sample_id in 4:ncol(data_matrix))
{
    values = data_matrix[,..sample_id][[1]]
    data_matrix[,sample_id] = values*scale
}
data_matrix[,scale:=NULL]

for(index in 2:length(chromosomes))
{
    current_chrom  = chromosomes[index]
    H3K27ac_matrix = fread(file = sprintf("%s/temp/signal/%s/Raw.%s.data.gz", data_folder, curr_tissue, current_chrom), sep='\t', header=T, stringsAsFactors=FALSE, nThread = cores)
    data_scalings  = fread(file = sprintf("%s/temp/normalize/Consensus.%s.scale.data.gz", data_folder, current_chrom), sep='\t', header=T, stringsAsFactors=FALSE, nThread = cores)
    data_scalings  = data_scalings[,..fields]

    data_scalings[,chromosome:=current_chrom]
    colnames(data_scalings) = c("region","scale","chrom")
    setcolorder(data_scalings, c("region","chrom","scale"))

    work_matrix = merge(data_scalings,H3K27ac_matrix, by = "region", sort = F)
    scale = as.numeric(work_matrix$scale)

    for(sample_id in 4:ncol(work_matrix))
    {
	values = work_matrix[,..sample_id][[1]]
	work_matrix[,sample_id] = values*scale
    }
    work_matrix[,scale:=NULL]
    merge = list(data_matrix,work_matrix)
    data_matrix = rbindlist(merge,use.names = TRUE)
    cat("Merge with",current_chrom,"\n")
}

index_fields = data_matrix[,c("region","chrom")]
data_matrix[,chrom:=NULL]
data_matrix[,region:=NULL]
accessions = colnames(data_matrix)

data_matrix = as.matrix(data_matrix)
colnames(data_matrix) = accessions

cat("\nMaximums Before:\n")
print(sprintf("%.3f",colMaxs(data_matrix, value = T, parallel = T)))

new_matrix = normalize.quantiles(data_matrix, copy=T, keep.names=T)
colnames(new_matrix) = accessions

cat("Done quantiles\n")
cat("\nMaximums after:\n")
print(sprintf("%.3f",colMaxs(new_matrix, value = T, parallel = T)))

complete_data = cbind(index_fields,new_matrix)
setDT(complete_data)
cat("\n")

for(curr_chrom in chromosomes)
{
    chrom_data = complete_data[ chrom == curr_chrom]
    chrom_data[,chrom:=NULL]
    fwrite(chrom_data, file=sprintf("Quantiles.%s.norm.data.gz",curr_chrom), sep='\t', col.names=T, quote = F, nThread = 3, compress = 'gzip')
    system(sprintf("mv -f Quantiles.%s.norm.data.gz %s/temp/signal/%s/Quantiles.%s.data.gz", curr_chrom, data_folder, curr_tissue, curr_chrom))
    cat("Saved",curr_chrom,"\n")
}
proc.time()
