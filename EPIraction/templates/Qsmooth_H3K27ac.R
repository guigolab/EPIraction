#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(qsmooth))
suppressPackageStartupMessages(library(Rfast))
options(width=180)
cores = 3

data_folder  = "!{data_folder}"
tissues_list = "!{tissues}"
chromosome   = "!{chrom}"
tissues_list = sub('\\[','',tissues_list)
tissues_list = sub('\\]','',tissues_list)
tissues      = strsplit(tissues_list,", ")[[1]]

scalings = fread(file = sprintf("%s/temp/normalize/Consensus.%s.scale.data.gz", data_folder, chromosome), sep='\t', header=T, stringsAsFactors=FALSE, nThread = cores)

current_tissue = tissues[1]
data_matrix = fread(file = sprintf("%s/temp/signal/%s.%s.data.gz", data_folder, tissues[1], chromosome), sep='\t', header=T, stringsAsFactors=FALSE, nThread = cores)
setDT(data_matrix)

fields = c("region",current_tissue)

current_scalings = scalings[,..fields]
colnames(current_scalings) = c("region","scale")

data_matrix = merge(current_scalings,data_matrix, by = "region", sort = F)
scale = as.numeric(data_matrix$scale)
for(sample_id in 3:ncol(data_matrix))
{
    values = data_matrix[,..sample_id][[1]]
    data_matrix[,sample_id] = values*scale
}
data_matrix[,scale:=NULL]
group_factor = rep(current_tissue,ncol(data_matrix)-1)
cat("Added",current_tissue,"\n")

for(index in 2:length(tissues))
{
    current_tissue = tissues[index]
    work_matrix = fread(file = sprintf("%s/temp/signal/%s.%s.data.gz", data_folder, current_tissue, chromosome), sep='\t', header=T, stringsAsFactors=FALSE, nThread = cores)
    setDT(work_matrix)
    fields = c("region",current_tissue)
    current_scalings = scalings[,..fields]
    colnames(current_scalings) = c("region","scale")

    work_matrix = merge(current_scalings,work_matrix, by = "region", sort = F)
    scale = as.numeric(work_matrix$scale)
    for(sample_id in 3:ncol(work_matrix))
    {
	values = work_matrix[,..sample_id][[1]]
	work_matrix[,sample_id] = values*scale
    }
    work_matrix[,scale:=NULL]

    data_matrix = merge(data_matrix,work_matrix, by = "region", sort = F)
    group_factor = c(group_factor,rep(current_tissue,ncol(work_matrix)-1))
    cat("Added",tissues[index],"\n")
}
regions = data_matrix$region
data_matrix[,region:=NULL]
data_matrix = as.matrix(data_matrix)
statistics = cbind(accession = colnames(data_matrix),tissue = group_factor, before_Max = sprintf("%.6f",colMaxs(data_matrix, value = T, parallel = T)), before_Mean = sprintf("%.6f",colmeans(data_matrix, parallel = T)))

dat_qs = qsmooth(object = data_matrix, group_factor = group_factor, window = 0.01)
normalized = qsmoothData(dat_qs)
colnames(normalized) = colnames(data_matrix)
report_matrix = cbind(region = regions,normalized)
report_matrix = as.data.table(report_matrix)
fwrite(report_matrix,file="Current.matrix.gz", sep='\t', col.names=T, quote = F, nThread = cores, compress = "gzip")

statistics = cbind(statistics, after_Max = sprintf("%.6f",colMaxs(normalized, value = T, parallel = T)), after_Mean = sprintf("%.6f",colmeans(normalized, parallel = T)))
statistics = as.data.table(statistics)
fwrite(statistics,file="Statistics.data", sep='\t', col.names=T, quote = F)
system(sprintf("mv Current.matrix.gz %s/data/yyy.H3K27ac.%s.quantiles.gz",data_folder,chromosome))
system(sprintf("mv Statistics.data %s/data/yyy.H3K27ac.%s.statistics",data_folder,chromosome))

proc.time()
