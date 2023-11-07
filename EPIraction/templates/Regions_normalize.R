#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(Rfast))
suppressPackageStartupMessages(library(data.table))
options(width=180)
cores = 4

data_folder  = "!{data_folder}"
tissues_list = "!{tissues}"

tissues_list = sub('\\[','',tissues_list)
tissues_list = sub('\\]','',tissues_list)
tissues      = sort(strsplit(tissues_list,", ")[[1]])

#########################################################
### Here I form "activity_matrix" for all my tissues. ###
#########################################################
curent_tissue = tissues[1]
data_matrix   = fread(file = sprintf("%s/temp/tissues/%s.regions.data.gz", data_folder, curent_tissue), sep='\t', header=T, stringsAsFactors=FALSE, nThread = cores, tmpdir="./")

activity_matrix = data_matrix[,c("name","Activity")]
activity_matrix[Activity==0, Activity:= NA]
colnames(activity_matrix)[ncol(activity_matrix)] = curent_tissue

cat("Load ",curent_tissue,"\n")

for(index in 2:length(tissues))
{
    curent_tissue = tissues[index]
    data_matrix   = fread(file = sprintf("%s/temp/tissues/%s.regions.data.gz", data_folder, curent_tissue), sep='\t', header=T, stringsAsFactors=FALSE, nThread = cores, tmpdir="./")
    activity_current  = data_matrix[,c("name","Activity")]
    activity_current[Activity==0, Activity:= NA]

    colnames(activity_current)[ncol(activity_current)] = curent_tissue
    activity_matrix   = merge(activity_matrix,activity_current, by = "name", sort = F)
    cat("Load ",curent_tissue,"\n")
}
cat("Load all tissues\n\n")

##################################################
### Here I normalize.quantiles Activity matrix ###
##################################################
regions = activity_matrix$name
activity_matrix[,name:=NULL]
activity_matrix = as.matrix(activity_matrix)
colnames(activity_matrix) = tissues

new_activity_matrix = normalize.quantiles(activity_matrix, copy=T, keep.names=T)
new_activity_matrix = as.data.frame(cbind(regions,new_activity_matrix))
setDT(new_activity_matrix)
fwrite(new_activity_matrix, file="Activity.data.gz", sep='\t', col.names=T, quote = F, nThread = cores, compress = 'gzip', na="0.0")
rm(new_activity_matrix)

activity_matrix = fread(file = "Activity.data.gz", sep='\t', header=T, stringsAsFactors=FALSE, nThread = cores, tmpdir="./")
activity_matrix[,regions:=NULL]
activity_matrix = as.matrix(activity_matrix)
colnames(activity_matrix) = tissues
activity_max  =  colMaxs(activity_matrix, parallel = T, value = T)
activity_mean = colmeans(activity_matrix, parallel = T)

statistics = data.table(tissue = tissues, mean_activity = sprintf("%.6f",activity_mean), max_activity = sprintf("%.6f",activity_max))
print(statistics)
fwrite(statistics, file = "Activity.log", sep='\t', col.names=T, quote = F)

system(sprintf("mv Activity.data.gz %s/temp/tissues",data_folder))
system(sprintf("cp Activity.log     %s/temp/tissues",data_folder))

warnings()
proc.time()
