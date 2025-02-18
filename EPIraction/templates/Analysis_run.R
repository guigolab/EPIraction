#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Rfast))
setDTthreads(threads = 1)

options(width=190)
options(scipen=999)
set.seed(123)

data_folder = "!{data_folder}"
cur_tissue  = "!{tissue}"
task_index  = "!{task_index}"
chromosome  = "!{chromosome}"
genes_list  = "!{list}"

file_out = "Pairs.data"
write(paste(c("human_ensembl","enhancer","tissue_abc","closest_abc","complete_abc"), collapse='\t'), file=file_out)

log_data = "Log.analysis"
write(sprintf("Analysis report for %s\n",task_index), file = log_data)

for(gene_id in strsplit(genes_list,";")[[1]])
{
##################################################
### If we have "$gene_id.input.matrix.gz" file ###
##################################################
    if(file.exists(sprintf("%s/Genes/%s/%s.input.matrix.gz", data_folder, chromosome, gene_id)) == FALSE)
    {
        cat(sprintf("absent %s.enhancers.input.gz\n\n", gene_id))
        next
    }

    cat(sprintf("Do %s gene in %s\n", gene_id, cur_tissue))
    write(sprintf("Do %s gene in %s", gene_id, cur_tissue),file=log_data, append=T)

################################
### I load the training data ###
################################
    feature_data = fread(file = sprintf("%s/Genes/%s/%s.input.matrix.gz",data_folder,chromosome,gene_id), sep='\t', header=T, stringsAsFactors=FALSE, blank.lines.skip=T, tmpdir="./")
    cat("load features\n")

################################################################
### I split training data into "tissue_data" and "rest_data" ###
################################################################
    total_features = ncol(feature_data) - 2
    total_samples  = nrow(feature_data)

    tissue_data    = feature_data[ tissue == cur_tissue]
    tissue_tissues = tissue_data$tissue
    tissue_exp     = tissue_data$expression
    tissue_data[,tissue:=NULL]
    tissue_data[,expression:=NULL]
    tissue_matrix  = as.matrix(tissue_data)
    tissue_samples = nrow(tissue_matrix)
    all_features   = colnames(tissue_matrix)

    rest_data      = feature_data[ tissue != cur_tissue]
    rest_tissues   = rest_data$tissue
    rest_exp       = rest_data$expression
    rest_data[,tissue:=NULL]
    rest_data[,expression:=NULL]
    rest_matrix    = as.matrix(rest_data)
    rest_samples   = nrow(rest_matrix)

    feature_data[,tissue:=NULL]
    feature_data[,expression:=NULL]
    complete_matrix = as.matrix(feature_data)
    complete_abc    = colmeans(complete_matrix, parallel = FALSE)

########################################################
### I select enhancers that are active in our tissue ###
########################################################
    tissue_abc    = colmeans(tissue_matrix, parallel = FALSE)
    susbtantial   = tissue_abc > 0
    work_features = sum(susbtantial)
    if(work_features < 10)
    {
	cat("We need at least 10 features with nonzero ABC scores\n")
	write("We need at least 10 features with nonzero ABC scores", file=log_data, append=T)

	report_data = data.table(gene_id = gene_id, enhancer = all_features, score = tissue_abc, closest = sprintf("%.10f",0), rest = sprintf("%.10f",complete_abc))
	report_data = report_data[score > 0]
	setorder(report_data,-score)
	report_data$score   = sprintf("%.10f",report_data$score)

	all_features = report_data$enhancer
	max_char = max(nchar(all_features))
	for(index in 1:length(all_features))
	{
	    while(nchar(all_features[index]) < max_char)
	    {
		all_features[index] = paste(c(" ", all_features[index]), sep='', collapse='')
	    }
	}
	report_data$enhancer = all_features
	
	fwrite(report_data, file=file_out, sep='\t', quote=F, append=T, col.names=F)
	write("####",file=file_out, append=T)
	write("conventional_ABC\n####",file=log_data, append=T)
	cat("conventional_ABC\n####\n")
	next
    }
    cat(sprintf("Consider %d features of %d total\n",work_features,total_features))
    write(sprintf("Consider %d features of %d total",work_features,total_features), file=log_data, append=T)

    tissue_matrix = tissue_matrix[,susbtantial]
    all_features  = colnames(tissue_matrix)
    tissue_abc    = tissue_abc[susbtantial]
    rest_matrix   = rest_matrix[,susbtantial]
    complete_abc  = complete_abc[susbtantial]

###############################################################################################
### I calculate the pearson correlation between current tissue "tissue_abc" and all samples ###
###############################################################################################
    sample_correlation = vector(,nrow(rest_matrix))
    for(index in 1:nrow(rest_matrix))
    {
	opponent_values = as.numeric(rest_matrix[index,])
	if(sd(opponent_values) == 0)
	{
	    sample_correlation[index] = 0
	    next
	}
	analysis = cor.test(opponent_values, tissue_abc)
	sample_correlation[index] = analysis$estimate
    }
    sample_data = data.table(tissue = rest_tissues, expression = rest_exp, correlation = sample_correlation)
    aggregation = sample_data[, .(expression = mean(expression), correlation = mean(correlation)), by=tissue]
    setorder(aggregation,-correlation)
    aggregation = aggregation[correlation > 0.85 & expression > 1 ]
    if(nrow(aggregation) > 3)
	aggregation = aggregation[1:3,]

#################################################################################
### If there is no tissues with correlated enhancers and gene being expressed ###
#################################################################################
    if(nrow(aggregation) == 0)
    {
	cat("There is no tissues with correlated enhancers and gene being expressed\n")
	write("There is no tissues with correlated enhancers and gene being expressed", file=log_data, append=T)

	report_data = data.table(gene_id = gene_id, enhancer = all_features, score = tissue_abc, closest = sprintf("%.10f",0), rest = sprintf("%.10f",complete_abc))
	report_data = report_data[score > 0]
	setorder(report_data,-score)
	report_data$score   = sprintf("%.10f",report_data$score)

	all_features = report_data$enhancer
	max_char = max(nchar(all_features))
	for(index in 1:length(all_features))
	{
	    while(nchar(all_features[index]) < max_char)
	    {
		all_features[index] = paste(c(" ", all_features[index]), sep='', collapse='')
	    }
	}
	report_data$enhancer = all_features

	fwrite(report_data, file=file_out, sep='\t', quote=F, append=T, col.names=F)
	write("####",file=file_out, append=T)
	write("conventional_ABC\n####",file=log_data, append=T)
	cat("conventional_ABC\n####\n")
	next
    }
    cat("correlated tissues:\n")
    write("correlated tissues:", file=log_data, append=T)
    print(aggregation)
    fwrite(aggregation, file=log_data, append=T, sep='\t', row.names=T, quote=F)

##########################################################################
### I calculate average ABC score across correlated tissues and report ###
##########################################################################
    closest_tissues = aggregation$tissue
    closest_samples = rest_tissues %in% closest_tissues

    closest_matrix = rest_matrix[closest_samples,]
    closest_abc    = colmeans(closest_matrix, parallel = FALSE)

    report_data = data.table(gene_id = gene_id, enhancer = all_features, score = tissue_abc, closest = sprintf("%.10f",closest_abc), rest = sprintf("%.10f",complete_abc))
    setorder(report_data,-score)
    report_data$score   = sprintf("%.10f",report_data$score)

    all_features = report_data$enhancer
    max_char = max(nchar(all_features))
    for(index in 1:length(all_features))
    {
	while(nchar(all_features[index]) < max_char)
	{
	    all_features[index] = paste(c(" ", all_features[index]), sep='', collapse='')
	}
    }
    report_data$enhancer = all_features

    fwrite(report_data, file=file_out, sep='\t', quote=F, append=T, col.names=F)
    write("####",file=file_out, append=T)
    write("####",file=log_data, append=T)
}
Sys.sleep(5)
system("gzip -f Pairs.data")
system(sprintf("mv -f Pairs.data.gz %s/temp/batches/%s.%s.pairs.data.gz", data_folder, cur_tissue, task_index))
system(sprintf("mv -f Log.analysis %s/temp/batches/%s.%s.analysis.log", data_folder, cur_tissue, task_index))
warnings()
proc.time()
