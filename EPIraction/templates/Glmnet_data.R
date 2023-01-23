#!/usr/bin/env Rscript

data_folder   = "!{data_folder}"
chromosome    = "!{chrom}"
genes_list    = "!{list}"
samples_index = "!{samples_index}"

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Rfast))
suppressPackageStartupMessages(library(nsprcomp))

options(width=180)
options(scipen=999)
minimum_signal    = 0.02
nsprcomp_tol      = 0.01

expressed_genes = fread(file = sprintf("%s/data/yyy.Expressed.genes.matrix",data_folder), sep='\t', header=T, stringsAsFactors=FALSE, blank.lines.skip=T)
working_data    = fread(file = sprintf("%s/data/yyy.H3K27ac.%s.matrix.gz",data_folder,chromosome), sep='\t', header=T, stringsAsFactors=FALSE, blank.lines.skip=T)
scaling_data    = fread(file = sprintf("%s/data/xxx.Scalings.%s.data.gz",data_folder,chromosome), sep='\t', header=T, stringsAsFactors=FALSE, blank.lines.skip=T)
samples         = fread(file = sprintf("%s/files/%s",data_folder,samples_index), sep='\t', header=T, stringsAsFactors=FALSE)
tissues         = unique(samples$tissue)

cat("Load data\n\n")

for(gene_id in strsplit(genes_list,";")[[1]])
{
    cat("Do",gene_id,"\n")

##########################################################################
##### I get scaling-per-sample matrix from scaling-per-tissue matrix #####
##########################################################################
    scaling_matrix = scaling_data[promoter == gene_id]
    enhancer_list  = scaling_matrix$enhancer
    scaling_matrix[,promoter:=NULL]
    scaling_matrix[,distance:=NULL]
    scaling_transpose = data.table::transpose(scaling_matrix, make.names = "enhancer", keep.names = "tissue")

    scaling_matrix = merge(samples, scaling_transpose, by = "tissue", sort = F)
    setorder(scaling_matrix, "index")
    scaling_matrix[,tissue:=NULL]
    scaling_matrix[,group:=NULL]
    scaling_matrix[,index:=NULL]
    scaling_matrix[,accession:=NULL]
    scaling_matrix[,file:=NULL]
    scaling_matrix[,title:=NULL]
    scaling_matrix[,signal:=NULL]
    scaling_matrix[,input:=NULL]
    setcolorder(scaling_matrix,enhancer_list)

#######################################################################
##### I load H3K27ac enhancer matrix and multiply it by scaligns. #####
#######################################################################
    enhancer_data   = working_data[region %in% enhancer_list]
    enhancer_matrix = data.table::transpose(enhancer_data, make.names = "region", keep.names = "samples")
    enhancer_matrix[,samples:=NULL]
    setcolorder(enhancer_matrix,enhancer_list)
    enhancer_matrix = enhancer_matrix*scaling_matrix
    enhancer_matrix = as.matrix(enhancer_matrix)
    cat("Complete matrix:   ",nrow(enhancer_matrix),"versus",ncol(enhancer_matrix),"\n")
    enhancer_list = colnames(enhancer_matrix)

###################################################################################################
##### I select nonzero enhancers. I need at least two columns (enahncers) to make matrix in R #####
###################################################################################################
    col_means = Rfast::colmeans(enhancer_matrix, parallel = FALSE)
    nonzero   = col_means > 0
    if(sum(nonzero) <= 1)
    {
	cat("No valid enhancers\n")
	next
    }
    enhancer_matrix = enhancer_matrix[,nonzero]
    cat("Nonzero matrix:    ",nrow(enhancer_matrix),"versus",ncol(enhancer_matrix),"\n")
    enhancer_list = colnames(enhancer_matrix)

#################################################
###  I select only substantial enahncers here ###
#################################################
    for(enh_index in 1:ncol(enhancer_matrix))
    {
	enhancer_values = enhancer_matrix[,enh_index]

	zeroing_values  = vector()
	for(tissue in tissues)
	{
	    mask = samples$tissue == tissue
	    tissue_enhancers = enhancer_values[mask]
	    if(mean(tissue_enhancers) < minimum_signal)
	    {
		zeroing_values = c(zeroing_values,rep(0,length(tissue_enhancers)))
	    }
	    else {
		zeroing_values = c(zeroing_values,tissue_enhancers)
	    }
	}
	enhancer_matrix[,enh_index] = zeroing_values
    }
    col_means = Rfast::colmeans(enhancer_matrix, parallel = FALSE)
    nonzero   = col_means > 0
    if(sum(nonzero) <= 1)
    {
	cat("No valid enhancers\n")
	next
    }
    enhancer_matrix = enhancer_matrix[,nonzero]
    cat("Substantial matrix:",nrow(enhancer_matrix),"versus",ncol(enhancer_matrix),"\n")
    enhancer_list = colnames(enhancer_matrix)

#########################################################################################################
### I add promoter H3K27ac to the features matrix. It helps a lot when no strong enhancers exist. In  ###
### this case promoter forms the most variable first PCA component, and give a frame for variation of ###
### weak enahncers. If strong enhancers exist, non-negative PCA merge them with promoter values.      ###
#########################################################################################################
    promoter_data = working_data[region == gene_id]
    promoter_data[,region:=NULL]
    promoter_data = as.numeric(promoter_data)
    if(mean(promoter_data) == 0)
    {
	cat("No promoter signal\n")
	next
    }
    enhancer_matrix = cbind(promoter_data,enhancer_matrix)
    colnames(enhancer_matrix)[1] = "pseudo_promoter_feature_id"

####################################
###  I run Nonnegative PCA here  ###
####################################
    set.seed(127)
    analysis  = nsprcomp( x = enhancer_matrix, nneg = TRUE, scale. = FALSE, center = FALSE, nrestart = 5, em_tol = 0.001, em_maxiter = 100000, verbosity = 0, tol = nsprcomp_tol)
    Rotations = analysis$rotation

    for(enh_index in 1:nrow(Rotations))
    {
	enhancer_row = Rotations[enh_index,]
	if(sum(enhancer_row) > 1)
	{
	    Rotations[enh_index,] = enhancer_row/sum(enhancer_row)
	}
    }
    Rotations = apply(Rotations, 1:2, function(x) sprintf("%.8f", x))

#################################
###  I from expression vector ###
#################################
    expressed = expressed_genes[ensembl_gene == gene_id]
    if(nrow(expressed) != 1)
    {
	cat("Not expressed gene",gene_id,"\n")
	next
    }
    tissue_expression = data.table(tissue = colnames(expressed)[-c(1,2,3)], expression = sprintf("%.7f",(as.vector(expressed))[-c(1,2,3)]))
    samples_expression = merge(samples, tissue_expression, by = "tissue", sort = F)
    setorder(samples_expression, "index")

###########
## Done ###
###########
    report_folder = sprintf("%s/Glmnet/%s", data_folder, chromosome)

    report_data = cbind(enhancer = rownames(Rotations), Rotations)
    fwrite(as.data.table(report_data), file = sprintf("%s/%s.rotations.matrix.gz", report_folder, gene_id), sep='\t', quote=F, col.names=T, compress="gzip")

    report_data = cbind(tissue = samples_expression$tissue, expression = samples_expression$expression, enhancer_matrix)
    fwrite(as.data.table(report_data), file = sprintf("%s/%s.enhancers.matrix.gz", report_folder, gene_id), sep = '\t', quote=F, compress="gzip")
    cat("Done ",gene_id,"\n\n")
}
proc.time()
