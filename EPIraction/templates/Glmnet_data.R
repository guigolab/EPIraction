#!/usr/bin/env Rscript

data_folder   = "!{data_folder}"
chromosome    = "!{chrom}"
genes_list    = "!{list}"
samples_index = "!{samples_index}"

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Rfast))

options(width=180)
options(scipen=999)
minimum_signal  = 0.05

expressed_genes = fread(file = sprintf("%s/data/yyy.Expressed.genes.matrix",data_folder), sep='\t', header=T, stringsAsFactors=FALSE, blank.lines.skip=T, tmpdir="./")
working_data    = fread(file = sprintf("%s/data/yyy.Activity.%s.matrix.gz",data_folder,chromosome), sep='\t', header=T, stringsAsFactors=FALSE, blank.lines.skip=T, tmpdir="./")
contacts_data   = fread(file = sprintf("%s/data/xxx.Contacts.%s.data.gz",data_folder,chromosome), sep='\t', header=T, stringsAsFactors=FALSE, blank.lines.skip=T, tmpdir="./")
samples         = fread(file = sprintf("%s/files/%s",data_folder,samples_index), sep='\t', header=T, stringsAsFactors=FALSE, tmpdir="./")
tissues         = unique(samples$tissue)

cat("Load data\n\n")

for(gene_id in strsplit(genes_list,";")[[1]])
{
    cat("Do",gene_id,"\n")
    if(file.exists(sprintf("%s/Glmnet/%s/%s.enhancers.matrix.gz", data_folder, chromosome, gene_id)) == TRUE)
    {
	cat(sprintf("already %s/Glmnet/%s/%s.enhancers.matrix.gz\n\n", data_folder, chromosome, gene_id))
#	next
    }

##########################################################################
##### I get contact-per-sample matrix from contact-per-tissue matrix #####
##########################################################################
    contact_matrix = contacts_data[promoter == gene_id]
    enhancer_list  = contact_matrix$enhancer
    contact_matrix[,promoter:=NULL]
    contact_matrix[,distance:=NULL]
    contact_transpose = data.table::transpose(contact_matrix, make.names = "enhancer", keep.names = "tissue")

    contact_matrix = merge(samples, contact_transpose, by = "tissue", sort = F)
    setorder(contact_matrix, "index")
    contact_matrix[,tissue:=NULL]
    contact_matrix[,index:=NULL]
    contact_matrix[,accession:=NULL]
    contact_matrix[,file:=NULL]
    contact_matrix[,title:=NULL]
    contact_matrix[,signal:=NULL]
    contact_matrix[,input:=NULL]
    setcolorder(contact_matrix,enhancer_list)

#######################################################################
##### I load H3K27ac enhancer matrix and multiply it by scaligns. #####
#######################################################################
    enhancer_data   = working_data[region %in% enhancer_list]
    enhancer_matrix = data.table::transpose(enhancer_data, make.names = "region", keep.names = "samples")
    enhancer_matrix[,samples:=NULL]
    setcolorder(enhancer_matrix,enhancer_list)
    enhancer_matrix = enhancer_matrix*contact_matrix
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
###  I select only substantial enhancers here ###
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

###########################################################################
### I add promoter H3K27ac and expression to enhancer matrix for report ###
###########################################################################
    promoter_data = working_data[region == gene_id]
    promoter_data[,region:=NULL]
    promoter_data = as.numeric(promoter_data)
    if(mean(promoter_data) == 0)
    {
	cat("No promoter signal\n")
	next
    }
    enhancer_matrix = cbind(promoter_data,enhancer_matrix)
    colnames(enhancer_matrix)[1] = "promoter_feature"

    expressed = expressed_genes[ensembl_gene == gene_id]
    if(nrow(expressed) != 1)
    {
	cat("Not expressed gene",gene_id,"\n")
	next
    }

    tissue_expression  = data.table(tissue = colnames(expressed)[-c(1,2,3)], expression = as.vector(expressed)[-c(1,2,3)])
    samples_expression = merge(samples, tissue_expression, by = "tissue", sort = F)
    setorder(samples_expression, "index")
    expression  = as.numeric(samples_expression$expression)

    enhancer_matrix = cbind(expression,enhancer_matrix)
    colnames(enhancer_matrix)[1] = "expression"

###########
## Done ###
###########
    report_folder = sprintf("%s/Glmnet/%s", data_folder, chromosome)

    report_data = cbind(tissue = samples_expression$tissue, enhancer_matrix)
    fwrite(as.data.table(report_data), file = sprintf("%s/%s.input.matrix.gz", report_folder, gene_id), sep = '\t', quote=F, compress="gzip")
    cat("Done ",gene_id,"\n\n")
}
proc.time()
