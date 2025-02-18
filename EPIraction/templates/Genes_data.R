#!/usr/bin/env Rscript

data_folder   = "!{data_folder}"
chromosome    = "!{chrom}"
genes_list    = "!{list}"
samples_index = "!{samples_index}"

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Rfast))

options(width=180)
options(scipen=999)

expressed_genes = fread(file = sprintf("%s/data/www.Expressed.genes.matrix",data_folder), sep='\t', header=T, stringsAsFactors=FALSE, blank.lines.skip=T, tmpdir="./")
working_data    = fread(file = sprintf("%s/data/yyy.Activity.%s.matrix.gz",data_folder,chromosome), sep='\t', header=T, stringsAsFactors=FALSE, blank.lines.skip=T, tmpdir="./")
contacts_data   = fread(file = sprintf("%s/data/xxx.Contacts.%s.data.gz",data_folder,chromosome), sep='\t', header=T, stringsAsFactors=FALSE, blank.lines.skip=T, tmpdir="./")
samples         = fread(file = sprintf("%s/files/%s",data_folder,samples_index), sep='\t', header=T, stringsAsFactors=FALSE, tmpdir="./")
tissues         = unique(samples$tissue)
samples_tissue  = samples$tissue

cat("Load data\n\n")
for(gene_id in strsplit(genes_list,";")[[1]])
{
    cat("Do",gene_id,"\n")
    if(file.exists(sprintf("%s/Genes/%s/%s.input.matrix.gz", data_folder, chromosome, gene_id)) == TRUE)
    {
	cat(sprintf("already %s/Genes/%s/%s.input.matrix.gz\n\n", data_folder, chromosome, gene_id))
	next
    }

##########################################################################
##### I get contact-per-sample matrix from contact-per-tissue matrix #####
##########################################################################
    contact_matrix = contacts_data[promoter == gene_id]
    enhancer_list  = contact_matrix$enhancer
    contact_matrix[,promoter:=NULL]
    contact_matrix[,distance:=NULL]
    contact_transpose = data.table::transpose(contact_matrix, make.names = "enhancer", keep.names = "tissue")

    contact_matrix  = merge(samples, contact_transpose, by = "tissue", sort = F)
    setorder(contact_matrix, "index")
    contact_matrix[,tissue:=NULL]
    contact_matrix[,index:=NULL]
    contact_matrix[,accession:=NULL]
    contact_matrix[,file:=NULL]
    contact_matrix[,title:=NULL]
    contact_matrix[,signal:=NULL]
    contact_matrix[,input:=NULL]
    setcolorder(contact_matrix,enhancer_list)

###############################################################
##### I load Activity matrix and multiply it by contacts. #####
###############################################################
    enhancer_data   = working_data[region %in% enhancer_list]
    enhancer_matrix = data.table::transpose(enhancer_data, make.names = "region", keep.names = "samples")
    enhancer_matrix[,samples:=NULL]
    setcolorder(enhancer_matrix,enhancer_list)
    enhancer_matrix = enhancer_matrix*contact_matrix
    enhancer_matrix = as.matrix(enhancer_matrix)
    cat("Complete matrix:   ",nrow(enhancer_matrix),"versus",ncol(enhancer_matrix),"\n")
    enhancer_list = colnames(enhancer_matrix)

######################################################################################
##### I select nonzero features. I need at least two columns to make matrix in R #####
######################################################################################
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

###########################################
### I add expression to features matrix ###
###########################################
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

    enhancer_matrix = cbind(expression = expression,enhancer_matrix)

###########
## Done ###
###########
    report_folder = sprintf("%s/Genes/%s", data_folder, chromosome)

    report_data = cbind(tissue = samples_expression$tissue, enhancer_matrix)
    fwrite(as.data.table(report_data), file = sprintf("%s/%s.input.matrix.gz", report_folder, gene_id), sep = '\t', quote=F, compress="gzip")
    cat("Done ",gene_id,"\n\n")
}
proc.time()
