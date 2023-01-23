#!/usr/bin/env Rscript

data_folder = "!{data_folder}"
cur_tissue  = "!{tissue}"
task_index  = "!{task_index}"
chromosome  = "!{chromosome}"
genes_list  = "!{list}"
alpha_value = 0
upper_coef  = 0.1

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(glmnet))
options(width=500)
options(scipen=999)

file_out = "Pairs.data"
write(paste(c("human_ensembl","enhancer","impact"), collapse='\t'), file=file_out)

#log_data = "Genes.log"
#write(sprintf("Coeffcient report for %s\n",task_index), file = log_data)

for(gene_id in strsplit(genes_list,";")[[1]])
{
    if(file.exists(sprintf("%s/Glmnet/%s/%s.enhancers.matrix.gz", data_folder, chromosome, gene_id)) == FALSE)
    {
        cat(sprintf("absent %s.enhancers.matrix.gz\n\n", gene_id))
        next
    }
    if(file.exists(sprintf("%s/Glmnet/%s/%s.rotations.matrix.gz", data_folder, chromosome, gene_id)) == FALSE)
    {
        cat(sprintf("absent %s.rotations.matrix.gz\n\n", gene_id))
        next
    }

#######################################################
##### I am loading enhancer and rotation matrices #####
#######################################################
    rotation_data  = fread(file = sprintf("%s/Glmnet/%s/%s.rotations.matrix.gz",data_folder,chromosome,gene_id), sep='\t', header=T, stringsAsFactors=FALSE, blank.lines.skip=T)
    enhancer_list  = rotation_data$enhancer
    rotation_data[,enhancer:=NULL]
    rotation_data  = as.matrix(rotation_data)

    enhancer_data  = fread(file = sprintf("%s/Glmnet/%s/%s.enhancers.matrix.gz",data_folder,chromosome,gene_id), sep='\t', header=T, stringsAsFactors=FALSE, blank.lines.skip=T)
    samples_tissue = enhancer_data$tissue
    enhancer_data[,tissue:=NULL]

    expression = as.numeric(enhancer_data$expression)
    enhancer_data[,expression:=NULL]
    setcolorder(enhancer_data,enhancer_list)

    enhancer_data = as.matrix(enhancer_data)

    cat(sprintf("Do %s gene in %s\n", gene_id, cur_tissue))
    cat("Enhancer matrix:",nrow(enhancer_data),"versus",ncol(enhancer_data),"\n")
    cat("Rotation matrix:",nrow(rotation_data),"versus",ncol(rotation_data),"\n")

##########################################################################################
### I form samples weights vector, five-fold cross-validation folds and transform TPMs ###
##########################################################################################
    tissue_mask   = samples_tissue == cur_tissue
    tissue_mask   = tissue_mask*50/sum(tissue_mask)
    rest_mask     = samples_tissue != cur_tissue
    rest_mask     = rest_mask*50/sum(rest_mask)
    final_weights = tissue_mask + rest_mask
    final_folds   = 1 + (0:(nrow(enhancer_data)-1))%%5
    penalty       = c(rep(1,ncol(rotation_data)),0)

    expression    = 2.763932 + (expression*10)**0.5

###############################################################################################
#### I have to introcude upper limit for ceoefficients, otherwise weak components overgrow ####
###############################################################################################
    upper_limit   = c(rep(upper_coef, ncol(rotation_data)),100)
    lower_limit   = c(rep(0.0,        ncol(rotation_data)),0)

############################################################################
##### Knowing rotation matrix, I get PCA, add Intercept and run Glmnet #####
############################################################################
    Components = enhancer_data %*% rotation_data
    Intercept  = rep(1,nrow(Components))
    Combined   = as.matrix(cbind(Components,Intercept))
    names(Combined)[ncol(Combined)] = "Intrcpt"
    cat("Glmnet matrix:  ",dim(Combined)[1],"versus",dim(Combined)[2],"\n")

    analysis    = cv.glmnet(x = Combined, y = expression, weights = final_weights, family = "gaussian", intercept = F, alpha = alpha_value, standardize = F, foldid = final_folds, grouped=FALSE, nlambda = 300, lower.limit = lower_limit, upper.limit = upper_limit, penalty.factor = penalty)
    best_lambda = analysis$lambda.min
    best_model  = analysis$glmnet.fit
    best_index  = analysis$index[1]
    rsquare     = best_model$dev.ratio[best_index]
    cat("Did Glmnet: lambda =",best_lambda,"rsquare =",rsquare,"\n")

    cs = as.matrix(coef(analysis, s = "lambda.min"))
    print(cs)

#########################################################################
##### Here I calculate enhancer importance. It is as a drop of mean #####
##### predicted expression after setting this enahncer to zero      #####
#########################################################################
    report_data       = data.table(gene_id = gene_id, enhancer = enhancer_list, impact = 0)
    tissue_mask       = samples_tissue == cur_tissue
    tissue_samples    = sum(tissue_mask)
    tissue_data       = enhancer_data[tissue_mask,]
    tissue_expression = expression[tissue_mask]
    tissue_intercept  = Intercept[tissue_mask]

##################################################
## I use complete model predictions as baseline ##
##################################################
    tissue_components  = tissue_data %*% rotation_data
    tissue_combined    = as.matrix(cbind(tissue_components,tissue_intercept))
    prediction_data    = as.numeric(predict(best_model, tissue_combined, s = best_lambda))
    predicted_baseline = mean(prediction_data)

    cat(sprintf("Average %s expression = %.6f, predicted expression = %.6f\n",cur_tissue,mean(tissue_expression),mean(predicted_baseline)))

#################################################################
## I set each enahncer to zero, get PCA components and predict ##
#################################################################
    for(enh_index in 1:nrow(report_data))
    {
	work_matrix     = tissue_data
	work_matrix[,enh_index] = rep(0,tissue_samples)
	work_comp       = work_matrix %*% rotation_data
	work_combined   = as.matrix(cbind(work_comp,tissue_intercept))

	prediction_data = as.numeric(predict(best_model, work_combined, s = best_lambda))
	predicted_mean  = mean(prediction_data)

	if(predicted_mean < 0)
	    predicted_mean = 0

	if(predicted_mean > predicted_baseline)
	    predicted_mean = predicted_baseline

	report_data[enh_index,"impact"] = 1 - predicted_mean/predicted_baseline
    }
    setorder(report_data,-impact)

    report_data = report_data[ impact > 0.00001 ]
    if(nrow(report_data) == 0)
    {
	cat("No enhancers found\n\n")
	next
    }
    report_data$impact = sprintf("%.7f", report_data$impact)
    cat("Found ",nrow(report_data),"enahncers with non-zero impact\n\n")

    fwrite(report_data, file=file_out, sep='\t', quote=F, append=T, col.names=F)
    write("####",file=file_out, append=T)
}
system("gzip -f Pairs.data")
system(sprintf("mv Pairs.data.gz %s/temp/Glmnet/%s.%s.pairs.data.gz", data_folder, cur_tissue, task_index))

warnings()
proc.time()
