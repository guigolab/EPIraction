#!/usr/bin/env Rscript

data_folder  = "!{data_folder}"
cur_tissue   = "!{tissue}"
task_index   = "!{task_index}"
chromosome   = "!{chromosome}"
genes_list   = "!{list}"
aaaaa        = !{transform_aaa}
bbbbb        = !{transform_bbb}

alpha_value = 0
upper_coef  = 1
pseudocount = 0.01

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(Rfast))
options(width=10000)
options(scipen=999)

if(file.exists(sprintf("%s/temp/Glmnet/%s.%s.glmnet.log", data_folder, cur_tissue, task_index)) == TRUE)
{
    cat(sprintf("Already done %s.%s.glmnet.log\n", cur_tissue, task_index))
    q()
}


log_transformation = function(TPM_value)
{
    log_expression = log2(pseudocount + TPM_value)*aaaaa + bbbbb
    if(log_expression < 0)
	return(0)
    return(log_expression)
}

file_out = "Pairs.data"
write(paste(c("human_ensembl","enhancer","impact"), collapse='\t'), file=file_out)

log_data = "Log.glmnet"
write(sprintf("Glmnet report for %s\n",task_index), file = log_data)

for(gene_id in strsplit(genes_list,";")[[1]])
{
    if(file.exists(sprintf("%s/Glmnet/%s/%s.input.matrix.gz", data_folder, chromosome, gene_id)) == FALSE)
    {
        cat(sprintf("absent %s.enhancers.input.gz\n\n", gene_id))
        next
    }
    cat(sprintf("Do %s gene in %s\n\n", gene_id, cur_tissue))
    write(sprintf("Do %s gene in %s", gene_id, cur_tissue),file=log_data, append=T)

##########################################
##### I am loading the training data #####
##########################################
    enhancer_data  = fread(file = sprintf("%s/Glmnet/%s/%s.input.matrix.gz",data_folder,chromosome,gene_id), sep='\t', header=T, stringsAsFactors=FALSE, blank.lines.skip=T)
    samples_tissue = enhancer_data$tissue
    enhancer_data[,tissue:=NULL]

    expression = as.numeric(enhancer_data$expression)
    enhancer_data[,expression:=NULL]

    all_features  = colnames(enhancer_data)
    enhancer_data = as.matrix(enhancer_data)
    colnames(enhancer_data) = all_features

##########################################################################################
### I form samples weights vector, five-fold cross-validation folds and transform TPMs ###
##########################################################################################
    tissue_mask   = samples_tissue == cur_tissue
    tissue_mask   = tissue_mask*50/sum(tissue_mask)
    rest_mask     = samples_tissue != cur_tissue
    rest_mask     = rest_mask*50/sum(rest_mask)
    final_weights = tissue_mask + rest_mask
    final_folds   = 1 + (0:(nrow(enhancer_data)-1))%%5
    response      = sapply(expression,log_transformation)

    tissue_mask     = samples_tissue == cur_tissue
    tissue_samples  = sum(tissue_mask)
    tissue_data     = enhancer_data[tissue_mask,]
    tissue_response = response[tissue_mask]

    average_tissue_activity = colmeans(tissue_data, parallel = FALSE)
    average_response        = mean(tissue_response)
    sum_tissue_activity     = sum(average_tissue_activity)
#    cat(sprintf("Transformed Expression = %.6f\n",average_response))
#    cat(sprintf("Sum of all activities  = %.6f\n",sum_tissue_activity))
    if(average_response < sum_tissue_activity)
    {
	correction = sum_tissue_activity/average_response
	cat(  sprintf("Response correction = %.6f\n",correction))
	write(sprintf("Response correction = %.6f",correction),file=log_data, append=T)
	response = response*correction
	tissue_response = response[tissue_mask]
    }

##################################################################################################
#### I introcude upper and lower limits for ceoefficients, otherwise weak components overgrow ####
##################################################################################################
    upper_limit   = rep(upper_coef, ncol(enhancer_data))
    upper_limit[1] = 1
    lower_limit   = rep(0,          ncol(enhancer_data))

##########################
##### Here I glmnet  #####
##########################
    analysis    = cv.glmnet(x = enhancer_data, y = response, weights = final_weights, family = "gaussian", intercept = F, alpha = alpha_value, standardize = F, foldid = final_folds, grouped=FALSE, nlambda = 300, lower.limit = lower_limit, upper.limit = upper_limit)
    best_lambda = analysis$lambda.min
    best_model  = analysis$glmnet.fit
    best_index  = analysis$index[1]
    rsquare     = best_model$dev.ratio[best_index]

#    cs = predict(analysis, s = best_lambda, type="coefficients")
#    print(cs)
#    fwrite(data.table(variable = rownames(cs), formatSpMatrix(cs, digits=10, col.names=T, zero.print="0.0", align = "fancy")), file = log_data, append = T, sep = "\t", quote=F, col.names=T)
#    write("#",file=log_data, append=T)

#######################################################################################################################################
##### Here I calculate enhancer importance. It is as a drop of mean predicted response value after setting this enahncer to zero  #####
#######################################################################################################################################
    report_data = data.table(gene_id = gene_id, enhancer = all_features, impact = 0)

################################################################################
## I use complete model to predict the gene expression to form baseline value ##
################################################################################
    prediction_data    = as.numeric(predict(best_model, tissue_data, s = best_lambda))
    predicted_baseline = mean(prediction_data)

    cat(sprintf("Average %s response = %.6f, predicted response = %.6f\n",cur_tissue,mean(tissue_response),mean(predicted_baseline)))
    write(sprintf("\nAverage %s response = %.6f, predicted response = %.6f",cur_tissue,mean(tissue_response),mean(predicted_baseline)),file=log_data, append=T)

##########################################################
## I set promoter and each enhancer to zero and predict ##
##########################################################
    work_matrix = tissue_data
    work_matrix[,1] = rep(0,tissue_samples)

    prediction_data = as.numeric(predict(best_model, work_matrix, s = best_lambda))
    predicted_mean  = mean(prediction_data)

    if(predicted_mean < 0)
    predicted_mean = 0
    if(predicted_mean > predicted_baseline)
	predicted_mean = predicted_baseline
    impact = 1 - predicted_mean/predicted_baseline
    cat(sprintf("Promoter impact is %.6f\n",impact))

    write(sprintf("Promoter impact is %.6f\n##",impact), file=log_data, append=T)
    write(sprintf("%s\tpromoter_feature_data\t%.7f",gene_id,impact), file=file_out, append=T)

    for(enh_index in 2:nrow(report_data))
    {
	work_matrix     = tissue_data
	work_matrix[,enh_index] = rep(0,tissue_samples)
	prediction_data = as.numeric(predict(best_model, work_matrix, s = best_lambda))
	predicted_mean  = mean(prediction_data)

	if(predicted_mean < 0)
	    predicted_mean = 0

	if(predicted_mean > predicted_baseline)
	    predicted_mean = predicted_baseline

	report_data[enh_index,"impact"] = 1 - predicted_mean/predicted_baseline
    }
    setorder(report_data,-impact)

    report_data = report_data[ impact > 0.000001 ]
    if(nrow(report_data) == 0)
    {
	cat("No enhancers found\n\n")
	next
    }
    report_data$impact = sprintf("%.7f", report_data$impact)
    cat("Found",nrow(report_data),"enahncers with non-zero impact\n\n")

    fwrite(report_data, file=file_out, sep='\t', quote=F, append=T, col.names=F)
    write("####",file=file_out, append=T)
}
system("gzip -f Pairs.data")
system(sprintf("mv Pairs.data.gz %s/temp/Glmnet/%s.%s.pairs.data.gz", data_folder, cur_tissue, task_index))
system(sprintf("mv Log.glmnet %s/temp/Glmnet/%s.%s.glmnet.log", data_folder, cur_tissue, task_index))
warnings()
proc.time()
