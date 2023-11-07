// Use DSL2
nextflow.enable.dsl=2

data_folder   = params.data_folder
tissues_index = params.tissues_index
samples_index = params.samples_index
regions       = params.regions
chromosomes   = channel.of("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX").toList()
chrom_batches = channel.of("chr1, chr13, chr17, chr19", "chr14, chr2, chrX", "chr16, chr21, chr3, chr6", "chr12, chr4, chr5", "chr11, chr20, chr7, chr8", "chr10, chr15, chr18, chr22, chr9").toList()

transform_aaa = 20.94289549
transform_bbb = 70.34457472

annotation_done = "yes"
normalize_done  = "yes"
activity_done   = "yes"

samples_done    = "yes"
samples_merge   = "yes"
pairs_done      = "yes"

hic_consensus   = "yes"
hic_done        = "yes"

first_qc        = '/dev/null'

contacts_done   = "yes"
contacts_merge  = "yes"

data_status     = "yes"
glmnet_done     = "yes"

// Print log information
log.info """
Welcome to EPIraction 1.2, nextflow implementation

data_folder : $data_folder
tissues     : $tissues_index
samples     : $samples_index
regions     : $regions
projectDir  : $projectDir
"""

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**** Working on consensus Open, Cofactor and H3K27ac to develop quantile-normalizes consensus activities for each tissue ****/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * This process takes tissue data from tissues_index file and annotate regions
 * with molecular information from corresponding BigWig files and RNA-seq data
 */
process Regions_annotate
{
    cache false

    input:
	val data_folder
	val regions
	tuple val(index), val(tissue), val(RNA_seq), val(Open), val(Cofactor), val(HiC)

    output:
	val tissue

    shell:
	data_folder = "$data_folder"
	regions     = "$regions"
	tissue      = "$tissue"
	RNA_seq     = "$RNA_seq"
	Open        = "$Open"
	Cofactor    = "$Cofactor"
	transform_aaa = "$transform_aaa"
	transform_bbb = "$transform_bbb"
	template "Regions_annotate.pl"
}
/*
 * This process quantile-normalizes Activity data
 */
process Regions_normalize
{
    cache false

    input:
	val data_folder
	val tissues
	val status

    output:
	val "done"

    shell:
	data_folder = "$data_folder"
	tissues     = "$tissues"
	template "Regions_normalize.R"
}

/*
 * This process wrote quantile-normalizes activities for enhancers
 */
process Regions_final
{
    cache false

    input:
	val data_folder
	tuple val(index), val(tissue), val(RNA_seq), val(Open), val(Cofactor), val(HiC)
	val status

    output:
	val tissue

    shell:
	data_folder = "$data_folder"
	tissue      = "$tissue"
	template "Regions_final.R"
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**** Using quantile-normalizes consensus activities to assign sample-specific activities for each sample and make full matrices ****/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * This process makes Activity matrix for every individual sample within tissue by chromosome in temp folder
 */
process Samples_activity
{
    cache false

    input:
	val data_folder
	val samples
	tuple val(tissue), val(chromosomes)
	val status

    output:
	val tissue

    shell:
	data_folder   = "$data_folder"
	samples_index = "$samples"
	regions       = "$regions"
	tissue        = "$tissue"
	chromosomes   = "$chromosomes"
	template "Samples_activity.pl"
}
/*
 * This process merges all sample-specific activity data by chromosome
 */
process Samples_merge
{
    cache false

    input:
	val data_folder
	val chrom
	val tissues
	val status

    output:
	val chrom

    shell:
	data_folder = "$data_folder"
	chrom       = "$chrom"
	tissues     = "$tissues"
	template "Samples_merge.R"
}

/*
 * This process generates promoter-enhancer pairs
 */
process Make_pairs
{
    cache false

    input:
	val data_folder
	val regions
	val chrom

    output:
	val chrom

    shell:
	data_folder = "$data_folder"
	regions     = "$regions"
	chrom       = "$chrom"
	template "Pairs.pl"
}

///////////////////////////////////////////////////////////////
/**** The following processes dump and normalize HiC data ****/
///////////////////////////////////////////////////////////////
/*
 * These two processes dumps HiC data. For all promoters versus all 2000 nt bins
 */
process HiC_consensus
{
    cache false

    input:
	val data_folder
	val regions
	val chrom

    output:
	val "done"

    shell:
	data_folder  = "$data_folder"
	chrom        = "$chrom"
	regions      = "$regions"
	genome_file  = "$projectDir/files/Human.bed"
	blanks_pl    = "$projectDir/bin/blanks.pl"

	template "HiC_consensus.pl"
}
process HiC_dump
{
    cache false

    input:
	val data_folder
	val regions
	tuple val(sample), val(chrom)
	val status

    output:
	val v_return

    shell:
	data_folder  = "$data_folder"
	sample       = "$sample"
	chrom        = "$chrom"
	regions      = "$regions"
	v_return     = "$sample.$chrom"
	genome_file  = "$projectDir/files/Human.bed"
	blanks_pl    = "$projectDir/bin/blanks.pl"

	template "HiC_dump.pl"
}


/*
 * This is first QC. Need region annotation for each tissue, all promoter-enhancer pairs, dump and merge Activity data and HiC dump
 */
process First_qc
{
    cache false

    input:
	val data_folder
	val tissues
	val hic_samples
	val chromosomes
	val activity_done
	val samples_merge
	val pairs_done
	val hic_done
	val hic_consensus

    output:
	path "First_QC.log"

    shell:
	data_folder = "$data_folder"
	tissues     = "$tissues"
	hic_samples = "$hic_samples"
	chromosomes = "$chromosomes"
	template "First_qc.pl"
}

/***********------------------********  Second Block ******------------------------*******/
//////////////////////////////////////////////////////////////////////////////////////////////
/**** Working on Activity and Hi-C contact data to develop their product for each tissue ****/
//////////////////////////////////////////////////////////////////////////////////////////////
/*
 * This process annotates all promoter-enhancer pairs with HiC data for every tissue
 */
process Contacts_tissues
{
    cache false

    input:
	val data_folder
	tuple val(index), val(tissue), val(RNA_seq), val(Open), val(Cofactor), val(HiC)
	val chromosomes
	path qc_first

    output:
	val tissue

    shell:
	data_folder  = "$data_folder"
	tissue       = "$tissue"
	regions      = "$regions"
	HiC_samples  = "$HiC"
	chromosomes  = "$chromosomes"
	HiC_baseline = "$projectDir/files/HiC.baseline"
	upstream     = "$projectDir/files/Upstream.smooth"
	downstream   = "$projectDir/files/Downstream.smooth"
	template "Contacts_tissues.pl"
}
/*
 * This process aggregates tissue-spcific contact data into matrices for each chromosome
 */
process Contacts_merge
{
    cache false

    input:
	val data_folder
	val tissues
	val regions
	val chrom
	val t_status

    output:
	val chrom

    shell:
	data_folder = "$data_folder"
	tissues     = "$tissues"
	chrom       = "$chrom"
	regions     = "$regions"
	template "Contacts_merge.pl"
}

//////////////////////////////////////////////////////////////////////////////
/**** Use of "Activity by Contact data" to develop input data for Glmnet ****/
//////////////////////////////////////////////////////////////////////////////
/*
 * This process scans "$tissue.regions.data" files, aggregates all expressed genes,
 * splits this list by chromosome into the blocks and returns them as a file
 */
process Glmnet_expressed
{
    cache false

    input:
	val data_folder
	val tissues
	val status

    output:
	path "Expressed.genes.tsv"

    shell:
	data_folder = "$data_folder"
	tissues     = "$tissues"
	blocks      = 300
	template "Glmnet_expressed.pl"
}

/*
 * This process takes block of genes generated by "Genes_expressed". For every gene it generates ABC matrix,
 * simplifies this matrix to avoid small values and generates "$gene_id.input.matrix.gz" for linear modelling
*/
process Glmnet_data
{
    cache false

    input:
	val data_folder
	val samples_index
	tuple val(group), val(list), val(score)

    output:
	val "done"

    shell:
	data_folder   = "$data_folder"
	chrom         = "$group"
	list          = "$list"
	samples_index = "$samples_index"
	template "Glmnet_data.R"
}

/*
 * This process scans all "$tissue.regions.data" files one by one. For every tissue it obtains 
 * all expressed genes, splits this list into blocks and returns this list to "stdout"
 */
process Glmnet_tissues
{
    cache false

    input:
	val data_folder
	val tissues
	val chromosomes
	val status

    output:
	path "Tissue.genes.tsv"

    shell:
	data_folder = "$data_folder"
	tissues     = "$tissues"
	blocks      = 100
	chromosomes = "$chromosomes"
	template "Glmnet_tissues.pl"
}

/*
 * This process runs Glmnet
 */
process Glmnet_run
{
    cache false

    input:
	val data_folder
	tuple val(task_index), val(tissue), val(chromosome), val(list)

    output:
	val "done"

    shell:
	data_folder   = "$data_folder"
	tissue        = "$tissue"
	task_index    = "$task_index"
	chromosome    = "$chromosome"
	list          = "$list"
	transform_aaa = "$transform_aaa"
	transform_bbb = "$transform_bbb"
	template "Glmnet_run.R"
}

/*
 * This process agregates predicted promoter-enhancer pairs and makes final report data
 */
process Report_pairs
{
    cache false

    input:
	val data_folder
	tuple val(index), val(tissue), val(RNA_seq), val(Open), val(Cofactor), val(HiC)
	val status

    output:
	val tissue

    shell:
	data_folder  = "$data_folder"
	tissue       = "$tissue"
	genome_file  = "$projectDir/files/Human.genome"
	interact     = "$projectDir/files/interact.as"
	HiC_baseline = "$projectDir/files/HiC.baseline"
	template "Report.pl"
}

workflow
{
    indexData   = channel.fromPath(data_folder+"/files/"+tissues_index) | splitCsv(header:true, sep: '\t') | toList
    All_tissues = indexData.flatMap{ row -> row.tissue }                                                   | toList
    HiC_samples = indexData.flatMap{ row -> row.HiC    } | flatMap{  it -> it.split(";") } | unique        | toList

//#####################################################################################
//### Here I normalize Open chromatin, Cofactor and H3K27ac at the level of tissues ###
//#####################################################################################
    annotation_done = Regions_annotate(data_folder,regions,indexData.flatMap())     | collect
    normalize_done  = Regions_normalize(data_folder,All_tissues,annotation_done)    | collect
    activity_done   = Regions_final(data_folder,indexData.flatMap(),normalize_done) | collect

//#############################################################################
//### Here I write Activities for each enhancer in each individual H3K27ac  ###
//#############################################################################
    samples_done  = Samples_activity(data_folder,samples_index,All_tissues.flatMap().combine(chrom_batches.flatMap()),activity_done) | collect
    samples_merge = Samples_merge(data_folder,chromosomes.flatMap(),All_tissues,samples_done)                                        | collect

//####################################################################
//### Here I write all gene-enhancer pairs, dump Hi-C data and QC  ###
//####################################################################
    hic_consensus  = HiC_consensus(data_folder,regions,chromosomes.flatMap())                                         | collect
    hic_done       = HiC_dump(data_folder,regions,HiC_samples.flatMap().combine(chromosomes.flatMap()),hic_consensus) | collect

    pairs_done    = Make_pairs(data_folder,regions,chromosomes.flatMap()) | collect
    first_qc      = First_qc(data_folder, All_tissues, HiC_samples, chromosomes, activity_done, samples_merge, pairs_done, hic_done, hic_consensus) | collect

//########################################################################
//### Here I form tissue-specific Contact and Activity-by-Contact data ###
//########################################################################
    contacts_done  = Contacts_tissues(data_folder,indexData.flatMap(),chromosomes,first_qc)              | collect
    contacts_merge = Contacts_merge(data_folder,All_tissues,regions,chromosomes.flatMap(),contacts_done) | collect

//##############################################
//### Here I form Glmnet data and run Glmnet ###
//##############################################
    All_genes   = Glmnet_expressed(data_folder,All_tissues,contacts_merge) | splitCsv(header:true, sep: '\t')           | toList
    data_status = Glmnet_data(data_folder, samples_index, All_genes.flatMap())                                          | collect
    Tis_genes   = Glmnet_tissues(data_folder, All_tissues, chromosomes, data_status) | splitCsv(header:true, sep: '\t') | toList
    glmnet_done = Glmnet_run(data_folder,Tis_genes.flatMap())               | collect

//##############################
//### Here I report the data ###
//##############################
    report_done = Report_pairs(data_folder,indexData.flatMap(),glmnet_done) | collect
}
