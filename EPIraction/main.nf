// Use DSL2
nextflow.enable.dsl=2

data_folder   = params.data_folder
tissues_index = params.tissues_index
samples_index = params.samples_index
regions       = params.regions
chromosomes   = channel.of("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX").toList()
chrom_batches = channel.of("chr1, chr13, chr17, chr19", "chr14, chr2, chrX", "chr16, chr21, chr3, chr6", "chr12, chr4, chr5", "chr11, chr20, chr7, chr8", "chr10, chr15, chr18, chr22, chr9").toList()

////////////////////////////////////////
/* Pipeline workflow status variables */
////////////////////////////////////////
minimum_expression = 1
genes_blocks       = 100
analysis_blocks    = 300
version            = "v1.6"

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
analysis_done   = "yes"


// Print log information
log.info """
Welcome to EPIraction 1.6, nextflow implementation

data_folder : $data_folder
tissues     : $tissues_index
samples     : $samples_index
regions     : $regions
projectDir  : $projectDir
"""

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**** Working on tissue-consensus Open, Cofactor and H3K27ac to develop tissue-consensus activities for each tissue ****/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
	tuple val(index), val(tissue), val(short), val(RNA_seq), val(CTCF), val(Open), val(Cofactor), val(HiC), val(Open_min), val(H3K27ac_min), val(CTCF_min)

    output:
	val tissue

    shell:
	data_folder = "$data_folder"
	regions     = "$regions"
	tissue      = "$tissue"
	RNA_seq     = "$RNA_seq"
	Open        = "$Open"
	Cofactor    = "$Cofactor"
	CTCF        = "$CTCF"
	open_min    = "$Open_min"
	H3K27ac_min = "$H3K27ac_min"
	CTCFc_min   = "$CTCF_min"
	genome_file = "$projectDir/files/Human.bed"
	expression  = "$minimum_expression"
	template "Regions_annotate.pl"
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**** Using tissue-consensus activities to assign sample-specific activities for each sample and make full matrices ****/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
 * These four processes dumps HiC data. For all promoters versus all 2000 nt bins
 */
process HiC_consensus_dump
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

	template "HiC_consensus_dump.pl"
}
process HiC_consensus_contacts
{
    cache false

    input:
	val data_folder
	val regions
	val chromosomes
	val hic_dump

    output:
	val "done"

    shell:
	data_folder  = "$data_folder"
	regions      = "$regions"
	genome_file  = "$projectDir/files/Human.bed"
	HiC_baseline = "$projectDir/files/HiC.baseline"
	chromosomes  = "$chromosomes"

	template "HiC_consensus_contacts.pl"
}

process HiC_dump
{
    cache false

    input:
	val data_folder
	val regions
	tuple val(sample), val(chrom)

    output:
	val v_return

    shell:
	data_folder  = "$data_folder"
	sample       = "$sample"
	chrom        = "$chrom"
	regions      = "$regions"
	v_return     = "$sample.$chrom"
	genome_file  = "$projectDir/files/Human.bed"

	template "HiC_dump.pl"
}
process HiC_contacts
{
    cache false

    input:
	val data_folder
	val regions
	val sample
	val chromosomes
	val hic_dump

    output:
	val "done"

    shell:
	data_folder  = "$data_folder"
	sample       = "$sample"
	regions      = "$regions"
	genome_file  = "$projectDir/files/Human.bed"
	HiC_baseline = "$projectDir/files/HiC.baseline"
	chromosomes  = "$chromosomes"

	template "HiC_contacts.pl"
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
	val hic_contacts
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
	tuple val(index), val(tissue), val(short), val(RNA_seq), val(CTCF), val(Open), val(Cofactor), val(HiC), val(Open_min), val(H3K27ac_min), val(CTCF_min)
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
	HiC_upper    = "$projectDir/files/HiC.upper"
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

/////////////////////////////////////////////////////////////////////////////////
/**** Use of Activity_by_Contact data to develop input data for predictions ****/
/////////////////////////////////////////////////////////////////////////////////
/*
 * This process aggregates all expressed genes and splits this list by chromosome into the blocks and returns them as a file
 */
process Genes_expressed
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
	blocks      = genes_blocks
	minimum_exp = "$minimum_expression"
	template "Genes_expressed.pl"
}

/*
 * This process takes block of genes generated by "Genes_expressed" and generates ABC matrix for all tissues
*/
process Genes_data
{
    cache false

    input:
	val data_folder
	val samples_index
	tuple val(chrom), val(list), val(score)

    output:
	val "done"

    shell:
	data_folder   = "$data_folder"
	chrom         = "$chrom"
	list          = "$list"
	samples_index = "$samples_index"
	template "Genes_data.R"
}

/*
 * This process scans all "$tissue.regions.data" files one by one. For every tissue it obtains
 * all expressed genes, splits this list into blocks and returns this list for further analysis
*/
process Analysis_tissues
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
	blocks      = analysis_blocks
	chromosomes = "$chromosomes"
	minimum_exp = "$minimum_expression"
	template "Analysis_tissues.pl"
}

/*
 * This process runs Analysis
*/

process Analysis_run
{
    cache false

    input:
	val data_folder
	tuple val(task_index), val(tissue), val(chromosome), val(list)

    output:
	val "done"

    shell:
	data_folder  = "$data_folder"
	tissue       = "$tissue"
	task_index   = "$task_index"
	chromosome   = "$chromosome"
	list         = "$list"
	template "Analysis_run.R"
}

/////////////////////////////
/**** Making the report ****/
/////////////////////////////
/*
 * This process uses convetional ABC algorithm to make EPIraction_ABC predictions
 */
process EPIraction_ABC
{
    cache false

    input:
	val data_folder
	tuple val(index), val(tissue), val(short), val(RNA_seq), val(CTCF), val(Open), val(Cofactor), val(HiC), val(Open_min), val(H3K27ac_min), val(CTCF_min)
	val status

    output:
	val tissue

    shell:
	data_folder  = "$data_folder"
	tissue       = "$tissue"
	minimum_exp  = "$minimum_expression"
	template "EPIraction_ABC.pl"
}

/*
 * This process agjust ABC scores for cooccurrence and report
 */
process Report_pairs
{
    cache false

    input:
	val data_folder
	tuple val(index), val(tissue), val(short), val(RNA_seq), val(CTCF), val(Open), val(Cofactor), val(HiC), val(Open_min), val(H3K27ac_min), val(CTCF_min)
	val chromosomes
	val status

    output:
	val tissue

    shell:
	data_folder  = "$data_folder"
	tissue       = "$tissue"
	genome_file  = "$projectDir/files/Human.genome"
	interact     = "$projectDir/files/interact.as"
	HiC_baseline = "$projectDir/files/HiC.baseline"
	minimum_exp  = "$minimum_expression"
	chromosomes  = "$chromosomes"
	version      = "$version"

	template "Report.pl"
}

workflow
{
    indexData   = channel.fromPath(data_folder+"/files/"+tissues_index) | splitCsv(header:true, sep: '\t') | toList
    All_tissues = indexData.flatMap{ row -> row.tissue }                                                   | toList
    HiC_samples = indexData.flatMap{ row -> row.HiC    } | flatMap{  it -> it.split(";") } | unique        | toList

//################################################################################################
//### Here I develop tissue-consensus Activities based on Open chromatin, Cofactor and H3K27ac ###
//################################################################################################
    activity_done = Regions_annotate(data_folder,regions,indexData.flatMap()) | collect

//#############################################################################
//### Here I write Activities for each enhancer in each individual H3K27ac  ###
//#############################################################################
    samples_done  = Samples_activity(data_folder,samples_index,All_tissues.flatMap().combine(chrom_batches.flatMap()),activity_done) | collect
    samples_merge = Samples_merge(data_folder,chromosomes.flatMap(),All_tissues,samples_done)                                        | collect

//####################################################################
//### Here I write all gene-enhancer pairs, dump Hi-C data and QC  ###
//####################################################################
    pairs_done    = Make_pairs(data_folder,regions,chromosomes.flatMap()) | collect

    hic_consensus = HiC_consensus_dump(data_folder,regions,chromosomes.flatMap())                      | collect
    hic_done      = HiC_dump(data_folder,regions,HiC_samples.flatMap().combine(chromosomes.flatMap())) | collect

    hic_consensus = HiC_consensus_contacts(data_folder,regions,chromosomes,hic_consensus)              | collect
    hic_done      = HiC_contacts(data_folder,regions,HiC_samples.flatMap(),chromosomes,hic_done)       | collect

    first_qc      = First_qc(data_folder, All_tissues, HiC_samples, chromosomes, activity_done, samples_merge, pairs_done, hic_done, hic_consensus) | collect

//#########################################################################
//### Here I form tissue-specific Contact and Activity-by-Contact data. ###
//#########################################################################
    contacts_done  = Contacts_tissues(data_folder,indexData.flatMap(),chromosomes,first_qc)              | collect
    contacts_merge = Contacts_merge(data_folder,All_tissues,regions,chromosomes.flatMap(),contacts_done) | collect

//###############################################################
//### Here I form ABC matrixes for each gene and run analysis ###
//###############################################################
    All_genes   = Genes_expressed(data_folder,All_tissues,contacts_merge) | splitCsv(header:true, sep: '\t') | toList
    data_status = Genes_data(data_folder, samples_index, All_genes.flatMap())                                | collect

    Tis_genes     = Analysis_tissues(data_folder, All_tissues, chromosomes, data_status) | splitCsv(header:true, sep: '\t') | toList
    analysis_done = Analysis_run(data_folder,Tis_genes.flatMap())                                                           | collect

//##############################
//### Here I report the data ###
//##############################
    abc_done    = EPIraction_ABC(data_folder,indexData.flatMap(),analysis_done)           | collect
    report_done = Report_pairs(data_folder,indexData.flatMap(),chromosomes,analysis_done) | collect
}
