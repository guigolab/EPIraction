// Use DSL2
nextflow.enable.dsl=2

data_folder   = params.data_folder
tissues_index = params.tissues_index
samples_index = params.samples_index
regions       = params.regions
chromosomes   = channel.of("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX").toList()
t_status = "yes"
c_status = "yes"
q_status = "yes"
n_status = "yes"
m_status = "yes"
a_status = "yes"
p_status = "yes"
h_status = "yes"

s_status = "yes"
d_status = "yes"
g_status = "yes"
first_qc = '/dev/null'

// Print log information
log.info """
Welcome to EPIraction 1.2, nextflow implementation

data_folder : $data_folder
tissues     : $tissues_index
samples     : $samples_index
regions     : $regions

"""

/*
 * This process takes tissue data from tissues_index file and annotate regions
 * with molecular information from corresponding BigWig files and RNA-seq data
 */
process Annotate_regions
{
    cache false

    input:
	val data_folder
	val regions
	tuple val(index), val(tissue), val(RNA_seq), val(Open), val(Cofactor), val(HiC)

    output:
	val data_folder
	val tissue

    shell:
	data_folder = "$data_folder"
	regions     = "$regions"
	tissue      = "$tissue"
	RNA_seq     = "$RNA_seq"
	Open        = "$Open"
	Cofactor    = "$Cofactor"
	template "Annotate.pl"
}

/*
 * This process calculates activities for enhancers
 */
process Activity_regions
{
    cache false

    input:
	val data_folder
	val tissue

    output:
	val tissue

    shell:
	data_folder = "$data_folder"
	tissue      = "$tissue"
	template "Activity.pl"
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


/*
 * This process makes H3K27ac matrix for every
 * tissue and chromosome in temp folder
 */
process H3K27ac
{
    cache false

    input:
	val data_folder
	val samples
	val regions
	val chromosomes
	val tissue

    output:
	val tissue

    shell:
	data_folder   = "$data_folder"
	samples_index = "$samples"
	regions       = "$regions"
	tissue        = "$tissue"
	chromosomes   = "$chromosomes"
	template "H3K27ac.pl"
}

/*
 * This process makes H3K27ac matrix for tissue 
 * H3K27ac consensus BigWigs in temp folder
 */
process Consensus
{
    cache false

    input:
	val data_folder
	val regions
	val tissues
	val chromosome

    output:
	val chromosome

    shell:
	data_folder = "$data_folder"
	regions     = "$regions"
	tissues     = "$tissues"
	chromosome  = "$chromosome"
	template "Consensus.pl"
}

/*
 * This process runs quantile normalization on consensus tissue H3K27ac mextrix
 */
process Quantiles_consensus {
    cache false

    input:
	val data_folder
	val regions
	val chromosomes
	val status

    output:
	val "done"

    shell:
	data_folder = "$data_folder"
	regions     = "$regions"
	chromosomes = "$chromosomes"
	template "Quantiles.R"
}

/*
 * This process calculates rescalings for every region in every tissue
 * based on quantile normalization results
 */
process Tissues_rescale {
    cache false

    input:
	val data_folder
	val chromosome
	val status

    output:
	val "done"

    shell:
	data_folder = "$data_folder"
	chromosome  = "$chromosome"
	template "Tissues_rescale.pl"
}


/*
 * This process normalize H3K27ac data within each tissue
 */
process Normalize_H3K27ac {
    cache false

    input:
	val data_folder
	val tissue
	val chromosomes
	val status_s
	val status_n

    output:
	val tissue

    shell:
	data_folder = "$data_folder"
	chromosomes = "$chromosomes"
	tissue      = "$tissue"
	template "Normalize_H3K27ac.R"
}

/*
 * This process merges all normalized H3K27ac data by chromosome
 */
process Merge_H3K27ac {
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
	template "Merge_H3K27ac.R"
}


/*
 * This process dumps HiC data. For all promoters versus all enhancers within 1M.
 */
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
	data_folder = "$data_folder"
	sample      = "$sample"
	chrom       = "$chrom"
	regions     = "$regions"
	v_return    = "$sample.$chrom"

	template "HiC_dump.pl"
}

/*
 * This is first QC. Need redion annotation for each tissue, all promoter-enhancer
 * pairs, dump and merge H3K27ac data and HiC dump
 */
process First_qc
{
    cache false

    input:
	val data_folder
	val tissues
	val hic_samples
	val chromosomes
	val t_status
	val c_status
	val q_status
	val n_status
	val m_status
	val a_status
	val p_status
	val h_status

    output:
	path "First_QC.log"

    shell:
	data_folder = "$data_folder"
	tissues     = "$tissues"
	hic_samples = "$hic_samples"
	chromosomes = "$chromosomes"
	template "First_qc.pl"
}


/*
 * This process annotates all valid promoter-enhancer pairs with HiC data for every tissue. 
 * It selects the best contact probability from each associated HiC sample or baseline.
 * It annotates all valid promoter-enhancer with scalings factors.
 */
process Tissue_scalings
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
	data_folder = "$data_folder"
	tissue      = "$tissue"
	regions     = "$regions"
	HiC_samples = "$HiC"
	chromosomes = "$chromosomes"
	template "Tissue_scalings.pl"
}

/*
 * This process aggregates tissue-spcific scaling data into matrices for each chromosome
 */
process Scalings_merge
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
	template "Scalings_merge.pl"
}

/*
 * This process scans "$tissue.regions.data" files, aggregates all expressed genes,
 * splits this list by chromosome into the blocks and returns them as a file
 */
process Expressed_genes
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
	blocks      = 50
	template "Expressed_genes.pl"
}

/*
 * This process takes each block of genes generated by "Expressed_genes". For every gene it generates
 * enhancer matrix taking into account scalign factors, simplifies this matrix to avoid small values,
 * adds 5% promoter, performs PCA and keeps ony corelated components. Two files are generated for each gene:
 * "$gene_id.enhancers.matrix.gz" contains the matrix itself, "$gene_id.rotations.matrix.gz" contains rotations.
*/
process Glmnet_data
{
    cache false

    input:
	val data_folder
	val samples_index
	tuple val(group), val(list)

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
process Tissue_genes
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
	blocks      = 300
	chromosomes = "$chromosomes"
	template "Tissue_genes.pl"
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
	data_folder = "$data_folder"
	tissue      = "$tissue"
	task_index  = "$task_index"
	chromosome  = "$chromosome"
	list        = "$list"
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
	data_folder = "$data_folder"
	tissue      = "$tissue"
	template "Report.pl"
}

workflow
{
    indexData   = channel.fromPath(data_folder+"/files/"+tissues_index) | splitCsv(header:true, sep: '\t') | toList
    All_tissues = indexData.flatMap{ row -> row.tissue }                                                   | toList
    HiC_samples = indexData.flatMap{ row -> row.HiC    } | flatMap{  it -> it.split(";") } | unique        | toList

//    a_status = Annotate_regions(data_folder,regions,indexData.flatMap()) | Activity_regions       | collect
//    t_status = H3K27ac(data_folder,samples_index,regions,chromosomes,All_tissues.flatMap())       | collect
//    c_status = Consensus(data_folder,regions,All_tissues,chromosomes.flatMap())                   | collect
//    q_status = Quantiles_consensus(data_folder,regions,chromosomes,c_status)                      | collect
//    n_status = Tissues_rescale(data_folder,chromosomes.flatMap(),q_status)                        | collect
//    t_status = Normalize_H3K27ac(data_folder,All_tissues.flatMap(),chromosomes,t_status,n_status) | collect
//    m_status = Merge_H3K27ac(data_folder,chromosomes.flatMap(),All_tissues,t_status)              | collect
//    p_status = Make_pairs(data_folder,regions,chromosomes.flatMap())                              | collect
//    h_status = HiC_dump(data_folder,regions,HiC_samples.flatMap().combine(chromosomes.flatMap())) | collect
//    first_qc = First_qc(data_folder, All_tissues, HiC_samples, chromosomes, t_status, c_status, q_status, n_status, m_status, a_status, p_status, h_status) | collect

//    s_status = Tissue_scalings(data_folder,indexData.flatMap(),chromosomes,first_qc)          | collect
//    m_status = Scalings_merge(data_folder,All_tissues,regions,chromosomes.flatMap(),s_status) | collect

//    All_genes = Expressed_genes(data_folder,All_tissues,m_status) | splitCsv(header:true, sep: '\t') | toList
//    d_status  = Glmnet_data(data_folder, samples_index, All_genes.flatMap())                         | collect

//    Tis_genes = Tissue_genes(data_folder, All_tissues, chromosomes, d_status) | splitCsv(header:true, sep: '\t') | toList
//    g_status  = Glmnet_run(data_folder,Tis_genes.flatMap())            | collect
    f_status  = Report_pairs(data_folder,indexData.flatMap(),g_status) | collect
}
