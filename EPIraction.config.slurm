
params {
    data_folder   = '/users/project/encode_005982/EPIraction.Nextflow.2025'
    tissues_index = 'EPIraction.tissues.data'
    samples_index = 'EPIraction.H3K27ac.data'
    regions       = 'EPIraction.regions.bed'
}

executor {
    name = 'slurm'
    queueSize = 100
    pollInterval = '10 sec'
    submitRateLimit = '10 sec'
}

process {
    beforeScript = 'module load R/4.3.3-gfbf-2023b Java/21.0.2 GSL/2.7-GCC-13.2.0'

///////////
    withName: Regions_annotate {
	cpus     = 2
	memory   = '6.GB'
	time     = '1h'
	array    = 100
    }

///////////
    withName: Samples_activity {
	cpus     = 2
	memory   = '6.GB'
	time     = '1h'
	array    = 500
    }

    withName: Samples_merge {
	cpus     = 2
	memory   = '6.GB'
	time     = '1h'
	array    = 25
    }

///////////
    withName: Make_pairs {
	cpus     = 2
	memory   = '2.GB'
	time     = '1h'
	array    = 25
    }

///////////
    withName: HiC_consensus_dump {
	cpus     = 2
	memory   = '10.GB'
	time     = '1h'
	array    = 25
    }

    withName: HiC_dump {
	cpus     = 2
	memory   = '10.GB'
	time     = '1h'
	array    = 2000
    }

    withName: HiC_consensus_contacts {
	cpus     = 2
	memory   = '4.GB'
	time     = '1h'
    }

    withName: HiC_contacts {
	cpus     = 2
	memory   = '4.GB'
	time     = '1h'
	array    = 100
    }


///////////
    withName: First_qc {
	cpus     = 1
	memory   = '2.GB'
	time     = '1h'
    }

///////////
    withName: Contacts_tissues {
	cpus     = 2
	memory   = '4.GB'
	time     = '2h'
	array    = 100
    }
    withName: Contacts_merge {
	cpus     = 3
	memory   = '20.GB'
	time     = '2h'
	array    = 25
    }

///////////
    withName: Genes_expressed {
	cpus     = 2
	memory   = '5.GB'
	time     = '2h'
    }

    withName: Genes_data {
	cpus     = 2
	memory   = '15.GB'
	time     = '1h'
	array    = 4000
    }

///////////
    withName: Analysis_tissues {
	cpus     = 1
	memory   = '2.GB'
	time     = '1h'
    }
    withName: Analysis_run {
	cpus     = 1
	memory   = '2.GB'
	time     = '1h'
	array    = 5000
	errorStrategy = { task.exitStatus == 143 ? 'retry' : 'terminate' }
	maxRetries = 3
    }

///////////
    withName: EPIraction_ABC {
	cpus     = 2
	memory   = '8.GB'
	time     = '1h'
	array    = 100
    }

    withName: Report_pairs {
	cpus     = 2
	memory   = '8.GB'
	time     = '1h'
	array    = 100
    }

    clusterOptions = { task.time <= 1.h ? '--qos=vshort' : 
	(task.time <= 3.h ? '--qos=shorter' : 
        (task.time <= 6.h ? '--qos=short' : 
        (task.time <= 12.h ? '--qos=normal' : 
        (task.time <= 24.h ? '--qos=long' : 
        (task.time <= 48.h ? '--qos=vlong' : '--qos=marathon' ))))) }
}
