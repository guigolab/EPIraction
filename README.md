# Welcome to EPIraction 1.6, nextflow implementation.

Our pipeline needs a lot of computations. In order to run it efficiently and smoothly we wrap all the scripts in the nextflow environment to run them into the HPC cluster.
To run this pipeline you need to download 1,179 Gb of bigWig files and 443 Gb of Hi-C data. You can find links to these data below.
Please estimate your resource properly before downloading the data and attempting to run the pipeline.<br><br>

## Softwares:
Linux operating system with bash, perl, sort and other default executables.<br>
nextflow - version 24.10.3 or above. That supports DSL2 and arrays https://www.nextflow.io<br>
bedtools - any more or less recent version https://bedtools.readthedocs.io/en/latest<br>
bigWigAverageOverBed - https://github.com/ucscGenomeBrowser/kent<br>
bedToBigBed - https://github.com/ucscGenomeBrowser/kent<br>
R - starting from 4.0.0 https://cran.r-project.org <br>
pigz - any working version https://github.com/madler/pigz<br><br>

## R libraries:
R.utils 2.12.2 - https://cran.r-project.org/web/packages/R.utils<br>
data.table 1.14.6 - https://cran.r-project.org/web/packages/data.table<br>
Rfast 2.0.6 - https://cran.r-project.org/web/packages/Rfast<br>

## Obtaining EPIraction data:
You need to create a specific folder ("$work_folder") that would store all the input data, all the processed data and all the predictions.<br>
You need to download this file: https://public-docs.crg.es/rguigo/Data/rnurtdinov/EPIraction.data/EPIraction.tar.gz<br>
put it into $work_folder and extract the content:<br>```tar xzvf EPIraction.tar.gz```<br><br>
You need to download all bigWig files from: https://public-docs.crg.es/rguigo/Data/rnurtdinov/EPIraction.data/BigWig/ into "$work_folder/BigWig"<br>
You need to download all Hi-C files from: https://public-docs.crg.es/rguigo/Data/rnurtdinov/EPIraction.data/HiC.files/ into "$work_folder/HiC.files"<br>
You need to download all RSEM files from: https://public-docs.crg.es/rguigo/Data/rnurtdinov/EPIraction.data/RSEM/ into "$work_folder/RSEM"<br><br>
You need to run:<br>```$work_folder/QC_files.pl```<br><br>to test if all necessary files are present.<br><br>

## Downloading and running the pipeline:
We recommend to download the pipeline into a folder that is different from "$work_folder" above:<br>
```git clone https://github.com/guigolab/EPIraction.git```<br><br>
You need to edit the EPIraction.config.slurm file. Put the actual full path to "$work_folder" into<br>
```data_folder  = '*****'``` variable<br><br>
Edit the "executor" section:<br>
Set "name" value, consult https://www.nextflow.io/docs/latest/executor.html<br>
Other parameters are not so essential, consult https://www.nextflow.io/docs/latest/config.html#scope-executor<br><br>
Edit the "process" section:<br>
If you use modules like we or conda envirement please use "beforeScript" variable.
Please clarify if the job parameters "cpu", "memory" and "time" are compatible with your executor and have correct format.<br><br>

Run the pipeline:<br>
```nextflow -c EPIraction.config.slurm run EPIraction```<br>

The complete run of the pipeline takes several hours. Please be sure that your Linux administrators do not kill this nextflow process.<br>
Nextflow pipeline stores and runs its own stuff within the "work" folder. However, the main location of all temporary and final files is "$work_folder".<br>
The scripts are internally configured to search for proper intermediate files. Nextflow caching is disabled, do not rely on it.<br><br>
All the output files will be stored within the "$work_folder/report" folder.
