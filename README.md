Welcome to EPIraction 1.2, nextflow implementation.

At the moment we are distributing the complete functional pipeline with small demo data. They contain four tissues: Blood.Lymphoid.GM12878, Brain.Cortex.Primary, Mesenchyme.Stem and Muscle.Smooth.Artery. 
We will release the complete data for 83 tissues after our publication is accepted.
However all the predicted enhancer-gene links are available at Encode portal

Our pipeline needs a lot of computations. In order to run it efficiently and smoothly we wrap all the scripts in nextflow. We recommend creating a specific conda environment (EPIraction) for this project and install all software within.

Software:
nextflow - any version that supports DSL2 suits https://www.nextflow.io
bedtools - any more or less recent version https://bedtools.readthedocs.io/en/latest
bigWigAverageOverBed - https://github.com/ucscGenomeBrowser/kent
R - starting from 4.0.0 https://cran.r-project.org 

R libraries:
R.utils 2.12.2 - https://cran.r-project.org/web/packages/R.utils
data.table 1.14.6 - https://cran.r-project.org/web/packages/data.table
glmnet 4.1-6 - https://cran.r-project.org/web/packages/glmnet
Rfast 2.0.6 - https://cran.r-project.org/web/packages/Rfast
nsprcomp 0.5.1-2 - https://cran.r-project.org/web/packages/nsprcomp
r-wcorr 1.9.5 - https://cran.r-project.org/web/packages/wCorr

Bioconductor R libraries (install them directly from R):
BiocManager::install("preprocessCore", configure.args="--disable-threading")

The easy way to create conda environment:
conda create -n EPIraction python=3.7
conda activate EPIraction
conda install -c bioconda ucsc-bigwigaverageoverbed nextflow bedtools
conda install -c conda-forge r-base r-r.utils r-data.table r-glmnet r-rfast r-biocmanager r-nsprcomp r-wcorr

#Use this command to reset openssl:
conda install -c conda-forge --force-reinstall openssl=1.1.1

#Bioconductor R libraries (install it directly from R):
BiocManager::install("preprocessCore", configure.args="--disable-threading")
