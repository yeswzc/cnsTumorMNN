# cnsTumorMNN
Query nearest neighbors in 'Capper et al.. Nature. 2018' for central nervous tumors

# How to install
### Install LUMP
devtools::install_github("yeswzc/LUMP")
### Install minfi
BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("IlluminaHumanMethylationEPICmanifest")
install.packages("https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-14.tar.gz", repos = NULL, type = "source")

devtools::install_github('mwsill/RFpurify')

BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("minfi")

