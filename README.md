# cnsTumorMNN
Query nearest neighbors in 'Capper et al.. Nature. 2018' for central nervous tumors
R/3.6 is suggested.

# How to install
#### There are many packages required before installations. The packages include but not limitted to umap, Rtsne, minfi, LUMP, RFpurify, dplyr, caret, glmnet, plotly, RANN, knitr, kableExtra, RColorBrewer...
### Install LUMP (not been well tested yet)
devtools::install_github("yeswzc/LUMP")
### Install cnsTumorMNN, an imperfect package
install.packages("cnsTumorMNN_0.1.0.tar.gz", repos = NULL, type = "source")

devtools::install_github("yeswzc/cnsTumorMNN")

# How to Run the code
After install finished.
Download 'pipeline.v5.R' and 'Generate_HTMLreport.v5.Rmd', then you can run
Rscript /path/to/pipeline.v5.R /path/to/idat_directory/ output_prefix
