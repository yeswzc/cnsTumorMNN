# cnsTumorMNN
* Query nearest neighbors in 'Capper et al.. Nature. 2018' for central nervous tumors. Meanwhile, umap, tsne will also be plotted.

* R/3.6 is suggested.

* (The package is not been well tested yet.)

# How to install
* There are many packages required before installations. The packages include but not limitted to umap, Rtsne, minfi, LUMP, RFpurify, dplyr, caret, glmnet, plotly, RANN, knitr, kableExtra, RColorBrewer...
### Install LUMP 
devtools::install_github("yeswzc/LUMP")
### Install cnsTumorMNN, an imperfect package
install.packages("cnsTumorMNN_0.1.0.tar.gz", repos = NULL, type = "source")

#devtools::install_github("yeswzc/cnsTumorMNN") #This may not work.

# How to Run the code
After install finished.
Download 'pipeline.v5.R' and 'Generate_HTMLreport.v5.Rmd', then you can run <br> </br>
`Rscript /path/to/pipeline.v5.R /path/to/idat_directory/ output_prefix`


# Reference: 
Wu et al.. Impact of the methylation classifier and ancillary methods on CNS tumor diagnostics. Neuro-oncology, 2021 (https://doi.org/10.1093/neuonc/noab227).
