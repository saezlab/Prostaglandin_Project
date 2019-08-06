# Prostaglandin E2 Signaling Pathways in T cell Subsets

Prostaglandin E2 signaling networks in T cells revealed through a systems approach.

## Purpose

Here we provide documentation and the analysis scripts used for generating the PHONEMeS network models of Prostaglandin signalling pathways in CD4 and CD8 T-cell subsets.

## Organisation of the repository

+ [Data](https://github.com/saezlab/Prostaglandin_Project/tree/master/Data): Here we provide the large-scale phospho-proteomics dataset for the CD4 and CD8 cells.

+ [Background-Network](https://github.com/saezlab/Prostaglandin_Project/tree/master/Background-Network): Here we prepare the background network which is used for contextualizing the CD4 & CD8 signalling pathways through PHONEMeS.

+ [PHONEMeS-Inputs](https://github.com/saezlab/Prostaglandin_Project/tree/master/PHONEMeS-Inputs): Here we prepare the PHONEMeS inputs objects from the data.

+ [PHONEMeS-Analysis](https://github.com/saezlab/Prostaglandin_Project/tree/master/PHONEMeS-Analysis): Here we perform the PHONEMeS analysis to train our models.

+ [Results](https://github.com/saezlab/Prostaglandin_Project/tree/master/Results/): Here we store the analysis results.

## Requirements

All network modelling steps were performed in [R](https://www.rstudio.com/) v3.5.1 and visualised using [Cytoscape](https://cytoscape.org/) v3.3.

For running the PHONEMeS analysis, the user must obtain the cplex [license](https://www.ibm.com/products/ilog-cplex-optimization-studio?S_PKG=CoG&cm_mmc=Search_Google-_-Data+Science_Data+Science-_-WW_IDA-_-+IBM++CPLEX_Broad_CoG&cm_mmca1=000000RE&cm_mmca2=10000668&cm_mmca7=9041989&cm_mmca8=kwd-412296208719&cm_mmca9=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_&cm_mmca10=267798126431&cm_mmca11=b&mkwid=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_%7C470%7C135655&cvosrc=ppc.google.%2Bibm%20%2Bcplex&cvo_campaign=000000RE&cvo_crid=267798126431&Matchtype=b&gclid=Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB)
from IBM, which is for free for academic purposes. Once the user obtains the license, they should store the cplex executable file to the 
working directory where PHONEMeS is running (here: on the [Analysis](https://github.com/saezlab/EDN_phospho/tree/master/Analysis) directory).

Other prerequisites include downloading and installing the following `R` package dependencies:

+ [PHONEMeS](https://saezlab.github.io/PHONEMeS/)
+ [igraph](https://igraph.org/r/)
+ [BioNet](https://bioconductor.org/packages/release/bioc/html/BioNet.html)
+ [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
+ [XML](https://cran.r-project.org/web/packages/XML/index.html)

The Analysis was run under MAC OS X El Captain machine. Please note that the results might slightly change when running the analysis from machines with different OS.

## License

Distributed under the [GNU GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.html).
