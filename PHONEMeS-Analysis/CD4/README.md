# PHONEMeS Analysis

This directory contains scripts necessary to perform the PHONEMeS analysis. PHONEMeS has been reformulated as an Integer Linear Programming (ILP) problem and this new version is the one used for this study. **For performing th ILP analysis, users must first obtain a [CPLEX license](https://www.ibm.com/products/ilog-cplex-optimization-studio?S_PKG=CoG&cm_mmc=Search_Google-_-Data+Science_Data+Science-_-WW_IDA-_-+IBM++CPLEX_Broad_CoG&cm_mmca1=000000RE&cm_mmca2=10000668&cm_mmca7=9041989&cm_mmca8=kwd-412296208719&cm_mmca9=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_&cm_mmca10=267798126431&cm_mmca11=b&mkwid=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_%7C470%7C135655&cvosrc=ppc.google.%2Bibm%20%2Bcplex&cvo_campaign=000000RE&cvo_crid=267798126431&Matchtype=b&gclid=Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB) (free for academic purposes) and store on this working directory the cplex executable file.**

## Steps to run the analysis

Separate analysis are run for the two different cell types: [*CD4*](https://github.com/saezlab/Prostaglandin_Project/tree/master/PHONEMeS-Analysis/CD4) and [*CD8*](https://github.com/saezlab/Prostaglandin_Project/tree/master/PHONEMeS-Analysis/CD8). On [*/CD4*](https://github.com/saezlab/Prostaglandin_Project/tree/master/PHONEMeS-Analysis/CD4) and [*/CD8*](https://github.com/saezlab/Prostaglandin_Project/tree/master/PHONEMeS-Analysis/CD8) are expalined the details of the analysis.

## Public directory

The [*/Public*](https://github.com/saezlab/Prostaglandin_Project/tree/master/PHONEMeS-Analysis/Public) directory contains all the finctions where the ILP problem is formulated and which are needed to run the analysis.
