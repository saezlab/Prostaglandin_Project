# PHONEMeS-Inputs

 Here we build the input objects from the data for the PHONEMeS analysis. The input objects are generated simply by running the R scripts `buildObject_CD4.R` and `buildObject_CD8.R` for *CD4* and *CD8* respectively.
 
 ## Steps to build the inputs
 
 Below we describe the main steps followed to build the PHONEMeS inputs.

+ We read the Mass-Spect data `CD4_New.xlsx` and `CD8_New.xlsx` containing the site level information (fold changes and significance).
+ We generate a table mapping the quantified phospho-proteins from uniprot to gene identifiers.
+ We assign a threshold value `pValThresh = 0.05` to the differential abundance.
+ Based on the data tables (`CD4_New.xlsx` and `CD8_New.xlsx`), we generate the GMM objects containing a list of matrices for each site. Each matrix contains different values and labels assigned to each measurement at all time-points as below:
+ The `Indiv` columns contain the scores we assign to each of the measurements by computing the log2 value of the ratio between the significance of a measurement at a specific time-point (as found in the data tables) over the `pValThresh` value. Sites with a `pValue < pValThresh` on a specific time-point will be assigned a negative score, otherwise they will be assigned a positive score.
+ The `clus` columns contain information about significance status of the measurement. If `pValue < pValThresh`, then to the data-point it will be assigned a Perturbation `P` status, otherwise the measurement will take a Control `C` status.
+ The `FCvCaPval` contain the `pValue` significance values for each measurement.
+ The `status` columns show whether a measurement has been observed to be higly regulated (`OK`) or not (`FP`). We do not apply fold change (FC) thrreshods for the PHONEMeS analysis.
+ The `FCvC` column contains the log2 fold changes of activated sites compared to control for each experimetnal condition.

## Generating and storing the inputs
We build two separate data objects `dataGMM_CD4.RData` and `dataGMM_CD8.RData` and we keep them stored on the local directory. These inputs objects are generated simply by executing the `buildObject_CD4.R` and `buildObject_CD8.R` scripts.
