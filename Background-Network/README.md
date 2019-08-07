# Background network

PHONEMeS is a method that trains a prior knowledge network to phosphorylation data. To adapt this method to GCPR signalling the background network (representing our prior knowledge), was build as follows: First, the background kinase-to-substrate (K-S) network obtained from Omnipath. Then, this was complemented with a subset of all directed and signed protein-protein interactions (PPIs) obtained also from Omnipath which are associated with GPCR signalling as defined by Reactome (Pathway R-HSA-372790). Closely related heterotrimeric G-Protein subunits and kinase isoforms were grouped in the prior knowledge network.

## Generating and storing the background network

Our background network is generated simply by executing the `buildBN.R` script and which applies all the steps described below. The output of the analysis is the `allD.RData` object which contains the list of all 26367 interactions used for the training of the models for both the cell lines.

## Steps to build the background network

Below we describe the main resources used to build our background network.

+ We start by adding all the Omnipath K-S interactions as loaded by the `import_Omnipath_PTMS()` finction of **OmnipathR** package.
+ We identify the protein associated with GPCR signalling in Reactome (Pathway R-HSA-372790). This list of proteins is the one stored in inst/geneset.txt file.
+ From the list of protein interactions of Omnipath (as loaded by the `import_Omnipath_Interactions()` function of **OmnipathR** package), we select only those that are signed and directed and which are involving the proteins associated to GPCR signalling.
+ We add interactions from canonical knowledge of the Prostaglandin signalling system as described in the lines `[65-16]` of the `buildBN.R` script.
+ We group the isoforms of several G-Proteins and other protein kinases as described in the lines `[164-339]` based on their functional similarities.


## Side Note

+ We have used the Omnipath version of August 2019 to build our background network.
