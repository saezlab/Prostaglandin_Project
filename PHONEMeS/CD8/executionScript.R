# Load packages
library(BioNet)
library(igraph)
library(PHONEMeS)
library(hash)
library(dplyr)

# Call the PHONEMeS functions
source("../Public/buildDataMatrix.R")
source("../Public/ilpFunctions.R")
source("../Public/buildDataObject.R")
source("../Public/build_Nw.R")
source("../Public/build_PKN.R")

load(file = "../buildBN/allD.RData")

load(file = "../build-Object/CD8/dataGMM_CD8.RData")

GMM.ID$S.cc <- GMM.ID$dataID
GMM.res.noFC <- GMM
GMM.res <- GMM.wFC
GMM.res.ID <- GMM.ID

#Make the data objects that will be needed
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=GMM.res.noFC, IDmap=GMM.res.ID, resFC=GMM.res)

conditions <- list()
for(i in 1:nrow(GMM[[1]])){
  
  conditions[[length(conditions)+1]] <- rownames(GMM[[1]])[i]
  
}

names(conditions) <- rownames(GMM[[i]])

targets.P<-list(cond1=c("PTGER1"), cond2=c("PTGER2"), cond3=c("PTGER3"), cond4=c("PTGER4"), 
                cond5=c("PTGER1", "PTGER2", "PTGER3", "PTGER4"))

targets <- targets.P

sifAll <- matrix(data = , nrow = 1, ncol = 3)
colnames(sifAll) <- c("Source", "f50", "Target")

for(ii in 1:6){
  
  if(ii < 6){
    
    experiments=ii
    
    source("../Public/runPHONEMeS.R")
    resultsMulti <- runPHONEMeS(targets.P = targets.P, conditions = conditions, dataGMM = dataGMM, experiments = experiments, bg = bg, solver = "cplex", nSolutions = 100, nK = "all", timelimit = 600, case = ii)
    resultsMulti <- removeRedundantNodes(resultsSIF1 = resultsMulti)
    write.table(x = resultsMulti, file = paste0("CD8_cplex_exp", ii, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    sifAll <- unique(rbind(sifAll, as.matrix(resultsMulti)))
    
    #Assigning attributes for better visualization
    source("../Public/assignAttributes.R")
    nodesAttributes <- assignAttributes(sif = sifAll[-1, ], dataGMM = dataGMM, targets = targets.P)
    write.table(x = nodesAttributes, file = paste0("nodesAttributes_", ii, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
  } else {
    
    experiments=1:5
    
    source("../Public/runPHONEMeS.R")
    resultsMulti <- runPHONEMeS(targets.P = targets.P, conditions = conditions, dataGMM = dataGMM, experiments = experiments, bg = bg, solver = "cplex", nSolutions = 100, nK = "all", timelimit = 600, case = 5)
    resultsMulti <- removeRedundantNodes(resultsSIF1 = resultsMulti)
    write.table(x = resultsMulti, file = paste0("CD8_cplex_combined.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    sifAll <- unique(rbind(sifAll, as.matrix(resultsMulti)))
    
    #Assigning attributes for better visualization
    source("../Public/assignAttributes.R")
    nodesAttributes <- assignAttributes(sif = sifAll[-1, ], dataGMM = dataGMM, targets = targets.P)
    write.table(x = nodesAttributes, file = paste0("nodesAttributes_combined.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
  }
  
}