##
# Function to group the nodes
groupNodes = function(sif = sif, nodesAttributes = nodesAttributes, sifName = "sif_grouped.txt", attribName = "attrib_grouped.txt"){
  
  idx = c()
  for(ii in 1:nrow(nodesAttributes)){
    
    if(is.na(nodesAttributes[ii, 2])){
      
      if(length(strsplit(x = as.character(nodesAttributes[ii, 1]), split = "_")[[1]]) > 1){
        
        idx = c(idx, ii)
        
      }
      
    }
    
  }
  
  species = nodesAttributes$Species[idx]
  speciesProt = c()
  residues = c()
  mapping = matrix(data = , nrow = length(species), ncol = 2)
  mapping[, 1] = species
  for(ii in 1:length(species)){
    speciesProt = c(speciesProt, strsplit(x = species[ii], split = "_")[[1]][1])
    residues = c(residues, strsplit(x = species[ii], split = "_")[[1]][2])
  }
  uSpeciesProt = unique(speciesProt)
  sites = matrix(data = , nrow = length(uSpeciesProt), ncol = 2)
  sites[, 1] = uSpeciesProt
  for(ii in 1:nrow(sites)){
    
    idx = which(speciesProt==sites[ii, 1])
    sites[ii, 2] = paste0(sites[ii, 1], "_")
    for(jj in 1:length(idx)){
      sites[ii, 2] = paste0(sites[ii, 2], residues[idx[jj]], ";")
    }
    sites[ii, 2] = substr(sites[ii, 2], 1, nchar(sites[ii, 2])-1)
    mapping[idx, 2] = sites[ii, 2]
    
  }
  
  for(ii in 1:nrow(nodesAttributes)){
    
    if(as.character(nodesAttributes[ii, 1])%in%mapping[, 1]){
      
      idx = which(mapping[, 1]==as.character(nodesAttributes[ii, 1]))
      nodesAttributes[ii, 1] = mapping[idx, 2]
      
    }
    
  }
  nodesAttributes = unique(nodesAttributes)
  write.table(x = nodesAttributes, file = attribName, quote = FALSE, sep = "\t", row.names = FALSE)
  
  #
  for(ii in 1:nrow(sif)){
    
    idx = which(mapping[, 1]==as.character(sif[ii, 1]))
    if(length(idx)>0){
      
      sif[ii, 1] = mapping[idx, 2]
      
    }
    
    idx = which(mapping[, 1]==as.character(sif[ii, 3]))
    if(length(idx)>0){
      
      sif[ii, 3] = mapping[idx, 2]
      
    }
    
  }
  
  for(ii in 1:nrow(sif)){
    
    ss = as.character(sif[ii, 1])
    tt = as.character(sif[ii, 3])
    
    idx1 = which(sif$Source==ss)
    idx2 = which(sif$Target==tt)
    idx = intersect(x = idx1, y = idx2)
    sif$Weight[idx] = as.character(mean(as.numeric(sif$Weight[idx])))
    
  }
  
  sif = unique(sif)
  
  idx2rem = c()
  for(ii in 1:nrow(sif)){
    
    if(length(strsplit(x = as.character(sif[ii, 1]), split = "_")[[1]])>1){
      ss1 = strsplit(x = as.character(sif[ii, 1]), split = "_")[[1]][1]
      ss2 = strsplit(x = as.character(sif[ii, 1]), split = "_")[[1]][2]
      if(ss1==ss2){idx2rem=c(idx2rem, ii)}
    }
    
    if(length(strsplit(x = as.character(sif[ii, 3]), split = "_")[[1]])>1){
      ss1 = strsplit(x = as.character(sif[ii, 3]), split = "_")[[1]][1]
      ss2 = strsplit(x = as.character(sif[ii, 3]), split = "_")[[1]][2]
      if(ss1==ss2){idx2rem=c(idx2rem, ii)}
    }
    
  }
  
  write.table(x = sif, file = sifName, quote = FALSE, sep = "\t", row.names = FALSE)
  
}

##
# grouping nodes and saving resulting networks
library(readr)
for(ii in 1:5){
  
  sif1 = read_delim(paste0("CD4_cplex_exp", ii, ".txt"), "\t", escape_double = FALSE, trim_ws = TRUE)
  nodesAttributes = read_delim(paste0("nodesAttributes_", ii, ".txt"), "\t", escape_double = FALSE, trim_ws = TRUE)
  file.remove(paste0("CD4_cplex_exp", ii, ".txt"))
  file.remove(paste0("nodesAttributes_", ii, ".txt"))
  
  groupNodes(sif = sif1, nodesAttributes = nodesAttributes, sifName = paste0("../../Results/PHONEMeS/CD4/cplex_grouped_", ii, ".txt"), 
             attribName = paste0("../../Results/PHONEMeS/CD4/attributes_grouped_cplex_", ii, ".txt"))
  
}

##
# removing redundant remaining files
file.remove("clone0.log")
file.remove("clone1.log")
file.remove("clone2.log")
file.remove("clone3.log")
file.remove("cplex")
file.remove("cplex_exp_combined.txt")
file.remove("nodesAttributes_combined.txt")
file.remove("nodesAttributes.txt")
