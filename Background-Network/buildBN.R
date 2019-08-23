library(readr)
library(OmnipathR)

ptms = import_Omnipath_PTMS()
interactions = import_Omnipath_Interactions()

write.table(x = interactions, file = "ppi.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

geneset <- read_csv("geneset.txt")
geneset <- geneset[-1, ]

proteins = geneset$REACTOME_GPCR_DOWNSTREAM_SIGNALING

idx1 <- which(interactions$source_genesymbol%in%proteins)
idx2 <- which(interactions$target_genesymbol%in%proteins)

idx <- unique(c(idx1, idx2))

interactions <- interactions[idx, ]

idx = which(interactions$is_directed==1)

interactions = interactions[idx, ]

ss = interactions$is_stimulation + interactions$is_inhibition
idx = which(ss==1)

interactions = interactions[idx, ]

##
allD = matrix(data = , nrow = nrow(ptms), ncol = 8)
colnames(allD) = c("S.AC", "S.ID", "K.AC", "K.ID", "res", "pos", "SID", "S.cc")

allD[, 1] <- ptms$substrate_genesymbol
allD[, 2] <- ptms$substrate_genesymbol
allD[, 3] <- ptms$enzyme_genesymbol
allD[, 4] <- ptms$enzyme_genesymbol
allD[, 5] <- ptms$residue_type
allD[, 6] <- ptms$residue_offset
allD[, 8] <- paste0(ptms$substrate_genesymbol, "_", ptms$residue_type, ptms$residue_offset)

ppi = interactions

temp <- matrix(data = , nrow = 1, ncol = 8)
colnames(temp) <- colnames(allD)
for(ii in 1:nrow(ppi)){
  
  toBind <- matrix(data = , nrow = 1, ncol = 8)
  colnames(toBind) <- colnames(allD)
  
  toBind[1, 1] <- ppi$target_genesymbol[ii]
  toBind[1, 2] <- ppi$target_genesymbol[ii]
  toBind[1, 3] <- ppi$source_genesymbol[ii]
  toBind[1, 4] <- ppi$source_genesymbol[ii]
  toBind[1, 5] <- "R"
  toBind[1, 6] <- "1"
  toBind[1, 8] <- paste0(toBind[1, 2], "_R1")
  
  temp <- unique(rbind(temp, toBind))
  
  allD <- unique(rbind(allD, toBind))
  
}

##
# PTGER4 <- GNAS

toBind <- matrix(data = , nrow = 1, ncol = 8)
colnames(toBind) <- colnames(allD)

toBind[1, 1] <- "GNAS"
toBind[1, 2] <- "GNAS"
toBind[1, 3] <- "PTGER4"
toBind[1, 4] <- "PTGER4"
toBind[1, 5] <- "R"
toBind[1, 6] <- "1"
toBind[1, 8] <- paste0(toBind[1, 2], "_R1")

allD <- unique(rbind(allD, toBind))

##
# GNAS <- PRKACA

toBind <- matrix(data = , nrow = 1, ncol = 8)
colnames(toBind) <- colnames(allD)

toBind[1, 1] <- "PRKACA"
toBind[1, 2] <- "PRKACA"
toBind[1, 3] <- "GNAS"
toBind[1, 4] <- "GNAS"
toBind[1, 5] <- "R"
toBind[1, 6] <- "1"
toBind[1, 8] <- paste0(toBind[1, 2], "_R1")

allD <- unique(rbind(allD, toBind))

##
# GNAI1 <- PRKACA

toBind <- matrix(data = , nrow = 1, ncol = 8)
colnames(toBind) <- colnames(allD)

toBind[1, 1] <- "PRKACA"
toBind[1, 2] <- "PRKACA"
toBind[1, 3] <- "GNAI1"
toBind[1, 4] <- "GNAI1"
toBind[1, 5] <- "R"
toBind[1, 6] <- "1"
toBind[1, 8] <- paste0(toBind[1, 2], "_R1")

allD <- unique(rbind(allD, toBind))

##
# PTGER4 <- GNAI

toBind <- matrix(data = , nrow = 1, ncol = 8)
colnames(toBind) <- colnames(allD)

toBind[1, 1] <- "GNAI1"
toBind[1, 2] <- "GNAI1"
toBind[1, 3] <- "PTGER4"
toBind[1, 4] <- "PTGER4"
toBind[1, 5] <- "R"
toBind[1, 6] <- "1"
toBind[1, 8] <- paste0(toBind[1, 2], "_R1")

##
# Other suggestions from the file (from Anna-Mari)

kinases <- c("PTGER1", "GNAQ", "PLCB1", "PRKCE", "Pyk2FAK2", "Shc", "GRB2", "SOS1", "HRAS", "PIK3CA", "PIK3CA", "PDK1", "GNAQ", "BTK", "GNAQ", "LARG", "RHOA")
substrates <- c("GNAQ", "PLCB1", "PRKCE", "Pyk2FAK2", "Shc", "GRB2", "SOS1", "HRAS", "PIK3CA", "AKT1", "PDK1", "AKT1", "BTK", "PLCG1", "LARG", "RHOA", "ROCK1")

kinases <- c(kinases, c("PTGER2", "GNAS", "GNAS", "adenylate", "CAMP", "PDZGEF1", "HRAS", "RAF1", "MAP2K1", "MAP2K2", "CAMP", "CAMPGEF", "RAP1A", "BRAF", "BRAF", "RAP1A", "RAF1", "RAF1", "CAMP", "PRKACA", "PRKACA"))
substrates <- c(substrates, c("GNAS", "SRC", "adenylate", "CAMP", "PDZGEF1", "HRAS", "RAF1", "MAP2K1", "MAPK3", "MAPK3", "CAMPGEF", "RAP1A", "BRAF", "MAP2K1", "MAP2K2", "RAF1", "MAP2K1", "MAP2K2", "PRKACA", "GSK3A", "GSK3B"))

kinases <- c(kinases, c("PTGER3", "GNAI1", "GNAI1", "RAP1GAP1", "GNAI1", "SRC", "SRC", "SHC", "GRB2", "SOS1", "HRAS", "GNAI1", "GNA12", "GNA12", "PLCE1", "GNA12", "RASA2", "MRAS", "MRGEF", "RASA2", "RRAS", "RASA2", "TC21", "GNA12", "LBC", "LBC", "RAC1", "MEKK1", "MEK4", "MEK4", "LBC", "CDC42", "CDC42", "PAK1", "LIMK1", "LBC", "RHOA", "ROCK1", "PTGER3", "PTGER3", "PTGER3", "PLCB1", "PLCB1", "PLCG1", "PLCG1"))
substrates <- c(substrates, c("GNAI1", "PRKACA", "RAP1GAP1", "RAP1A", "SRC", "STAT3", "SHC", "GRB2", "SOS1", "HRAS", "RAF1", "GNA12", "BTK", "PLCE1", "RAP1A", "RASA2", "MRAS", "MRGEF", "RAP1A", "RRAS", "RAF1", "TC21", "PIK3CA", "LBC", "PRKACA", "RAC1", "MEKK1", "MEK4", "MAP2K1", "MAP2K2", "CDC42", "MEKK1", "PAK1", "LIMK1", "RHOGTP", "RHOA", "ROCK1", "RHOGTP", "GNAS", "GNAQ", "GNA12", "cPRKC", "CALM1", "cPRKC", "CALM1"))

kinases <- c(kinases, c("PTGER4", "PTGER3", "PTGER4", "GBG", "BTK", "GBG", "PLCB1", "GBG", "GBG", "GBG", "PTGER3", "PTGER4", "ARRB", "ARRB", "ARRB", "MAP3K5", "MAP2K4", "ARRB", "RAF1", "ARRB"))
substrates <- c(substrates, c("GNAS", "GBG", "GBG", "BTK", "PLCG1", "PLCB1", "IP3DAG", "PIK3CA", "adenylate", "SRC", "ARRB", "ARRB", "SRC", "MAPK10", "MAP3K5", "MAP2K4", "MAPK10", "RAF1", "MEK1", "EGFR"))

kinases <- c(kinases, c("EGFR", "PLC", "PKC", "PLC", "EGFR", "EGFR", "MTOR", "EGFR", "CBL", "NCK", "EGFR", "VAV1", "RAC1", "RAC1", "EGFR", "JAK1", "STAT1", "SRC"))
substrates <- c(substrates, c("PLC", "PKC", "IKK", "CAMK", "HRAS", "PIK3CA", "S6K", "CBL", "NCK", "PAK1", "VAV1", "RAC1", "MEKK1", "RHOA", "JAK1", "STAT1", "STAT3", "EGFR"))

kinases <- c(kinases, c("PTGER1", "PTGER2", "PTGER3", "PTGER4", "PTGER1", "PTGER2", "PTGER3", "PTGER3", "PTGER4", "PTGER4", "PTGER4", "GNAQ"))
substrates <- c(substrates, c("ARRB", "ARRB", "ARRB", "ARRB", "GNAQ", "GNAS", "GNAI", "GBG", "GNAS", "GNAI", "GBG", "PKC"))

kinases <- c(kinases, c("CALM1", "CALM1", "CALM1", "CALM1", "CALM1", "ARRB", "GNAQ"))
substrates <- c(substrates, c("CAMK2A", "CAMK2B", "CAMK2C", "CAMK2D", "CAMK2G", "CSNK2A1", "PLCG1"))

toBind <- matrix(data = , nrow = length(kinases), ncol = 8)
colnames(toBind) <- colnames(allD)

toBind[, 3] <- kinases
toBind[, 4] <- kinases
toBind[, 2] <- substrates
toBind[, 1] <- substrates
toBind[, 5] <- "R"
toBind[, 6] <- "1"
toBind[, 8] <- paste0(substrates, "_R1")

allD <- unique(rbind(allD, toBind))
##
# Fixing GNG12 & GNB3
allD[which(allD[, 2]=="GNG12"), 2] <- "GBG"
allD[which(allD[, 4]=="GNG12"), 4] <- "GBG"
allD[which(allD[, 8]=="GNG12_R1"), 8] <- "GBG_R1"

allD[which(allD[, 2]=="GNB1"), 2] <- "GBG"
allD[which(allD[, 4]=="GNB1"), 4] <- "GBG"
allD[which(allD[, 2]=="GNB2"), 2] <- "GBG"
allD[which(allD[, 4]=="GNB2"), 4] <- "GBG"
allD[which(allD[, 2]=="GNB3"), 2] <- "GBG"
allD[which(allD[, 4]=="GNB3"), 4] <- "GBG"
allD[which(allD[, 8]=="GNB1_R1"), 8] <- "GBG_R1"
allD[which(allD[, 8]=="GNB2_R1"), 8] <- "GBG_R1"
allD[which(allD[, 8]=="GNB3_R1"), 8] <- "GBG_R1"

allD <- unique(allD)

# allD[, 7] <- paste0("e", 1:nrow(allD))

allD <- as.data.frame(x = allD)
allD$S.AC <- as.character(allD$S.AC)
allD$S.ID <- as.character(allD$S.ID)
allD$K.AC <- as.character(allD$K.AC)
allD$K.ID <- as.character(allD$K.ID)
allD$res <- as.character(allD$res)
allD$pos <- as.character(allD$pos)
allD$SID <- as.character(allD$SID)
allD$S.cc <- as.character(allD$S.cc)

idx <- which(allD$S.ID%in%c("PTGER1", "PTGER2", "PTGER3", "PTGER4"))
allD <- allD[-idx, ]

##
# Grouping
# cPRKC
toReplace <- c("PRKCA", "PRKCB", "PRKCG")
idx1 <- which(allD$K.ID%in%toReplace)
idx2 <- which(allD$S.ID%in%toReplace)
allD$K.ID[idx1] <- "cPRKC"
allD$K.AC[idx1] <- "cPRKC"
allD$S.ID[idx2] <- "cPRKC"
allD$S.AC[idx2] <- "cPRKC"
for(ii in 1:length(idx2)){
  
  allD$S.cc[idx2[ii]] <- paste0("cPRKC", "_", strsplit(x = allD$S.cc[idx2[ii]], split = "_", fixed = TRUE)[[1]][2])
  
}

# nPRKC
toReplace <- c("PRKCD", "PRKCE", "PRKCH", "PRKCQ")
idx1 <- which(allD$K.ID%in%toReplace)
idx2 <- which(allD$S.ID%in%toReplace)
allD$K.ID[idx1] <- "nPRKC"
allD$K.AC[idx1] <- "nPRKC"
allD$S.ID[idx2] <- "nPRKC"
allD$S.AC[idx2] <- "nPRKC"
for(ii in 1:length(idx2)){
  
  allD$S.cc[idx2[ii]] <- paste0("nPRKC", "_", strsplit(x = allD$S.cc[idx2[ii]], split = "_", fixed = TRUE)[[1]][2])
  
}

# aPRKC
toReplace <- c("PRKCZ", "PRKCI")
idx1 <- which(allD$K.ID%in%toReplace)
idx2 <- which(allD$S.ID%in%toReplace)
allD$K.ID[idx1] <- "aPRKC"
allD$K.AC[idx1] <- "aPRKC"
allD$S.ID[idx2] <- "aPRKC"
allD$S.AC[idx2] <- "aPRKC"
for(ii in 1:length(idx2)){
  
  allD$S.cc[idx2[ii]] <- paste0("aPRKC", "_", strsplit(x = allD$S.cc[idx2[ii]], split = "_", fixed = TRUE)[[1]][2])
  
}

# PI3K
toReplace <- c("PIK3CA", "PIK3R1")
idx1 <- which(allD$K.ID%in%toReplace)
idx2 <- which(allD$S.ID%in%toReplace)
allD$K.ID[idx1] <- "PI3K"
allD$K.AC[idx1] <- "PI3K"
allD$S.ID[idx2] <- "PI3K"
allD$S.AC[idx2] <- "PI3K"
for(ii in 1:length(idx2)){
  
  allD$S.cc[idx2[ii]] <- paste0("PI3K", "_", strsplit(x = allD$S.cc[idx2[ii]], split = "_", fixed = TRUE)[[1]][2])
  
}

# MAP3K1
toReplace <- c("MEKK1", "MAPKKK1", "MEKK")
idx1 <- which(allD$K.ID%in%toReplace)
idx2 <- which(allD$S.ID%in%toReplace)
allD$K.ID[idx1] <- "MAP3K1"
allD$K.AC[idx1] <- "MAP3K1"
allD$S.ID[idx2] <- "MAP3K1"
allD$S.AC[idx2] <- "MAP3K1"
for(ii in 1:length(idx2)){
  
  allD$S.cc[idx2[ii]] <- paste0("MAP3K1", "_", strsplit(x = allD$S.cc[idx2[ii]], split = "_", fixed = TRUE)[[1]][2])
  
}

# MAP2K1
toReplace <- c("MEK1", "PRKMK1")
idx1 <- which(allD$K.ID%in%toReplace)
idx2 <- which(allD$S.ID%in%toReplace)
allD$S.ID[idx2] <- "MAP2K1"
allD$S.AC[idx2] <- "MAP2K1"
for(ii in 1:length(idx2)){
  
  allD$S.cc[idx2[ii]] <- paste0("MAP2K1", "_", strsplit(x = allD$S.cc[idx2[ii]], split = "_", fixed = TRUE)[[1]][2])
  
}

# MAP2K4
toReplace <- c("MEK4", "PRKMK1")
idx1 <- which(allD$K.ID%in%toReplace)
idx2 <- which(allD$S.ID%in%toReplace)
allD$K.ID[idx1] <- "MAP2K4"
allD$K.AC[idx1] <- "MAP2K4"
allD$S.ID[idx2] <- "MAP2K4"
allD$S.AC[idx2] <- "MAP2K4"
for(ii in 1:length(idx2)){
  
  allD$S.cc[idx2[ii]] <- paste0("MAP2K4", "_", strsplit(x = allD$S.cc[idx2[ii]], split = "_", fixed = TRUE)[[1]][2])
  
}

# CAMKii
toReplace <- c("CAMK2A", "CAMK2B", "CAMK2C", "CAMK2D", "CAMK2G", paste0("CAMK", 1:4), paste0("CAMKK", 1:2))
idx1 <- which(allD$K.ID%in%toReplace)
idx2 <- which(allD$S.ID%in%toReplace)
allD$K.ID[idx1] <- "CAMKii"
allD$K.AC[idx1] <- "CAMKii"
allD$S.ID[idx2] <- "CAMKii"
allD$S.AC[idx2] <- "CAMKii"
for(ii in 1:length(idx2)){
  
  allD$S.cc[idx2[ii]] <- paste0("CAMKii", "_", strsplit(x = allD$S.cc[idx2[ii]], split = "_", fixed = TRUE)[[1]][2])
  
}

##
# Fix ARRB residues
toReplace <- c("ARRB1", "ARRB2")
idx <- which(allD$S.ID%in%toReplace)
allD$S.ID[idx] <- "ARRB"
allD$S.AC[idx] <- "ARRB"
for(ii in 1:length(idx)){
  allD$S.cc[idx[ii]] <- paste0("ARRB_", strsplit(x = allD$S.cc[idx[ii]], split = "_", fixed = TRUE)[[1]][2])
}
idx = which(allD$K.ID%in%toReplace)
allD$K.AC[idx] = "ARRB"
allD$K.ID[idx] = "ARRB"

##
# Fix GRK's
toReplace <- paste0("GRK", 1:7)
idx <- which(allD$K.ID%in%toReplace)
allD$K.ID[idx] <- "GRK"
idx <- which(allD$S.ID%in%toReplace)
allD$S.ID[idx] <- "GRK"
for(ii in 1:length(idx)){
  allD$S.cc[idx[ii]] <- paste0("GRK_", strsplit(x = allD$S.cc[idx[ii]], split = "_", fixed = TRUE)[[1]][2])
}

##
# Fix GNAI's
idx <- which(grepl(pattern = "GNAI", x = allD$S.ID))
allD$S.ID[idx] = "GNAI"
idx <- which(grepl(pattern = "GNAI", x = allD$K.ID))
allD$K.ID[idx] = "GNAI"
idx <- which(grepl(pattern = "GNAI", x = allD$S.cc))
allD$S.cc[idx] = "GNAI_R1"

idx2rem = which(duplicated(allD[, c(4, 8)]))
allD = allD[-idx2rem, ]

connection = t(as.matrix(c("PRKACA", "PRKACA", "GNAI", "GNAI", "R", "1", NA, "PRKACA_R1")))
colnames(connection) = colnames(allD)
connection = as.data.frame(connection)

allD = unique(rbind(allD, connection))
##
# remove GNAS -> PI3K
idx1 <- which(allD$K.ID=="GNAS")
idx2 <- which(allD$S.ID=="PI3K")
idx <- intersect(x = idx1, y = idx2)
allD <- allD[-idx, ]

# remove ARRB -> EGFR
idx1 <- which(allD$K.ID=="ARRB")
idx2 <- which(allD$S.ID=="EGFR")
idx <- intersect(x = idx1, y = idx2)
allD <- allD[-idx, ]

allD$SID = paste0("e", 1:nrow(allD))

idx1 = which(x = grepl(pattern = "GNAI", x = allD$S.AC, fixed = TRUE))
allD$S.AC[idx1] = "GNAI"
allD$S.ID[idx1] = "GNAI"
idx2 = which(x = grepl(pattern = "GNAI", x = allD$K.AC, fixed = TRUE))
allD$K.AC[idx2] = "GNAI"
allD$K.ID[idx2] = "GNAI"

# remove redundant GNAS
idx1 = which(x = grepl(pattern = "GNAS", x = allD$S.AC))
idx2 = which(x = grepl(pattern = "PTGER", x = allD$K.AC))
idx = intersect(x = idx1, y = idx2)
allD = allD[-setdiff(x = idx1, y = idx), ]

# remove redundant GNAQ
idx1 = which(x = grepl(pattern = "GNAQ", x = allD$S.AC))
idx2 = which(x = grepl(pattern = "PTGER", x = allD$K.AC))
idx = intersect(x = idx1, y = idx2)
allD = allD[-setdiff(x = idx1, y = idx), ]

# remove redundant GNAI
idx1 = which(x = grepl(pattern = "GNAI", x = allD$S.AC))
idx2 = which(x = grepl(pattern = "PTGER", x = allD$K.AC))
idx = intersect(x = idx1, y = idx2)
allD = allD[-setdiff(x = idx1, y = idx), ]

# remove redundant GNA12
idx1 = which(x = grepl(pattern = "GNA12", x = allD$S.AC))
idx2 = which(x = grepl(pattern = "PTGER", x = allD$K.AC))
idx = intersect(x = idx1, y = idx2)
allD = allD[-setdiff(x = idx1, y = idx), ]

# Now grouping all the AKT's
allD$S.AC[which(allD$S.AC%in%paste0("AKT", 1:3))] = "AKT"
allD$S.ID[which(allD$S.ID%in%paste0("AKT", 1:3))] = "AKT"
allD$K.AC[which(allD$K.AC%in%paste0("AKT", 1:3))] = "AKT"
allD$K.ID[which(allD$K.ID%in%paste0("AKT", 1:3))] = "AKT"
for(ii in 1:nrow(allD)){
  
  pp = strsplit(x = allD$S.cc[ii], split = "_", fixed = TRUE)[[1]][1]
  ss = strsplit(x = allD$S.cc[ii], split = "_", fixed = TRUE)[[1]][2]
  
  if(pp%in%paste0("AKT", 1:3)){
    
    allD$S.cc[ii] = paste0("AKT_", ss)
    
  }
  
}
allD$SID = "e1"
allD = unique(allD)
allD$SID = paste0("e", 1:nrow(allD))

save(allD, file = "allD.RData")
