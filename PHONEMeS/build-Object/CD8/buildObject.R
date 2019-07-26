library(readxl)
CD8 <- read_excel("../CD8_New.xlsx")

library(readr)
uniprot <- read_delim("uniprot.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

rNames <- CD8$`Unique identifier`

###
# Preparing mapping table

data.IDmap <- matrix(, nrow = length(rNames), ncol = 3)
colnames(data.IDmap) <- c("dataID", "UPID", "S.cc")
data.IDmap[, 1] <- rNames
for(i in 1:nrow(data.IDmap)){
  
  protein <- CD8$`Gene names`[i]
  site <- CD8$`Amino acid`[i]
  position <- CD8$Position[i]
  
  idx <- which(uniprot$Entry==protein)
  
  if(length(idx) > 0){
    
    data.IDmap[i, 2] <- uniprot$`Entry name`[idx[1]]
    data.IDmap[i, 3] <- paste0(uniprot$`Entry name`[idx[1]], "_", site, position)
    
  } else {
    
    data.IDmap[i, 2] <- paste0(protein)
    data.IDmap[i, 3] <- paste0(protein, "_", site, position)
    
  }
  
}

data.IDmap[, 1] = data.IDmap[, 3]

##
ep1Val <- CD8[, 6:10]
ep2Val <- CD8[, 11:15]
ep3Val <- CD8[, 16:20]
ep4Val <- CD8[, 21:25]
pgeVal <- CD8[, 26:30]

pValThresh = 0.2

ep1Score <- log2(CD8$`Student's T-test q-value EP1_CT`/pValThresh)

ep2Score <- log2(CD8$`Student's T-test q-value EP2_CT`/pValThresh)

ep3Score <- log2(CD8$`Student's T-test q-value EP3_CT`/pValThresh)

ep4Score <- log2(CD8$`Student's T-test q-value EP4_CT`/pValThresh)

pgeScore <- log2(CD8$`Student's T-test q-value PGE2_CT`/pValThresh)

##
ep1.p <- CD8$`Student's T-test q-value EP1_CT`
ep1.lo <- ep1Score
ep1.c <- rep("C", length(ep1.p))
ep1.c[which(ep1.lo<0)] <- "P"
ep1.s <- rep("OK", length(ep1.p))
ep1.fc <- rowMeans(x = ep1Val)

##
ep2.p <- CD8$`Student's T-test q-value EP2_CT`
ep2.lo <- ep2Score
ep2.c <- rep("C", length(ep2.p))
ep2.c[which(ep2.lo<0)] <- "P"
ep2.s <- rep("OK", length(ep2.p))
ep2.fc <- rowMeans(x = ep2Val)

##
ep3.p <- CD8$`Student's T-test q-value EP3_CT`
ep3.lo <- ep3Score
ep3.c <- rep("C", length(ep3.p))
ep3.c[which(ep3.lo<0)] <- "P"
ep3.s <- rep("OK", length(ep3.p))
ep3.fc <- rowMeans(x = ep3Val)

##
ep4.p <- CD8$`Student's T-test q-value EP4_CT`
ep4.lo <- ep4Score
ep4.c <- rep("C", length(ep4.p))
ep4.c[which(ep4.lo<0)] <- "P"
ep4.s <- rep("OK", length(ep4.p))
ep4.fc <- rowMeans(x = ep4Val)

##
pge.p <- CD8$`Student's T-test q-value PGE2_CT`
pge.lo <- pgeScore
pge.c <- rep("C", length(pge.p))
pge.c[which(pge.lo<0)] <- "P"
pge.s <- rep("OK", length(pge.p))
pge.fc <- rowMeans(x = pgeVal)

###
# Preparing GMM objects
GMM<-vector("list", length = length(rNames))
names(GMM)<-data.IDmap[, 1]
for(i in 1:length(GMM)){
  GMM[[i]]<-rbind(
    c(as.character(ep1.lo[i]), as.character(ep1.c[i]), as.character(ep1.p[i]), as.character(ep1.s[i])),
    c(as.character(ep2.lo[i]), as.character(ep2.c[i]), as.character(ep2.p[i]), as.character(ep2.s[i])),
    c(as.character(ep3.lo[i]), as.character(ep3.c[i]), as.character(ep3.p[i]), as.character(ep3.s[i])),
    c(as.character(ep4.lo[i]), as.character(ep4.c[i]), as.character(ep4.p[i]), as.character(ep4.s[i])),
    c(as.character(pge.lo[i]), as.character(pge.c[i]), as.character(pge.p[i]), as.character(pge.s[i]))
  )
  colnames(GMM[[i]])<-c("Indiv", "clus","FCvCaPval","status")
  rownames(GMM[[i]])<-c("EP1", "EP2", "EP3", "EP4", "PGE2")
}

GMM.wFC<-vector("list", length = length(rNames))
names(GMM.wFC)<-data.IDmap[, 1]
for(i in 1:length(GMM.wFC)){
  GMM.wFC[[i]]<-rbind(
    c(as.character(ep1.lo[i]), as.character(ep1.c[i]), as.character(ep1.p[i]), as.character(ep1.s[i]), as.character(ep1.fc[i])),
    c(as.character(ep2.lo[i]), as.character(ep2.c[i]), as.character(ep2.p[i]), as.character(ep2.s[i]), as.character(ep2.fc[i])),
    c(as.character(ep3.lo[i]), as.character(ep3.c[i]), as.character(ep3.p[i]), as.character(ep3.s[i]), as.character(ep3.fc[i])),
    c(as.character(ep4.lo[i]), as.character(ep4.c[i]), as.character(ep4.p[i]), as.character(ep4.s[i]), as.character(ep4.fc[i])),
    c(as.character(pge.lo[i]), as.character(pge.c[i]), as.character(pge.p[i]), as.character(pge.s[i]), as.character(pge.fc[i]))
  )
  colnames(GMM.wFC[[i]])<-c("Indiv", "clus","FCvCaPval","status", "avgVal")
  rownames(GMM.wFC[[i]])<-c("EP1", "EP2", "EP3", "EP4", "PGE2")
}

###
# Saving GMM object as a list

GMM.ID <- data.IDmap
colnames(GMM.ID) <- c("dataID", "UPID", "S.cc")
GMM.ID <- as.data.frame(GMM.ID)

save(list=c("GMM.ID", "GMM","GMM.wFC"), file="dataGMM_CD8.RData")