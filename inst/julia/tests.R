fid_vertex <- function(VT1,CC1,VTsum,U,L,m,Dim,k,n){
  whichVT <- matrix(NA, nrow=m, ncol=2) 
  whichVT[,1] <- VTsum>=L[k]
  whichVT[,2] <- VTsum<=U[k]
  whichl <- which(whichVT[,1])  #vertices that satisfy lower constraint
  whichu <- which(whichVT[,2])  #vertices that satisfy upper constraint
  both <- which(whichVT[,1] & whichVT[,2])
  checkl <- which(!whichVT[,1])
  checku <- which(!whichVT[,2])
  CCtemp <- NULL  
  VTtemp <- NULL
  vert <- 0 
  CA <- as.matrix(CC1[,checkl])
  CB <- as.matrix(CC1[,whichl]) 
  if(length(checkl)>0){ #i.e. they do not all satisfy the lower constraint
    INT  <-  matrix(0,nrow=2*n,ncol=length(checkl)) 
    for(ll in 1:length(checkl)){ #check lower constraints first
      INT[CA[,ll],ll] <- 1 
    }
    for(ii in 1:length(whichl)){
      INT2 <- as.matrix(INT[CB[,ii],]) # ok si une seule ligne ??? 
      use <- which(colSums(INT2)==(Dim-1)) # sum sur colonne
      for(dd in use){
        inter <- CB[which(INT2[,dd]==1),ii]  #this will be intersection
        vert <- vert+1 
        CCtemp <- cbind(CCtemp, c(inter,k+n))  # need to add n indicating lower constraint
        lambda <- (L[k]-VTsum[whichl[ii]])/(VTsum[checkl[dd]]-VTsum[whichl[ii]]) 
        VTtemp <- cbind(VTtemp,lambda*VT1[,checkl[dd]]+(1-lambda)*VT1[,whichl[ii]]) 
      }
    }
  }
  CA <- as.matrix(CC1[,checku]) 
  CB <- as.matrix(CC1[,whichu]) 
  if(length(checku)>0){ #i.e. they do not all satisfy the lower constraint
    INT  <-  matrix(0,nrow=2*n,ncol=length(checku)) 
    for(ll in 1:length(checku)){ #check lower constraints first
      INT[CA[,ll],ll] <- 1 
    }
    for(ii in 1:length(whichu)){
      INT2 <- as.matrix(INT[CB[,ii],]) 
      use <- which(colSums(INT2)==(Dim-1)) # sum sur colonne
      for(dd in use){
        inter <- CB[which(INT2[,dd]==1),ii]  #this will be intersection
        vert <- vert+1 
        CCtemp <- cbind(CCtemp, c(inter,k))  # need to add n indicating lower constraint
        lambda <- (U[k]-VTsum[whichu[ii]])/(VTsum[checku[dd]]-VTsum[whichu[ii]]) 
        VTtemp <- cbind(VTtemp,lambda*VT1[,checku[dd]]+(1-lambda)*VT1[,whichu[ii]]) 
      }
    }
  }
  if(length(both)>0){ 
    for(ll in both){
      vert <- vert+1 
      CCtemp <- cbind(CCtemp,CC1[,ll]) 
      VTtemp <- cbind(VTtemp,VT1[,ll]) 
    }
  }
  return(list(VTtemp=VTtemp,CCtemp=CCtemp,vert=vert))
}


VT1 = rbind(
  c(1.0, 2, 3, 4, 5, 6, 7, 8),
  c(8, 7, 6, 5, 4, 3, 2, 1),
  c(5, 5, 6, 7, 7, 4, 5, 4)
)
CC1 = rbind(
  c(3, 3, 2, 2, 1, 1, 2, 2), # plante si head = 3 3
  c(2, 1, 2, 1, 2, 1, 2, 1),
  c(1, 1, 1, 1, 2, 2, 2, 2)
)
VTsum = c(1.0, 2, 4, 5, 9, 3, 5, 7)
U = c(6.0,6)
L = c(2.0,2)
m = 8
dim = 3
k = 2
n = 18

VT1 = rbind(
  c(5.665427347856155, 5.665427347856155, 5.665427347856156, 5.665427347856155, 5.665427347856155, 5.665427347856156),
  c(1.0570980348004602, 1.0570980348004604, 1.0570980348004608, 1.0570980348004602, 1.0570980348004604, 1.0570980348004608),
  c(3.9799878512017155e-16, 3.1211045504395e-16, 0.0, 3.9799878512017155e-16, 3.1211045504395e-16, -0.0)
)

VT1 = rbind(
    c(1.0, 2, 3, 4, 5, 6, 7, 8),
    c(8, 7, 6, 5, 4, 3, 2, 1),
    c(5, 5, 6, 7, 7, 4, 5, 4)
)

fid_vertex(VT1, CC1, VTsum, U, L, m, dim, k, n)
