#setwd("D:/Work/R/Fiducial")

# update: use assign() and get() instead of list
# done for VT - checking ok

# � faire: remplacer des listes avec aasign() et get()

# pour gram-schmidt taper "orthonormalization" dans RSiteSearch

#######################################################
#######    AUXILIARY FUNCTIONS    #####################
####################################################### 
rankM <- Matrix::rankMatrix  # or pracma::mrank
#rankM <- pracma::mrank
null <- function(M) MASS::Null(t(M))
#null <- function(M) t(as.matrix(pracma::nullspace(t(M))))
orth <- pracma::orth
#linsolve <- function(A,b) MASS::ginv(t(A)%*%A)%*%t(A)%*%b  
#linsolve <- function(A,b) qr.solve(A,b,tol=1e-5)
linsolve <- function(A,b) pracma::mldivide(A,b)
#linsolve <-  function(A,b) lsfit(A,b,intercept=FALSE)$coefficients

fid_sample <- function(VT2,VTsum,U,L){
  high <- which(VT2>0)
  low <- which(VT2<0)
  zero <- which(VT2==0)
  if((length(high)>0 & length(low)>0) | length(zero)>0){
    UU <- sign(U-VTsum)
    LL <- sign(L-VTsum)
    SS <- sign(VT2)
    zeronot <- which(VT2!=0)
    if(length(zero)==length(VT2)){
      MAX <- Inf
      MIN <- -Inf
      temp <- any(UU>0)*any(LL<0)
    }else{
      c1 <- !any(UU[zero]==-1) & !any(LL[zero]==-1)
      c2 <- !any(UU[zero]==1) & !any(LL[zero]==1)
      d1 <- !any(SS==-1)
      d2 <- !any(SS==1)
      if((d1 & c1)|(d2 & c2)){
        Us <- (U-VTsum[zeronot])/VT2[zeronot]
        Ls <- (L-VTsum[zeronot])/VT2[zeronot]
        MAX <- Inf
        MIN <- min(Us,Ls)
        temp <- 1-(atan(MIN)/pi+0.5)
      }else if((d2 & c1)|(d1 & c2)){
        Us <- (U-VTsum[zeronot])/VT2[zeronot]
        Ls <- (L-VTsum[zeronot])/VT2[zeronot]
        MAX <- max(Us,Ls)
        MIN <- -Inf
        temp <- atan(MAX)/pi+0.5 # il y avait une faute ici				
      }else{
        HUs <- (U-VTsum[high])/VT2[high]
        HLs <- (L-VTsum[high])/VT2[high]
        Hmax=max(HUs,HLs)
        Hmin=min(HUs,HLs)
        LUs <- (U-VTsum[low])/VT2[low]
        LLs <- (L-VTsum[low])/VT2[low]
        Lmax=max(LUs,LLs)
        Lmin=min(LUs,LLs)
        if(round(Lmin-Hmax,12)>=0){
          bpos <- -Inf
          tpos <- Hmax
          bneg <- Lmin
          tneg <- Inf
        }else if(round(Hmin-Lmax,12)>=0){
          bpos <- Hmin
          tpos <- Inf
          bneg <- -Inf
          tneg <- Lmax
        }else{
          bpos <- -Inf
          tpos <- Inf
          bneg <- -Inf
          tneg <- Inf
        }
        if(tpos==Inf){
          Pprob <- 1-(atan(bpos)/pi+0.5)
        }else{
          Pprob <- atan(tpos)/pi+0.5
        }
        if(tneg==Inf){
          Nprob <- 1-(atan(bneg)/pi+0.5)
        }else{
          Nprob <- atan(tneg)/pi+0.5
        }
        temp <- Pprob+Nprob
        Pprob <- Pprob/temp
        Nprob <- 1-Pprob
        if(runif(1)<=Nprob){
          MIN <- bneg
          MAX <- tneg
        }else{
          MIN <- bpos
          MAX <- tpos
        }
      }
    }
    y <- atan(MAX)/pi+0.5
    x <- atan(MIN)/pi+0.5
    u <- x+(y-x)*runif(1)
    ZZ <- tan(pi*(u-0.5))
    wt <- exp(-ZZ^2/2)*(1+ZZ^2)*temp
  }else{
    Us <- (U-VTsum)/VT2
    Ls <- (L-VTsum)/VT2
    MAX <- max(Us,Ls)
    MIN <- min(Us,Ls)
    y<-atan(MAX)/pi+.5; #cdf
    x<-atan(MIN)/pi+.5
    u<-x+(y-x)*runif(1)
    ZZ<-tan(pi*(u-.5)) #Inverse cdf
    wt<-exp(-ZZ^2/2)*(1+ZZ^2)*(y-x)
  }
  return(c(ZZ=ZZ,wt=wt))
}

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

#library(rgr)
inference <- function(vertex, weight, alpha=0.05){ 
  out <- rep(NA,4)
  names(out) <- c("mean","median","low","up")
  out[1] <- sum(vertex*weight)
  h <- cbind(vertex,weight)
  hsort <- gx.sort(h,1)
  hsum <- cumsum(hsort[,2])
  ci_u <- min(which(hsum>=1-alpha/2)) #upper confidence bound
  ci_l <- min(which(hsum>=alpha/2))   #lower confidence bound
  ci_m <- min(which(hsum>=0.5))
  out[3] <- hsort[ci_l,1]  #lower bound
  out[4] <- hsort[ci_u,1] #upper bound
  out[2] <- hsort[ci_m,1] #estimate
  out
}





######################################################
######## OUTPUT ######################################
######################################################
fid_nMLM <- function(dat, FE, RE, N, thresh = N/2){
  n <- nrow(dat)
  random_design  <-  1 
  Y  <-  dat 
  L <- Y[,1] 
  U <- Y[,2] 
  
  fe <- ncol(FE) 
  
  break_point <- 10
  
  ########--------SET-UP RANDOM EFFECTS
  # RE is declared.
  re <- ncol(RE)+1 
  E <- rep(NA,re) # E[i] number of levels of i-th random effect
  E[re] <- n 
  for(i in 1:(re-1)){
    E[i] <- length(levels(RE[,i])) 
  }
  ESUM <- cumsum(E)
  RE2  <-  cbind(RE,factor(1:n))   #Adds the error effect
  RE <-  NULL 
  for(i in 1:re){ #Builds an indicator RE matrix for the effects
    re_levels <- levels(RE2[,i])
    for(j in 1:E[i]){
      temp1 <- which(RE2[,i]==re_levels[j]) 
      temp2 <- rep(0,n) 
      temp2[temp1] <- 1 
      RE <- cbind(RE,temp2)
    }
  } #  
  
  Dim <- fe+re  #Dimension of the space
  repp <- 1 
  
  
  
  ########--------SET-UP ALGORITHM OBJECTS
  Z <-  weight <- vector(mode = "list", length = re)  #Particles and Weights
  ESS <- rep(N,n)    #Effective sample size
  VC <- rep(0,N)    #Number of vertices
  CC <- C <- vector(mode = "list", length = N)       # VT Verticies ; CC Constraints ; C initial constraints
  # je crois que seul C[[1]] intervient dans la suite 
  VTnames <- paste("VT.", 1:N, sep="")
  
  ########--------SAMPLE ALL Z's/SET-UP WEIGHTS
  
  A <- vector(mode = "list", length = N) 
  for(i in 1:re){  #Sample all the Z's and set-up the weights
    Z[[i]] <- matrix(rnorm(E[i]*N), ncol=N) # m�me nb de colonnes on pourrait empiler en matrice - oui mais pas dans la suite !!
    weight[[i]]  <- matrix(1, nrow=E[i], ncol=N)
    for(j in 1:N){
      A[[j]]= cbind(A[[j]], RE[ , (ESUM[i]-E[i]+1):ESUM[i]] %*% Z[[i]][,j])
    } # chaque A[[j]] est une matrice n x re (c'est une simulation des random effects)
  }     # donc on pourrait faire un tableau � 3 dimensions 
  
  #temp=[]; #Initial constraints selected
  for(i in 1:1){ # input A[[1]] uniquement      
    #AA <- cbind(rbind(FE,-FE), rbind(A[[i]],-A[[i]]))
    AA <- cbind(FE, A[[i]])	
    AT <- NULL 
    r <-  0   
    if(rankM(AA) != Dim){
      stop("Design is not of full rank")
    }else{
      for(k in 1:n){
        A_temp <- rbind(AT, AA[k,])  # prendre [FE A] au lieu de AA c'est pareil ! 
        if(rankM(A_temp)>r){
          AT <- A_temp
          r <- rankM(AT)
          C[[i]] <- c(C[[i]],k) 
        }
      }
    }
    K <- c(1:n)[!is.element(1:n,C[[1]])] 
    K_start <- C[[1]]  
    Z[[re]][K,] <- 0  #remove error particles not selected in intitialization
    # K will be the values left to be sampled
  } #C[[1]]=K_start sont les indices de Z[[residual]] qu'on garde 
  # construits avec une condtion sur des rangs - il y en a Dim (=rank(AA))
  # les lignes correspondantes de AA forment une matrice de rang Dim
  # => �a ne me semble pas utile ici de concat�ner [FE A] avec [-FE, -A] 
  
  
  
  ########--------FIND INITIAL VERTICES
  # input : C[[1]]=K_start
  # une ligne de USE = une combinaison de Dim indices tels que les lignes correspondantes de [[FE A] ; [-FE -A]] sont libres
  #
  #USE=t(matrix(rep(C[[1]],2^Dim),ncol=2^Dim))  #all combinations of I         ??? C{1,i} dans le code original => i=1 ?
  #for j=1:Dim
  #    if j==1
  #        USE(2^(Dim-1)+1:2^(Dim),1)=USE(1:2^(Dim-1),1)+n 
  #    else
  #        temp=2^(j-2)+1 
  #        while temp<=2^(Dim)
  #            for ii=1:2^(j-2)
  #                USE(temp,j)=USE(temp,j)+n 
  #                temp=temp+1 
  #            end
  #            temp=temp+2^(j-2) 
  #        end
  #    end
  #end
  #
  #L <- vector(mode = "list", length = Dim) 
  #for(i in 1:Dim){
  #	L[[i]] <- c(0,n)+C[[1]][i]
  #}
  USE <- as.matrix(expand.grid(lapply(1:Dim, function(i) c(0,n)+C[[1]][i]), KEEP.OUT.ATTRS=FALSE))
  # maintenant on calcule les sommets de Q1 dans V 
  b <- c(U,-L) 
  for(j in 1:N){
    V <- matrix(NA, nrow=Dim, ncol=2^(Dim)) 
    AA <- cbind(rbind(FE,-FE), rbind(A[[j]],-A[[j]]))
    #temp <- 0 # ??
    for(ii in 1:2^Dim){ #number of vertices
      II <- USE[ii,]  #constraints to use are by row
      V[,ii] <- solve(AA[II,],b[II])  # s�lection dans Q ? c'est la fonction V ; ici r�solution exacte ? oui
    }
    assign(VTnames[j], V) 
  }
  CC <- rep(list(t(USE)),N) #constraints are the same for all N particles
  VC <- rep(2^(Dim),N) # c'est (ncol(CC[[1]], ..., ncol(CC[[2]]) 
  
  
  ########--------MAIN ALGORITHM
  
  K_n=ceiling(length(K)/break_point); # d�coupage de K les valeurs � resampler - pourquoi ?? 
  K_temp=rep(list(NULL),2)
  for(i in seq_len(K_n-1)){
    K_temp[[i]] <- K[((i-1)*break_point+1):(i*break_point)]
  }
  K_temp[[K_n]] <- K[((K_n-1)*break_point+1):length(K)] 
  K1=NULL
  for(k_n in 1:K_n){
    K1 <- c(K1, K_temp[[k_n]])
    for(k in K_temp[[k_n]]){
      if(k_n>1){
        for(i in 1:re){
          Z[[i]] <- rbind(Z[[i]], matrix(rnorm(N*(E[i]-nrow(Z[[i]]))),ncol=N))  
        }
      }
      for(i in 1:N){
        m <- VC[i]
        VT1 <- get(VTnames[i]) # vertex values
        VT2 <- VT1[Dim,]
        VT1 <- VT1[-Dim,] 
        if(fe>0){
          Z1 <- FE[k,]
        }else{
          Z1 <- NULL
        }
        for(j in 1:re){
          if(RE2[k,j]==0){ # ne fait jamais 0 
            Z1 <- c(Z1,0)
          }else{ # ceci est correct
            Z1 <- c(Z1, Z[[j]][RE2[k,j],i])
          }
        }
        Z1 <- Z1[-Dim] # remove column to be sampled
        VTsum <- Z1%*%VT1
        fidsample <- fid_sample(VT2,VTsum,U[k],L[k]) ###Sample	
        ZZ <- fidsample["ZZ"]
        wt <- fidsample["wt"]
        Z[[re]][k,i] <- ZZ
        weight[[re]][k,i] <- wt
        #VTsum.back <- VTsum
        VTsum <- VTsum + ZZ*VT2
        VT1 <- get(VTnames[i]) 
        CC1 <- CC[[i]]
        fidvertex <- fid_vertex(VT1,CC1,VTsum,U,L,m,Dim,k,n)
        VC[i] <- fidvertex$vert # est-ce ncol(CCtemp) ? 
        CC[[i]] <- fidvertex$CCtemp
        assign(VTnames[i], fidvertex$VTtemp)
        if(fidvertex$vert==0){ #check
          weight[[re]][k,i] <- 0
        }
      }
      WT <- apply(weight[[re]],2,cumprod)  #only last re is restricted
      WT <- WT[nrow(WT),]
      if(sum(WT)==0){
        stop("Error: possible underflow")
      }
      WT <- WT/sum(WT)
      ESS[k] <- 1/crossprod(WT)
      #---------------------------------------------------------------Resample
      if(ESS[k]<thresh & k<max(K)){
        u <- rep(0,N) 
        N_sons <- rep(0,N)  
        # generate the cumulative distribution
        dist <- cumsum(WT) 
        aux <- runif(1)    # sample uniform rv in [0 1]
        u <- aux+c(0:(N-1)) 
        u <- u/N 
        j <- 1 
        for(i in 1:N){
          while(u[i]>dist[j]){ # qu'est-ce ?
            j <- j+1 
          }
          N_sons[j] <- N_sons[j]+1 
        }
        JJ <- unique(c(c(1:k),C[[1]])) #  K_start et les entiers <=k
        II <- c(1:n)[-JJ]
        ZZ <- rep(list(NULL),re)			
        VCVC <- rep(0,N) 
        CCCC <- rep(list(NULL),N) 
        VTVT <- rep(list(NULL),N)			
        for(i in 1:N){  ## simplifier �a for(i in which(N_sons[i]>0)){
          if(N_sons[i]>0){
            VCtemp <- VC[i]*rep(1,N_sons[i])
            Ztemp <- rep(list(NULL),re)
            VTtemp <- rep(list(NULL),N_sons[i])
            copy <- N_sons[i]-1  # to be resampled
            for(ii in 1:N_sons[i]){
              VTtemp[[ii]] <- get(VTnames[i])  # copy original vertices
            }
            for(ii in 1:re){
              Ztemp[[ii]] <- Z[[ii]][,rep(i,N_sons[i])]  #copy Z
            }			
            if(copy>0){
              for(rr in 1:repp){ # ? repp=1
                ord <- sample.int(re)  #Order to resample.  Each re will be resampled.
                for(kk in ord){
                  for(ii in 1:copy){
                    XX <- NULL
                    for(jj in 1:re){
                      XX <- cbind(XX,RE[,(ESUM[jj]-E[jj]+1):ESUM[jj]]%*%Ztemp[[jj]][,ii]) 
                    }
                    XX <- XX[JJ,-kk] #remove column of effect to be resampled
                    temp <- which(RE2[JJ,kk]!=0)  #find which levels of kk have been sampled          
                    ### ???? RE2 ne fait jamais 0 !!
                    Z1 <- Ztemp[[kk]][unique(RE2[JJ[temp],kk]),ii]  #Z being resampled
                    CO2 <- RE[JJ,(ESUM[kk]-E[kk]+1):ESUM[kk]] 
                    level0 <- which(colSums(abs(CO2))!=0)
                    CO2 <- CO2[,level0] #levels not sampled yet
                    #Z1==0 <=> Z's not sampled
                    Z00 <- which(Z1!=0) #These are the levels with Z for effect kk
                    Z1 <- Z1[Z00]
                    if(fe>0){
                      XX <- cbind(FE[JJ,],XX) 
                    }
                    MAT <- cbind(-XX,CO2) 						
                    if(rankM(MAT)<ncol(MAT)){
                      NUL <- null(MAT)  
                      n1 <- as.matrix(NUL[1:(nrow(NUL)-ncol(CO2)),])  # ? pas utilis� 
                      n2 <- as.matrix(NUL[(nrow(NUL)-ncol(CO2)+1):nrow(NUL),]) # n1, n2 : eta1 eta2 dans [JC] 
                      # non, n2 pas normalis� 
                      GS <- pracma::gramSchmidt(n2)
                      O2 <- GS$Q
                      O1 <- n1%*%solve(GS$R)
                      # O2 <- orth(n2) 
                      #B  <-  CO2%*%O2
                      #O1 <- linsolve(XX,B) # O1 en colonne dans Octave
                      a <- t(O2)%*%Z1 # a correspond � C dans [JC] mais avec n2 et pas O2
                      b <- sqrt(crossprod(Z1-O2%*%a))[1]
                      tau <- (Z1-O2%*%a)/b
                      rank.O2 <- rankM(O2) # ? = nrow(O2) ?
                      bb <- sqrt(rchisq(1,length(Z1)-rank.O2)) 
                      bbb  <-  b/bb	 # 					
                      aa <- rnorm(rank.O2) 
                      Ztemp[[kk]][Z00,ii] <- O2%*%aa+bb*tau 
                      vert <- c(1:Dim)[-(fe+kk)]
                      for(jj in 1:VC[i]){ 
                        check1 <- XX%*%VTtemp[[ii]][vert,jj]+VTtemp[[ii]][fe+kk,jj]*CO2%*%Z1 
                        VTtemp[[ii]][vert,jj] <- VTtemp[[ii]][vert,jj]- VTtemp[[ii]][fe+kk,jj]*O1%*%(bbb*aa-a) 
                        VTtemp[[ii]][fe+kk,jj] <- VTtemp[[ii]][fe+kk,jj]*b/bb
                        check2 <- XX%*%VTtemp[[ii]][vert,jj]+VTtemp[[ii]][fe+kk,jj]*CO2%*%(O2%*%aa+bb*tau)
                      }
                    }else{
                      b <- sqrt(crossprod(Z1))[1] 
                      tau <- Z1/b 
                      bb <- sqrt(rchisq(1,length(Z1))) 
                      Ztemp[[kk]][Z00,ii] <- bb*tau 
                      vert <- c(1:Dim)[-(fe+kk)] 
                      for(jj in 1:VC[i]){
                        VTtemp[[ii]][fe+kk,jj] <- VTtemp[[ii]][fe+kk,jj]*b/bb 
                      }
                    }
                  }
                }
              }
            }    
            for(ii in 1:re){
              ZZ[[ii]] <- cbind(ZZ[[ii]],Ztemp[[ii]])
            }
            VCVC[(sum(N_sons[1:i-1])+1):sum(N_sons[1:i])] <- VCtemp 
            d <- sum(N_sons[1:i-1])
            for(kk in 1:N_sons[i]){ 
              VTVT[[kk+d]] <- VTtemp[[kk]]
              CCCC[[kk+d]] <- CC[[i]]
            }
          }
        }
        Z <- ZZ
        for(j in 1:N){
          assign(VTnames[j], VTVT[[j]])
        }
        VC <- VCVC
        CC <- CCCC
        weight[[re]] <- matrix(1, nrow=E[re], ncol=N) #assign weights of error matrix to 1
      }
    } # ends reampling for k=K1
    #----------------------------------------------------determine signs
    signs=matrix(0, nrow=re, ncol=N)    
    for(i in 1:N){
      for(j in 1:re){
        if(all(get(VTnames[i])[fe+j,]>0)){ #i.e. all are positive
          signs[j,i] <- 1
        }
        else if(all(get(VTnames[i])[fe+j,]<0)){ #i.e. all are negative
          signs[j,i] <- -1 
        }
      }
    }
    #----------------------------------------------------FINAL RESAMPLE			
    ZZ <- rep(list(NULL),re)
    VTVT <- rep(list(NULL),N)
    n1 <- sort(c(K1,K_start)) # c'est 1:n ? non : pour le dernier k_n oui
    nn <- rep(list(NULL),re) 
    for(ii in 1:re){
      nn[[ii]] <- unique(RE2[n1,ii])  
    }
    lengths.nn <- sapply(nn, length)
    for(i in 1:N){
      Ztemp <- rep(list(NULL),re)
      VTtemp <- get(VTnames[i]) 
      for(ii in 1:re){
        Ztemp[[ii]] <- Z[[ii]][nn[[ii]],i]  #copy Z
      }
      ord <- sample.int(re)
      for(kk in ord){
        #CO=RE 
        XX <- NULL
        eff <- c(1:re)[-kk] 
        for(jj in eff){
          XX <- cbind(XX, RE[,(ESUM[jj]-E[jj]+1):(ESUM[jj]-E[jj]+lengths.nn[jj])]%*%Ztemp[[jj]]) 
        }
        #temp <- which(RE2[,kk]!=0)  #find which levels of kk have been sampled
        # ??? RE2 ne fait jamais 0 ! et temp n'intervient plus
        Z1 <- Ztemp[[kk]]  #Z being resampled
        CO2 <- RE[,(ESUM[kk]-E[kk]+1):(ESUM[kk]-E[kk]+lengths.nn[kk])]  # c'�tait E[:,kk] mais une seule ligne
        #level0 <- which(colSums(abs(CO2))==0) #levels not sampled yet
        # pour le final resample y'a rien 
        ## Z1==0 <=> Z's not sampled
        Z00 <- which(Z1!=0) #These are the levels with Z for effect kk
        #CO2 <- CO2[,-level0]
        Z1 <- Z1[Z00]
        if(fe>0){
          XX <- cbind(FE,XX)
        }
        MAT <- cbind(-XX,CO2)
        if(rankM(MAT)<ncol(MAT)){
          NUL <- null(MAT) 
          n1 <- as.matrix(NUL[1:(nrow(NUL)-ncol(CO2)),])  # ? pas utilis�
          n2 <- as.matrix(NUL[(nrow(NUL)-ncol(CO2)+1):nrow(NUL),])
          GS <- pracma::gramSchmidt(n2)
          O2 <- GS$Q
          O1 <- n1%*%solve(GS$R) # am�liorable forwardsolve(U, x=diag(10)) pour U triang inf 10x10 
          #O2 <- orth(n2)
          #B  <-  CO2%*%O2
          #O1 <- linsolve(XX,B)
          a <- t(O2)%*%Z1
          b <- sqrt(crossprod(Z1-O2%*%a))[1]
          tau <- (Z1-O2%*%a)/b
          rank.O2 <- rankM(O2) # ? = nrow(O2) ?
          bb <- sqrt(rchisq(1,length(Z1)-rank.O2))
          bbb  <-  b/bb	 #
          aa <- rnorm(rank.O2)
          Ztemp[[kk]][Z00] <- O2%*%aa+bb*tau
          vert <- c(1:Dim)[-(fe+kk)]
          for(jj in 1:VC[i]){
            VTtemp[vert,jj] <- VTtemp[vert,jj]- VTtemp[fe+kk,jj]*O1%*%(bbb*aa-a) 
            VTtemp[fe+kk,jj] <- VTtemp[fe+kk,jj]*b/bb 
          }
        }else{
          b <- sqrt(crossprod(Z1))[1]
          tau <- Z1/b
          bb <- sqrt(rchisq(1,length(Z1)))
          Ztemp[[kk]][Z00] <- bb*tau
          vert <- c(1:Dim)[-(fe+kk)]
          for(jj in 1:VC[i]){
            VTtemp[fe+kk,jj] <- VTtemp[fe+kk,jj]*b/bb
          }              
        }
      }    
      for(ii in 1:re){
        ZZ[[ii]] <- cbind(ZZ[[ii]],Ztemp[[ii]]) # cbind ou c ?
      }
      VTVT[[i]] <- VTtemp  
    }
    Z <- ZZ 
    #VT <- VTVT
    for(j in 1:N){
      assign(VTnames[j], VTVT[[j]])
    }
    #----------------------------------------------------flip negatives			
    for(i in 1:N){
      VTtemp.i <- get(VTnames[i])
      for(j in 1:re){ #only look at random effects
        if(signs[j,i]==-1){
          VTtemp.i[fe+j,] <- -get(VTnames[i])[fe+j,]  #only need to flip the negatives
          Z[[j]][,i] <- -Z[[j]][,i]             
        }
      }
      assign(VTnames[i], VTtemp.i)
    }    
    if(k_n == K_n){ #if finished pick coordinates			
      #pick the coordinates
      VT_end <- matrix(0, nrow=Dim, ncol=N) 
      for(i in 1:N){
        for(j in 1:Dim){
          if(runif(1)<=0.5){
            if(j<=fe){
              VT_end[j,i] <- min(get(VTnames[i])[j,])
            }else{
              VT_end[j,i] <- max(min(get(VTnames[i])[j,]),0) 
            }
          }else{
            if(j<=fe){
              VT_end[j,i] <- max(get(VTnames[i])[j,])
            }else{
              VT_end[j,i] <- max(max(get(VTnames[i])[j,]),0) 
            }
          }
        }
      }
      VERTEX <- VT_end 
      WEIGHT <- WT 
      p <- fe 
      r <- re 
    }
  }
  return(list(VERTEX=VERTEX,WEIGHT=WEIGHT))
}


