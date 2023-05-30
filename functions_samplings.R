## functions used to run the analyses

#load functions
library(lfmm);library(LEA);library(qvalue);library(qqman);
library(pcadapt);library(vcfR);library(gplots);library(raster);
library(fuzzySim);library(maps);library(BEDASSLE);library(qvalue);library(adegenet)
library(vegan);library(psych);library(geosphere); library(gradientForest);
library(magrittr);library(MaizePal);library(spacemakeR);library(hierfstat)





### functions individuals

create_input_ind <- function(input){
  load(input)
  #calculate hs
  hs <- basic.stats(data = input_genind)
  hs_pops <- colMeans(hs$Hs,na.rm = T)
  hs <- hs$overall["Hs"]
  print(paste("hs=",hs))
  
  frec <- input_genind$tab[,seq(1,ncol(input_genind$tab)-1,2)]
  colnames(frec) <- sapply(strsplit(colnames(frec),split = ".",fixed = T),"[",1)
  
  dir.create("geno")
  write.geno(frec,output.file = "geno/base.geno")
  write.lfmm(frec,output.file = "geno/base.lfmm")
  setwd("geno/")
  project =NULL
  project = snmf("base.geno",K = 1:10,entropy = TRUE,project="new")
  plot(project, col = "blue", pch = 19, cex = 1.2)
  print("select K")
  K <- scan(what=character(0),nlines=1,quiet=T)
  K <- as.numeric(K)
  project.missing = snmf("base.lfmm", K = K,
                         entropy = TRUE, repetitions = 1,
                         project = "new")
  impute(project.missing, "base.lfmm",
         method = "mode", K = K)
  input <- read.csv("base.lfmm_imputed.lfmm",sep = " ",header = F)
  rownames(input) <- rownames(frec)
  colnames(input) <- colnames(frec)
  input <- data.frame(pop=input_genind$pop,input)
  
  others <- input_genind$other
  
  if(length(grep("asc",names(others)))>0){
    names(others) <- sub("bio","bio_",sub(".asc","",names(others)))
  }
  setwd("..")
  if(names(others)[1]==c("longitude")){
    names(others)[1:2] <- c("X","Y")
  }
  others <- others[,c("X","Y",paste("bio_",1:19,sep = ""))]
  input_sp <- list(frec=input,env=others)
  save(input_sp,file = "input_sp.R")
  
} 


individual_sampling <- function(input,Ntot=c(1:10),nrep=10,N_pop=20){
  load(input)
  frec <- input_sp$frec; frec$pop <- as.vector(frec$pop)
  clim <- input_sp$env
  names(clim) <- sub(pattern = "bio",replacement = "pres.bio",names(clim))
  ras <- data.frame(raster::stack("~/Desktop/proyecto_bosques_niebla/capas_climaticas/wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp370_2061-2080.tif") %>% raster::extract(.,clim[,c("X","Y")]))
  names(ras)<-paste("bio_",1:19,sep = "")
  names(ras)<-paste("fut.",names(ras),sep = "")
  clim <- cbind(clim,ras)
  
  x <- prcomp(clim[,grep("pres",names(clim))],center = T,scale. =T )$x[,1]
  struct <- frec[,-1] 
  struct <- prcomp(struct)
  sdev <- struct$sdev[1:20]^2
  sdev <- data.frame(as.matrix(dist(sdev,diag = T,upper = T)))
  sdev <- c(sdev[1,2],sdev[2,3],sdev[3,4],sdev[4,5],
            sdev[5,6],sdev[6,7],sdev[7,8],sdev[8,9],
            sdev[9,10],sdev[10,11],sdev[11,12],sdev[12,13],
            sdev[13,14],sdev[14,15],sdev[15,16],sdev[16,17],
            sdev[17,18],sdev[18,19],sdev[19,20])
  val <- 0
  i <- 1
  while(val<=0){
    val <- sdev[i+1]-sdev[i]
    #print(val)
    i <- i+1
    if(is.na(val)){
      list_outputs[[r]]<-"did not work; no PCA value"
      names(list_outputs)[r] <- paste("rep",r,sep = "_")
      next
    }
  }
  k <- i 
  pc <- as.matrix(x)
  fun <- function(x){
    sum(x)/(2*length(x))
  }
  temp_frec <- data.frame(aggregate(x = frec[,-1],by=list(pop=frec$pop),FUN=fun),row.names = "pop")
  factors.lfmm = lfmm_ridge( Y = temp_frec, X =pc, K = k)
  mod.lfmm <- lfmm_test( Y = as.matrix(temp_frec), X =pc, lfmm = factors.lfmm, calibrate = "gif")
  pv.lfmm = mod.lfmm$calibrated.pvalue
  pv.lfmm[which(is.na(pv.lfmm))]<-1
  pv.lfmm <- as.vector(pv.lfmm)
  cand <- order(pv.lfmm,decreasing = F)[1:100]
  cand <- temp_frec[,cand]
  colnames(cand) <- paste("chr_",colnames(cand),sep = "")
  
  output <- data.frame(clim[,c("X","Y")])
  gf_data <- data.frame(clim[,grep("pres",names(clim))],cand)
  names(gf_data) <- sub("pres.","",names(gf_data))
  maxLevel <- nrow(gf_data)
  gf_model <- gradientForest(gf_data, predictor.vars=colnames(gf_data[,grep("bio",names(gf_data))]),
                             response.vars=colnames(gf_data[,grep("chr",names(gf_data))]), ntree=500, 
                             maxLevel=maxLevel, trace=F, corr.threshold=0.50)
  temp_pres <- clim[,grep("pres",names(clim))]
  temp_fut <- clim[,grep("fut",names(clim))]
  names(temp_pres) <- names(temp_fut) <- sub("pres.","",names(temp_pres))
  
  
  turn_pres <- predict(gf_model,temp_pres)
  turn_fut<- predict(gf_model,temp_fut)
  offset <- euclidian_distance(turn_pres,turn_fut)
  output$all <- offset
  
  # pops i -> N
  pops <- summary(factor(frec$pop))
  pops <- pops[order(pops,decreasing = T)]
  pops <- pops[1:20] %>% names()
  rmv <- gf_data[pops,grep("bio",names(gf_data))]
  rmv <- apply(rmv,MARGIN = 2,FUN = sd)
  rmv <- which(rmv==0)%>%names()
  if(length(rmv)>0){
    temp_pres <- temp_pres[-grep(rmv,names(temp_pres))]
    temp_fut <- temp_fut[names(temp_pres)]
  }
  for(n in Ntot){
    tryCatch(
      {
        print(paste("N=",n))
        for(r in 1:nrep){
          print(paste("missing",nrep-r,"reps"))
          temp_inds <- NULL
          for(i in 1:length(pops)){
            temp <- which(frec$pop==pops[i])
            if(length(temp)< n){
              t_ind <- temp
            }else{
              t_ind <- sample(temp,size = n,replace = F)
            }
            t_ind <- t_ind[order(t_ind)]
            temp_inds <- c(temp_inds,t_ind)
          }
          temp_frec <- frec[temp_inds,c("pop",sub("chr.","",names(cand)))]
          temp_frec <- data.frame(aggregate(x = temp_frec[,-1],by=list(pop=temp_frec$pop),FUN=fun),row.names = "pop")
          names(temp_frec)<-paste("chr.",names(temp_frec),sep = "")
          gf_data <- data.frame(temp_pres[rownames(temp_frec),],temp_frec)
          maxLevel <- nrow(gf_data)
          gf_model <- gradientForest(gf_data, predictor.vars=colnames(gf_data[,grep("bio",names(gf_data))]),
                                     response.vars=colnames(gf_data[,grep("chr",names(gf_data))]), ntree=500, 
                                     maxLevel=maxLevel, trace=F, corr.threshold=0.50)
          if(is.null(gf_model)){
            output$temp <- NA
            names(output)[ncol(output)]<-paste("N_",n,"_rep_",r,sep = "")
            next
          }
          turn_pres <- predict(gf_model,temp_pres[rownames(temp_frec),])
          turn_fut<- predict(gf_model,temp_fut[rownames(temp_frec),])
          offset <- euclidian_distance(turn_pres,turn_fut)
          output$temp <- NA
          output[rownames(temp_frec),"temp"] <- offset
          names(output)[ncol(output)]<-paste("N_",n,"_rep_",r,sep = "")
        }
      },
      error=function(e){
        an.error.occured <- TRUE
      }
    )
    
  }
  
  
  
  ## +1 ind pop
  
  
  # pops i -> N
  pops <- summary(factor(frec$pop))
  pops <- pops[order(pops,decreasing = T)]
  pops.1 <- pops[1:20] %>% names()
  pops.2 <- pops[21:length(pops)] %>% names()
  if(length(pops.2)>5){
    pops.2 <- pops.2[1:5]
  }
  for(n in 1:length(pops.2)){
    print(paste("N+=",n))
    for(r in 1:nrep){
      print(paste("missing",nrep-r,"reps"))
      temp_inds <- NULL
      for(i in 1:length(pops.1)){
        t_ind <- which(frec$pop==pops.1[i])
        t_ind <- t_ind[order(t_ind)]
        temp_inds <- c(temp_inds,t_ind)
      }
      pops.2.temp <- sample(pops.2,size = n,replace = F)
      for(i in 1:length(pops.2.temp)){
        temp <- which(frec$pop==pops.2.temp[i])
        t_ind <- sample(temp,1)
        temp_inds <- c(temp_inds,t_ind)
        
      }
      temp_inds <- temp_inds[order(temp_inds)]
      temp_frec <- frec[temp_inds,c("pop",sub("chr.","",names(cand)))]
      temp_frec <- data.frame(aggregate(x = temp_frec[,-1],by=list(pop=temp_frec$pop),FUN=fun),row.names = "pop")
      names(temp_frec)<-paste("chr.",names(temp_frec),sep = "")
      gf_data <- data.frame(temp_pres[rownames(temp_frec),],temp_frec)
      maxLevel <- nrow(gf_data)
      gf_model <- gradientForest(gf_data, predictor.vars=colnames(gf_data[,grep("bio",names(gf_data))]),
                                 response.vars=colnames(gf_data[,grep("chr",names(gf_data))]), ntree=500, 
                                 maxLevel=maxLevel, trace=F, corr.threshold=0.50)
      if(is.null(gf_model)){
        output$temp <- NA
        names(output)[ncol(output)]<-paste("N_plus_",n,"_rep",r,sep = "")
        next
      }
      turn_pres <- predict(gf_model,temp_pres[rownames(temp_frec),])
      turn_fut<- predict(gf_model,temp_fut[rownames(temp_frec),])
      offset <- euclidian_distance(turn_pres,turn_fut)
      output$temp <- NA
      output[rownames(temp_frec),"temp"] <- offset
      names(output)[ncol(output)]<-paste("N_plus_",n,"_rep",r,sep = "")
    }
    
  } 
  
  
  
  
  
  
  
  return(output)
  
}



individual_sampling_design <- function(input,Npop=1,Nind=1,nrep=10){
  load(input)
  frec <- input_sp$frec
  clim <- input_sp$env
  names(clim) <- sub(pattern = "bio",replacement = "pres.bio",names(clim))
  ras <- data.frame(raster::stack("~/Desktop/proyecto_bosques_niebla/capas_climaticas/wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp370_2061-2080.tif") %>% raster::extract(.,clim[,c("X","Y")]))
  names(ras)<-paste("bio_",1:19,sep = "")
  names(ras)<-paste("fut.",names(ras),sep = "")
  clim <- cbind(clim,ras)
  
  x <- prcomp(clim[,grep("pres",names(clim))],center = T,scale. =T )$x[,1]
  struct <- frec[,-1] 
  struct <- prcomp(struct)
  sdev <- struct$sdev[1:20]^2
  sdev <- data.frame(as.matrix(dist(sdev,diag = T,upper = T)))
  sdev <- c(sdev[1,2],sdev[2,3],sdev[3,4],sdev[4,5],
            sdev[5,6],sdev[6,7],sdev[7,8],sdev[8,9],
            sdev[9,10],sdev[10,11],sdev[11,12],sdev[12,13],
            sdev[13,14],sdev[14,15],sdev[15,16],sdev[16,17],
            sdev[17,18],sdev[18,19],sdev[19,20])
  val <- 0
  i <- 1
  while(val<=0){
    val <- sdev[i+1]-sdev[i]
    #print(val)
    i <- i+1
    if(is.na(val)){
      list_outputs[[r]]<-"did not work; no PCA value"
      names(list_outputs)[r] <- paste("rep",r,sep = "_")
      next
    }
  }
  k <- i 
  pc <- as.matrix(x)
  fun <- function(x){
    sum(x)/(2*length(x))
  }
  temp_frec <- data.frame(aggregate(x = frec[,-1],by=list(pop=frec$pop),FUN=fun),row.names = "pop")
  factors.lfmm = lfmm_ridge( Y = temp_frec, X =pc, K = k)
  mod.lfmm <- lfmm_test( Y = as.matrix(temp_frec), X =pc, lfmm = factors.lfmm, calibrate = "gif")
  pv.lfmm = mod.lfmm$calibrated.pvalue
  pv.lfmm[which(is.na(pv.lfmm))]<-1
  pv.lfmm <- as.vector(pv.lfmm)
  cand <- order(pv.lfmm,decreasing = F)[1:100]
  cand <- temp_frec[,cand]
  colnames(cand) <- paste("chr_",colnames(cand),sep = "")
  
  output <- data.frame(clim[,c("X","Y")])
  gf_data <- data.frame(clim[,grep("pres",names(clim))],cand)
  names(gf_data) <- sub("pres.","",names(gf_data))
  maxLevel <- nrow(gf_data)
  gf_model <- gradientForest(gf_data, predictor.vars=colnames(gf_data[,grep("bio",names(gf_data))]),
                             response.vars=colnames(gf_data[,grep("chr",names(gf_data))]), ntree=500, 
                             maxLevel=maxLevel, trace=F, corr.threshold=0.50)
  temp_pres <- clim[,grep("pres",names(clim))]
  temp_fut <- clim[,grep("fut",names(clim))]
  names(temp_pres) <- names(temp_fut) <- sub("pres.","",names(temp_pres))
  
  
  turn_pres <- predict(gf_model,temp_pres)
  turn_fut<- predict(gf_model,temp_fut)
  offset <- euclidian_distance(turn_pres,turn_fut)
  output$all <- offset
  pops <- unique(frec$pop)
  
  if(TRUE){
    
    for(r in 1:nrep){
      print(paste("missing",nrep-r,"reps"))
      
      temp_inds <- NULL
      pops_1 <- sample(pops,size = Npop,replace = F)
      pops_N <- setdiff(pops,pops_1)
      
      for(i in 1:length(pops_1)){
        temp <- which(frec$pop==pops_1[i])
        t_ind <- sample(temp,size = length(temp)-Nind,replace = F)
        temp_inds <- c(temp_inds,t_ind)
      }
      
      temp_frec <- frec[-temp_inds,c("pop",sub("chr.","",names(cand)))]
      summary(temp_frec$pop) %>% as.vector() %>% .[order(.)] %>% print()
      
      temp_frec <- data.frame(aggregate(x = temp_frec[,-1],by=list(pop=temp_frec$pop),FUN=fun),row.names = "pop")
      names(temp_frec)<-paste("chr.",names(temp_frec),sep = "")
      gf_data <- data.frame(temp_pres,temp_frec)
      maxLevel <- nrow(gf_data)
      gf_model <- gradientForest(gf_data, predictor.vars=colnames(gf_data[,grep("bio",names(gf_data))]),
                                 response.vars=colnames(gf_data[,grep("chr",names(gf_data))]), ntree=500, 
                                 maxLevel=maxLevel, trace=F, corr.threshold=0.50)
      if(is.null(gf_model)){
        output$temp <- NA
        names(output)[ncol(output)]<-paste("N_",n,"_rep_",r,sep = "")
        next
      }
      turn_pres <- predict(gf_model,temp_pres)
      turn_fut<- predict(gf_model,temp_fut)
      offset <- euclidian_distance(turn_pres,turn_fut)
      output$temp <- offset
      names(output)[ncol(output)]<-paste("Npop_m",Npop,"_Nind_m",Nind,"_rep_",r,sep = "")
    }
  }
  
  return(output)
  
}


individual_sampling_design_Ninds <- function(input,Npop=1,Nind=1,nrep=10,max_N=5){
  load(input)
  frec <- input_sp$frec
  clim <- input_sp$env
  names(clim) <- sub(pattern = "bio",replacement = "pres.bio",names(clim))
  ras <- data.frame(raster::stack("~/Desktop/proyecto_bosques_niebla/capas_climaticas/wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp370_2061-2080.tif") %>% raster::extract(.,clim[,c("X","Y")]))
  names(ras)<-paste("bio_",1:19,sep = "")
  names(ras)<-paste("fut.",names(ras),sep = "")
  clim <- cbind(clim,ras)
  
  x <- prcomp(clim[,grep("pres",names(clim))],center = T,scale. =T )$x[,1]
  struct <- frec[,-1] 
  struct <- prcomp(struct)
  sdev <- struct$sdev[1:20]^2
  sdev <- data.frame(as.matrix(dist(sdev,diag = T,upper = T)))
  sdev <- c(sdev[1,2],sdev[2,3],sdev[3,4],sdev[4,5],
            sdev[5,6],sdev[6,7],sdev[7,8],sdev[8,9],
            sdev[9,10],sdev[10,11],sdev[11,12],sdev[12,13],
            sdev[13,14],sdev[14,15],sdev[15,16],sdev[16,17],
            sdev[17,18],sdev[18,19],sdev[19,20])
  val <- 0
  i <- 1
  while(val<=0){
    val <- sdev[i+1]-sdev[i]
    #print(val)
    i <- i+1
    if(is.na(val)){
      list_outputs[[r]]<-"did not work; no PCA value"
      names(list_outputs)[r] <- paste("rep",r,sep = "_")
      next
    }
  }
  k <- i 
  pc <- as.matrix(x)
  fun <- function(x){
    sum(x)/(2*length(x))
  }
  temp_frec <- data.frame(aggregate(x = frec[,-1],by=list(pop=frec$pop),FUN=fun),row.names = "pop")
  factors.lfmm = lfmm_ridge( Y = temp_frec, X =pc, K = k)
  mod.lfmm <- lfmm_test( Y = as.matrix(temp_frec), X =pc, lfmm = factors.lfmm, calibrate = "gif")
  pv.lfmm = mod.lfmm$calibrated.pvalue
  pv.lfmm[which(is.na(pv.lfmm))]<-1
  pv.lfmm <- as.vector(pv.lfmm)
  cand <- order(pv.lfmm,decreasing = F)[1:100]
  cand <- temp_frec[,cand]
  colnames(cand) <- paste("chr_",colnames(cand),sep = "")
  
  output <- data.frame(clim[,c("X","Y")])
  gf_data <- data.frame(clim[,grep("pres",names(clim))],cand)
  names(gf_data) <- sub("pres.","",names(gf_data))
  maxLevel <- nrow(gf_data)
  gf_model <- gradientForest(gf_data, predictor.vars=colnames(gf_data[,grep("bio",names(gf_data))]),
                             response.vars=colnames(gf_data[,grep("chr",names(gf_data))]), ntree=500, 
                             maxLevel=maxLevel, trace=F, corr.threshold=0.50)
  temp_pres <- clim[,grep("pres",names(clim))]
  temp_fut <- clim[,grep("fut",names(clim))]
  names(temp_pres) <- names(temp_fut) <- sub("pres.","",names(temp_pres))
  
  
  turn_pres <- predict(gf_model,temp_pres)
  turn_fut<- predict(gf_model,temp_fut)
  offset <- euclidian_distance(turn_pres,turn_fut)
  output$all <- offset
  pops <- unique(frec$pop)
  
  if(TRUE){
    
    for(r in 1:nrep){
      print(paste("missing",nrep-r,"reps"))
      
      temp_inds <- NULL
      for(i in 1:length(pops)){
        t_ind <- sample(which(frec$pop==pops[i]),size = max_N,replace = F)
        t_ind <- t_ind[order(t_ind)]
        temp_inds <- c(temp_inds,t_ind)
      }
      temp_frec <- frec[temp_inds,c("pop",sub("chr.","",names(cand)))]
      
      pops_1 <- sample(pops,size = Npop,replace = F)
      pops_N <- setdiff(pops,pops_1)
      temp_inds <- NULL
      for(i in 1:length(pops_1)){
        temp <- which(temp_frec$pop==pops_1[i])
        t_ind <- sample(temp,size = length(temp)-Nind,replace = F)
        temp_inds <- c(temp_inds,t_ind)
      }
      
      temp_frec <- temp_frec[-temp_inds,c("pop",sub("chr.","",names(cand)))]
      summary(temp_frec$pop) %>% as.vector() %>% .[order(.)] %>% print()
      
      temp_frec <- data.frame(aggregate(x = temp_frec[,-1],by=list(pop=temp_frec$pop),FUN=fun),row.names = "pop")
      names(temp_frec)<-paste("chr.",names(temp_frec),sep = "")
      gf_data <- data.frame(temp_pres,temp_frec)
      maxLevel <- nrow(gf_data)
      gf_model <- gradientForest(gf_data, predictor.vars=colnames(gf_data[,grep("bio",names(gf_data))]),
                                 response.vars=colnames(gf_data[,grep("chr",names(gf_data))]), ntree=500, 
                                 maxLevel=maxLevel, trace=F, corr.threshold=0.50)
      if(is.null(gf_model)){
        output$temp <- NA
        names(output)[ncol(output)]<-paste("N_",n,"_rep_",r,sep = "")
        next
      }
      turn_pres <- predict(gf_model,temp_pres)
      turn_fut<- predict(gf_model,temp_fut)
      offset <- euclidian_distance(turn_pres,turn_fut)
      output$temp <- offset
      names(output)[ncol(output)]<-paste("Npop_",Npop,"_Nind_",Nind,"_maxN_",max_N,"_rep_",r,sep = "")
    }
  }
  
  return(output)
  
}

calculate_FST <- function(frec,output,loc,rdmloc=T){
  
  reps <- grep("offset",names(output)) %>% names(output)[.]
  mat_fst <- data.frame(matrix(NA,ncol = length(reps),nrow(output)))
  names(mat_fst) <- reps
  rownames(mat_fst) <- rownames(output)
  reps <- colMeans(output[,reps],na.rm = T)
  reps <- reps[!is.na(reps)] %>% names()
  p <- frec[,seq(1,ncol(frec)-1,2)]
  q <- frec[,seq(2,ncol(frec),2)]
  
  if(rdmloc){
    p <- p[,loc]
    q <- q[,loc]
  }
  
  for(i in 1:length(reps)){
    print(length(reps)-i)
    pops_t <- output[reps[i]]
    pops_t <- rownames(output)[which(!is.na(pops_t[,1]))]
    pt <- p[pops_t,]
    qt <- q[pops_t,]
    fst <- calculate.all.pairwise.Fst(allele.counts = pt,sample.sizes = pt+qt) %>% data.frame()
    names(fst) <- rownames(fst) <- pops_t
    fst[fst==0] <- NA
    fst <- colMeans(fst,na.rm = T)
    mat_fst[names(fst),reps[i]] <- fst
  }
  return(mat_fst)
}



#### function create input

create_input <- function(input="inputs/input_sp.R",stacks="~/Desktop/proyecto_bosques_niebla/capas_climaticas/",N_loc=1000,simplify=F,N_loc_sim=50000){
  load(input)
  frec <- input_sp$frec
  env <- input_sp$env
  
  pallele <- frec[,seq(1,ncol(frec)-1,2)]
  qallele <- frec[,seq(2,ncol(frec),2)]
  tot <- pallele+qallele
  frec <- pallele/tot
  na_loc <- colMeans(frec)
  if(length(which(is.na(na_loc)))>0){
    frec <- frec[,-which(is.na(na_loc))]
  }
  
  colnames(frec) <- sapply(strsplit(colnames(frec),".",fixed = T),"[",1)
  
  if(simplify==T){
    smp <- sample(ncol(tot),size = N_loc_sim,replace = F)
    smp <- smp[order(smp)]
    frec <- frec[,smp]
    tot <- tot[,smp]
  }
  
  hs <- 2*frec*(1-frec)
  hs <- rowMeans(hs,na.rm = T)

  if(!is.character(N_loc)){
    smp <- sample(ncol(tot),size = N_loc,replace = F)
    smp <- smp[order(smp)]
  }else{
    smp <- 1:ncol(tot)
  }
  
  # calculate fst
  fst <- calculate.all.pairwise.Fst(allele.counts = pallele[,smp],sample.sizes = tot[,smp])
  fst[fst==0]<-NA
  fst <- colMeans(fst,na.rm = T)
  
  stats <- data.frame(pop=rownames(frec),longitude=env$longitude,latitude=env$latitude,hs=hs,fst=fst)
  
  dirs <- list.files(stacks,full.names = T)
  temp <- dirs[grep("present",dirs)]
  present <- stack(temp)%>% raster::extract(.,stats[,c("longitude","latitude")])
  colnames(present) <- sub("present.","bio_",colnames(present))
  present <- present[,paste("bio_",1:19,sep = "")]
  colnames(present) <- paste(colnames(present),".present",sep = "")
  stats <- data.frame(stats,present)
  dirs <- dirs[-grep("present",dirs)]
  for(i in 1:length(dirs)){
    temp <- raster::stack(dirs[i]) %>% raster::extract(.,stats[,c("longitude","latitude")])
    nom <- strsplit(x = dirs[i],"//")[[1]][2]
    nom <- sub("wc2.1_2.5m_bioc_","",nom) %>% sub("tif","",.) %>% sub(".","@",.,fixed = T)
    nom <- sub("wc2_",nom,colnames(temp))
    colnames(temp) <- paste("project",nom,sep="_")
    stats <- cbind(stats,temp)
    print(length(dirs)-i)
    }

  temp <- names(stats)[grep("project",names(stats))]
  projections <- unique(sapply(strsplit(temp,"@",fixed = T),"[",1))
  for(i in projections){
    fut <- stats[,grep(i,names(stats))]
    env <- abs(fut-present)
    env <- rowSums(env)
    stats$temp <- env
    names(stats)[ncol(stats)]<-sub("project","env",i)
  }
  stats$env_mean <- rowMeans(stats[,grep("env_",names(stats))])
  
  centroid_niche <- stats[grep("present",names(stats))]
  pca <- data.frame(prcomp(centroid_niche,center = T,scale. = T)$x[,1:6] )
  niche_centroid <- colMeans(pca[,grep("PC",colnames(pca))])
  stats$niche_centroid <- mapply(FUN = function(x,y)(x-y)^2,x=pca[,grep("PC",colnames(pca))],y=niche_centroid) %>% rowSums()%>% sqrt()
  
  geo_centroid <- colMeans(stats[,c("longitude","latitude")])
  stats$geo_centroid <- mapply(FUN = function(x,y)(x-y)^2,x=stats[,c("longitude","latitude")],y=geo_centroid) %>% rowSums()%>% sqrt()
  
  
  pres <- stats[,grep(".present",names(stats))]
  coords <- stats[,c("longitude","latitude")]
  bios <- pres
  pc <- prcomp(bios,retx = T,center = T,scale. = T)$x[,1]
  
  # obtain K
  struct <- temp_frec <- frec
  struct <- prcomp(struct)
  sdev <- struct$sdev[1:20]^2
  sdev <- data.frame(as.matrix(dist(sdev,diag = T,upper = T)))
  sdev <- c(sdev[1,2],sdev[2,3],sdev[3,4],sdev[4,5],
            sdev[5,6],sdev[6,7],sdev[7,8],sdev[8,9],
            sdev[9,10],sdev[10,11],sdev[11,12],sdev[12,13],
            sdev[13,14],sdev[14,15],sdev[15,16],sdev[16,17],
            sdev[17,18],sdev[18,19],sdev[19,20])
  val <- 0
  i <- 1
  while(val<=0){
    val <- sdev[i+1]-sdev[i]
    #print(val)
    i <- i+1
    if(is.na(val)){
      test <- c(test,"NA_PCA_fail")
      next
    }
  }
  k <- i 
  # calculate lfmm
  pc <- as.matrix(pc)
  factors.lfmm = lfmm_ridge( Y = temp_frec, X =pc, K = k)
  mod.lfmm <- lfmm_test( Y = as.matrix(temp_frec), X =pc, lfmm = factors.lfmm, calibrate = "gif")
  pv.lfmm = mod.lfmm$calibrated.pvalue
  pv.lfmm[which(is.na(pv.lfmm))]<-1
  pv.lfmm <- as.vector(pv.lfmm)
  cand <- order(pv.lfmm,decreasing = F)[1:100]
  cand <- temp_frec[,cand]
  colnames(cand) <- paste("chr_",colnames(cand),sep = "")
  temp_pres <- stats[,grep(".present",names(stats))]
  names(temp_pres)<- sub(".present","",names(temp_pres))
  gf_data <- data.frame(temp_pres,cand)
  maxLevel <- nrow(gf_data)
  gf_model <- gradientForest(gf_data, predictor.vars=colnames(gf_data[,grep("bio",names(gf_data))]),
                             response.vars=colnames(gf_data[,grep("chr",names(gf_data))]), ntree=500, 
                             maxLevel=maxLevel, trace=F, corr.threshold=0.50)
  turn_pres <- predict(gf_model,temp_pres)
  
  offsets <- coords[,c("longitude","latitude")]
  temp <- names(stats)[grep("project",names(stats))]
  projections <- unique(sapply(strsplit(temp,"@",fixed = T),"[",1))
  
  for(i in 1:length(projections)){
    print(i)
    temp_fut <- stats[,grep(pattern = projections[i],names(stats))]
    names(temp_fut) <- paste("bio_",sapply(strsplit(names(temp_fut),split = "@",fixed = T),"[",2),sep = "")
    turn_fut<- predict(gf_model,temp_fut)
    offsets$temp <- euclidian_distance(turn_pres,turn_fut)
    names(offsets)[ncol(offsets)] <- sub("project","offset",projections[i])
  }
  corr <- cor(offsets[,-c(1:2)], method = 'spearman')
  colors = met.brewer(name="Hiroshige")
  hist(as.dist(corr,diag = F,upper = F),ylab="rank correlations",main="")
  corr[corr<0.5]<-NA
  succ <- round((length(which(is.na(corr)))/(length(which(!is.na(corr)))+length(which(is.na(corr)))))*100,2)
  
  plot(corr,col=colors[5:1],las=2,ylab="",xlab="",main=paste("Rank correlation (",succ,"% <0.5)",sep = ""))
  
  cmb <- names(offsets)[-c(1:3)]
  cmb <- combn(x = cmb,m = 2)
  for(i in 1:ncol(cmb)){
    temp <- cor.test(offsets[,cmb[1,i]],offsets[,cmb[2,i]],method = "spearman")
  }
  corr_warm <- corr#[-grep("2011-2040",rownames(corr)),-grep("2011-2040",rownames(corr))]
  meds <- NULL
  for(i in 1:nrow(corr_warm)){
    cr <- corr_warm[i,]
    cr <- cr[-which(cr==1)]
    meds <- c(meds,median(cr,na.rm = T))
    names(meds)[length(meds)]<- rownames(corr_warm)[i]
  }
  lay <- names(meds[order(meds,decreasing = T)][1])
  lay <- sub("offset_","",lay)
  output <- list(stats=stats,frec=frec,offsets=offsets,lay=lay,gf_model=gf_model)
  save(output,file = "results/stats.R")
  
  }


# function all offsets
all_offsets <- function(frec,stats,grid="wc2.1_2.5m_bioc_",year="2041-2060",year2="2050",path="~/Desktop/future_low_wc/"){
  #create output dataset
  pres <- stats[,grep(".present",names(stats))]
  coords <- stats[,c("longitude","latitude")]
  
  #files <- list.files(path = path,pattern = "tif",full.names = T)
  #files <- files[grep(year,files)]
  #l.layers <- NULL
  #for(i in 1:length(files)){
   # print(length(files)-i)
   # layer <- sapply(strsplit(files[i],"//"),"[",2)
   # layer <- sub(".tif","",sub(pattern = year,year2,sub(pattern = grid,"",layer)))
   # l.layers <- c(l.layers,layer)
   # ras <- data.frame(raster::stack(files[i]) %>% raster::extract(.,coords[,c("longitude","latitude")]))
   # names(ras) <- paste("bio_",1:19,".",layer,sep = "")
   # coords <- data.frame(coords,ras)
  #}
  #obtain environmental PC1 
  bios <- pres
  pc <- prcomp(bios,retx = T,center = T,scale. = T)$x[,1]
  
  # obtain K
  struct <- temp_frec <- frec
  struct <- prcomp(struct)
  sdev <- struct$sdev[1:20]^2
  sdev <- data.frame(as.matrix(dist(sdev,diag = T,upper = T)))
  sdev <- c(sdev[1,2],sdev[2,3],sdev[3,4],sdev[4,5],
            sdev[5,6],sdev[6,7],sdev[7,8],sdev[8,9],
            sdev[9,10],sdev[10,11],sdev[11,12],sdev[12,13],
            sdev[13,14],sdev[14,15],sdev[15,16],sdev[16,17],
            sdev[17,18],sdev[18,19],sdev[19,20])
  val <- 0
  i <- 1
  while(val<=0){
    val <- sdev[i+1]-sdev[i]
    #print(val)
    i <- i+1
    if(is.na(val)){
      test <- c(test,"NA_PCA_fail")
      next
    }
  }
  k <- i 
  # calculate lfmm
  pc <- as.matrix(pc)
  factors.lfmm = lfmm_ridge( Y = temp_frec, X =pc, K = k)
  mod.lfmm <- lfmm_test( Y = as.matrix(temp_frec), X =pc, lfmm = factors.lfmm, calibrate = "gif")
  pv.lfmm = mod.lfmm$calibrated.pvalue
  pv.lfmm[which(is.na(pv.lfmm))]<-1
  pv.lfmm <- as.vector(pv.lfmm)
  cand <- order(pv.lfmm,decreasing = F)[1:100]
  cand <- temp_frec[,cand]
  colnames(cand) <- paste("chr_",colnames(cand),sep = "")
  temp_pres <- stats[,grep(".present",names(stats))]
  names(temp_pres)<- sub(".present","",names(temp_pres))
  gf_data <- data.frame(temp_pres,cand)
  maxLevel <- nrow(gf_data)
  gf_model <- gradientForest(gf_data, predictor.vars=colnames(gf_data[,grep("bio",names(gf_data))]),
                             response.vars=colnames(gf_data[,grep("chr",names(gf_data))]), ntree=500, 
                             maxLevel=maxLevel, trace=F, corr.threshold=0.50)
  turn_pres <- predict(gf_model,temp_pres)
  l.layers <- gsub("-",".",l.layers)
  offsets <- coords[,c("longitude","latitude")]
  for(i in 1:length(l.layers)){
    print(i)
    temp_fut <- coords[,grep(pattern = l.layers[i],names(coords))]
    names(temp_fut) <- sapply(strsplit(names(temp_fut),split = ".",fixed = T),"[",1)
    turn_fut<- predict(gf_model,temp_fut)
    offsets$temp <- euclidian_distance(turn_pres,turn_fut)
    names(offsets)[ncol(offsets)] <- l.layers[i]
  }
  return(offsets)
  
}



#### create statistics table ####

# function random populations
populations_subsample <- function(frec,fut="project70_370",Npop="all",stats,reps=10,plot.pdf=TRUE,names="all",offsets){
  
  #create output dataset
  print(Npop)
  output <- data.frame(stats,N_pop=Npop)
  pres <- stats[,grep(".present",names(stats))]
  mean_nich <- mean(stats$niche_centroid)
  mean_geo <- mean(stats$geo_centroid)
  sum_niche <- sum(abs(stats$niche_centroid-mean_nich))/nrow(stats)
  sum_geo <- sum(abs(stats$geo_centroid-mean_geo))/nrow(stats)
  data.f.struct <- data.frame(pops="all",i=0,sum_niche=sum_niche,sd_niche=sd(stats$niche_centroid),sum_geo=sum_geo,sd_geo=sd(stats$geo_centroid))
  df_offsets <- data.frame(pops=rownames(offsets),real=offsets[,1],row.names = "pops")
  
  temp <- data.frame(matrix(data = NA,nrow = nrow(output),ncol = reps*2))
  names(temp)<-c(paste("offset",1:reps,sep = "_"),
                 paste("NOutlier",1:reps,sep = "_"))
  output <- data.frame(output,temp)
  list_outputs <- vector("list",reps)
  list_turnovers <- vector("list",reps)
  
  list_bios <- NULL
  if(plot.pdf){
    pdf(paste("results/",names,"_pca.pdf",sep = ""))
  }
  #run replicates
  for(r in 1:reps){
    tryCatch(
      {
        print(paste("missing",reps-r))
        #obtain environmental PC1 
        if(Npop=="all"){
          pops <- rownames(output)
        }else{
          pops <- rownames(output)[sample(1:nrow(output),size = Npop,replace = F)]
        }
        
        bios <- pres[pops,]
        pc <- prcomp(bios,retx = T,center = T,scale. = T)$x[,1]
        
        temp.struct <- stats[pops,c("niche_centroid","geo_centroid")]
        sum_niche <- sum(abs(temp.struct$niche_centroid-mean_nich))/nrow(temp.struct)
        sum_geo <- sum(abs(temp.struct$geo_centroid-mean_geo))/nrow(temp.struct)
        temp.struct <- data.frame(pops=paste0(pops,collapse = ","),i=r,
                                  sum_niche=sum_niche,
                                  sd_niche=sd(temp.struct$niche_centroid),
                                  sum_geo=sum_geo,
                                  sd_geo=sd(temp.struct$geo_centroid))
        
        data.f.struct <- rbind(data.f.struct,temp.struct)
        
        
        # obtain K
        struct <- temp_frec <- frec[pops,]
        struct <- prcomp(struct)
        sdev <- struct$sdev[1:20]^2
        sdev <- data.frame(as.matrix(dist(sdev,diag = T,upper = T)))
        sdev <- c(sdev[1,2],sdev[2,3],sdev[3,4],sdev[4,5],
                  sdev[5,6],sdev[6,7],sdev[7,8],sdev[8,9],
                  sdev[9,10],sdev[10,11],sdev[11,12],sdev[12,13],
                  sdev[13,14],sdev[14,15],sdev[15,16],sdev[16,17],
                  sdev[17,18],sdev[18,19],sdev[19,20])
        val <- 0
        i <- 1
        while(val<=0){
          val <- sdev[i+1]-sdev[i]
          #print(val)
          i <- i+1
          if(is.na(val)){
            list_outputs[[r]]<-"did not work; no PCA value"
            list_turnovers[[r]]<-"did not work; no PCA value"
            names(list_outputs)[r] <- paste("rep",r,sep = "_")
            names(list_turnovers)[r] <- paste("rep",r,sep = "_")
            next
          }
        }
        k <- i 
        if(plot.pdf){
          par(mfrow=c(1,2))
          plot(struct$x[,1:2],xlab="PC1",ylab="PC1")
          plot(struct$sdev[1:length(pops)]^2,xlab="PC",ylab="Variance Explained")
          points(k,struct$sdev[k]^2,type="h",col="blue")
          par(mfrow=c(1,1))
        }
        
        # calculate lfmm
        pc <- as.matrix(pc)
        factors.lfmm = lfmm_ridge( Y = temp_frec, X =pc, K = k)
        mod.lfmm <- lfmm_test( Y = as.matrix(temp_frec), X =pc, lfmm = factors.lfmm, calibrate = "gif")
        pv.lfmm = mod.lfmm$calibrated.pvalue
        pv.lfmm[which(is.na(pv.lfmm))]<-1
        pv.lfmm <- as.vector(pv.lfmm)
        if(sd(pv.lfmm)==0){
          list_outputs[[r]]<-"did not work; no candidate SNPs"
          list_turnovers[[r]] <- "did not work; no candidate SNPs"
          names(list_outputs)[r] <- paste("rep",r,sep = "_")
          names(list_turnovers)[r] <- paste("rep",r,sep = "_")
          
          next
        }
        cand <- order(pv.lfmm,decreasing = F)[1:100]
        loc <- paste("L",cand,sep = "")
        cand <- temp_frec[,cand]
        colnames(cand) <- paste("chr_",colnames(cand),sep = "")
        temp_pres <- stats[pops,grep(".present",names(stats))];temp_fut <- stats[pops,grep(fut,names(stats))]; temp_fut <- temp_fut[grep("project",names(temp_fut))]
        names(temp_pres)<- sub(".present","",names(temp_pres))
        names(temp_fut)<- paste("bio_",sapply(strsplit(names(temp_fut),"@",fixed = T),"[",2),sep="")
        gf_data <- data.frame(temp_pres,cand)
        maxLevel <- nrow(gf_data)
        gf_model <- gradientForest(gf_data, predictor.vars=colnames(gf_data[,grep("bio",names(gf_data))]),
                                   response.vars=colnames(gf_data[,grep("chr",names(gf_data))]), ntree=500, 
                                   maxLevel=maxLevel, trace=F, corr.threshold=0.50)
        if(is.null(gf_model)){
          list_outputs[[r]]<-"did not work; no gf loc"
          list_turnovers[[r]] <- "did not work; no gf loc"
          names(list_outputs)[r] <- paste("rep",r,sep = "_")
          names(list_turnovers)[r] <- paste("rep",r,sep = "_")
          next
        }
        
        turn_pres <- predict(gf_model,temp_pres)
        list_turnovers[[r]] <- turn_pres
        names(list_turnovers)[r] <-  paste("rep",r,sep = "_")
        turn_fut<- predict(gf_model,temp_fut)
        offset <- euclidian_distance(turn_pres,turn_fut)
        temp_all_pres <- stats[,grep(".present",names(stats))];temp_all_fut <- stats[,grep(fut,names(stats))]; temp_all_fut <- temp_all_fut[grep("project",names(temp_all_fut))]
        names(temp_all_pres) <- sub(".present","",names(temp_all_pres))
        names(temp_all_fut) <- paste("bio_",sapply(strsplit(names(temp_all_fut),"@"),"[",2),sep = "")
        turn_all_pres <- predict(gf_model,temp_all_pres)
        turn_all_fut <- predict(gf_model,temp_all_fut)
        offset_all <- euclidian_distance(turn_all_pres,turn_all_fut)
        
        df_offsets <- cbind(df_offsets,temp=offset_all); names(df_offsets)[ncol(df_offsets)]<-paste("sim_all_",r,sep = "")
        df_offsets$temp <- NA
        temp_off <- data.frame(df_offsets[,1:2],all=offset_all)
        df_offsets[pops,"temp"] <- temp_off[pops,"all"]; names(df_offsets)[ncol(df_offsets)]<-paste("sim_rep_",r,sep = "")
        
        NOutlier <- length(gf_model$result)
        loc_gfout <- gf_model$result
        L_Bios <- gf_model$overall.imp[order(gf_model$overall.imp,decreasing = T)]
        output[pops,paste("offset_",r,sep = "")] <- offset
        output[pops,paste("NOutlier_",r,sep = "")] <- NOutlier
        list_bios <- c(list_bios,names(L_Bios)[1:4])
        cats <- list(Bios=L_Bios,Loc_Lfmm=loc,Pops=pops,Loc_gf=loc_gfout)
        list_outputs[[r]] <- cats
        names(list_outputs)[r] <- paste("rep",r,sep = "_")
      },
      error=function(e){
        an.error.occured <- TRUE
      }
    )
  }
  
  
  if(plot.pdf==TRUE){
    dev.off()
  }
  
  off_means <- rowMeans(output[,grep("offset",names(output))],na.rm = T)
  output$off_means <- NA
  output[names(off_means),"off_means"] <- off_means
  
  off_all <- t(output[,grep("offset",names(output))])
  output$range_off <- NA
  for(i in 1:ncol(off_all)){
    temp <- range(off_all[,i],na.rm = T)
    if(is.infinite(temp)[1]){
      output[colnames(off_all)[i],"range_off"] <- NA
    }else{
      output[colnames(off_all)[i],"range_off"] <- temp[2]-temp[1]
    }
    
  }
  
  final <- list(output=output,list_outputs=list_outputs,list_bios=list_bios,data.f.struct=data.f.struct,df_offsets=df_offsets,list_turnovers=list_turnovers)
  return(final)
  
}


#function min pop

min_pops <- function(frec,Npop=5,stats,reps=10){
  #create output dataset
  print(Npop)
  pres <- stats[,grep(".present",names(stats))]
  test <- NULL
  #run replicates
  for(r in 1:reps){
    tryCatch(
      {
        print(paste("missing",reps-r))
        #obtain environmental PC1 
        pops <- rownames(stats)[sample(1:nrow(stats),size = Npop,replace = F)]
        bios <- pres[pops,]
        pc <- prcomp(bios,retx = T,center = T,scale. = T)$x[,1]

        # obtain K
        struct <- temp_frec <- frec[pops,]
        struct <- prcomp(struct)
        sdev <- struct$sdev[1:20]^2
        sdev <- data.frame(as.matrix(dist(sdev,diag = T,upper = T)))
        sdev <- c(sdev[1,2],sdev[2,3],sdev[3,4],sdev[4,5],
                  sdev[5,6],sdev[6,7],sdev[7,8],sdev[8,9],
                  sdev[9,10],sdev[10,11],sdev[11,12],sdev[12,13],
                  sdev[13,14],sdev[14,15],sdev[15,16],sdev[16,17],
                  sdev[17,18],sdev[18,19],sdev[19,20])
        val <- 0
        i <- 1
        while(val<=0){
          val <- sdev[i+1]-sdev[i]
          #print(val)
          i <- i+1
          if(is.na(val)){
            test <- c(test,"NA_PCA_fail")
            next
          }
        }
        k <- i 
        # calculate lfmm
        pc <- as.matrix(pc)
        factors.lfmm = lfmm_ridge( Y = temp_frec, X =pc, K = k)
        mod.lfmm <- lfmm_test( Y = as.matrix(temp_frec), X =pc, lfmm = factors.lfmm, calibrate = "gif")
        pv.lfmm = mod.lfmm$calibrated.pvalue
        pv.lfmm[which(is.na(pv.lfmm))]<-1
        pv.lfmm <- as.vector(pv.lfmm)
        if(sd(pv.lfmm)==0){
          test <- c(test,"NA_lfmm_fail")
          next
        }
        cand <- order(pv.lfmm,decreasing = F)[1:100]
        cand <- temp_frec[,cand]
        colnames(cand) <- paste("chr_",colnames(cand),sep = "")
        temp_pres <- stats[pops,grep(".present",names(stats))]
        names(temp_pres)<- sub(".present","",names(temp_pres))
        gf_data <- data.frame(temp_pres,cand)
        maxLevel <- nrow(gf_data)
        gf_model <- gradientForest(gf_data, predictor.vars=colnames(gf_data[,grep("bio",names(gf_data))]),
                                   response.vars=colnames(gf_data[,grep("chr",names(gf_data))]), ntree=500, 
                                   maxLevel=maxLevel, trace=F, corr.threshold=0.50)
        if(is.null(gf_model)){
          test <- c(test,"NA_fail_GF")
          next
        }
        test <- c(test,"success")
      },
      error=function(e){
        an.error.occured <- TRUE
      }
    )
  }

  return(test)
  
}



  
  
# function genetic offset
euclidian_distance <- function(proj_fut=future,pred_pres=present){
  num <- ncol(proj_fut)
  tot <- rep(0,nrow(proj_fut))
  for(i in 1:num ){
    sum <- (proj_fut[,i]-pred_pres[,i])^2
    tot <- tot+sum
  }
  tot <- sqrt(tot)
  return(tot)
}

CV <-function(x){
  cv <- sd(x,na.rm = T)/abs(mean(x,na.rm = T))*100
  return(cv)
}
# plot figures


fine_struct <- function(dataset,N=50){
  pops <- dataset
  all <- pops$all$output$off_means
  tests <- names(pops$pops_5_env$output)[grep("offset",names(pops$pops_5_env$output))]
  
  pops_5_env_high <- rowMeans(pops$pops_5_env$output[,head(tests,N)],na.rm = T)
  pops_5_geo_high <- rowMeans(pops$pops_5_geo$output[,head(tests,N)],na.rm = T)
  pops_5_env_low <- rowMeans(pops$pops_5_env$output[,tail(tests,N)],na.rm = T)
  pops_5_geo_low <- rowMeans(pops$pops_5_geo$output[,tail(tests,N)],na.rm = T)
  
  pops_10_env_high <- rowMeans(pops$pops_10_env$output[,head(tests,N)],na.rm = T)
  pops_10_geo_high <- rowMeans(pops$pops_10_geo$output[,head(tests,N)],na.rm = T)
  pops_10_env_low <- rowMeans(pops$pops_10_env$output[,tail(tests,N)],na.rm = T)
  pops_10_geo_low <- rowMeans(pops$pops_10_geo$output[,tail(tests,N)],na.rm = T)
  
  pops_5_env_high <- cor.test(x=all, y=pops_5_env_high, method = 'spearman')
  pops_5_geo_high <- cor.test(x=all, y=pops_5_geo_high, method = 'spearman')
  pops_5_env_low <- cor.test(x=all, y=pops_5_env_low, method = 'spearman')
  pops_5_geo_low <- cor.test(x=all, y=pops_5_geo_low, method = 'spearman')
  
  
  pops_10_env_high <- cor.test(x=all, y=pops_10_env_high, method = 'spearman')
  pops_10_geo_high <- cor.test(x=all, y=pops_10_geo_high, method = 'spearman')
  pops_10_env_low <- cor.test(x=all, y=pops_10_env_low, method = 'spearman')
  pops_10_geo_low <- cor.test(x=all, y=pops_10_geo_low, method = 'spearman')
  
  
  barplot(c(pops_5_env_high$estimate,pops_5_env_low$estimate,
            pops_5_geo_high$estimate,pops_5_geo_low$estimate,
            pops_10_env_high$estimate,pops_10_env_low$estimate,
            pops_10_geo_high$estimate,pops_10_geo_low$estimate),
          names=c("5_env_h","5_env_l","5_geo_h","5_geo_l",
                  "10_env_h","10_env_l","10_geo_h","10_geo_l"),las=2,ylab="Rho")
  
  
  barplot(c(pops_5_env_high$p.value,pops_5_env_low$p.value,
            pops_5_geo_high$p.value,pops_5_geo_low$p.value,
            pops_10_env_high$p.value,pops_10_env_low$p.value,
            pops_10_geo_high$p.value,pops_10_geo_low$p.value),
          names=c("5_env_h","5_env_l","5_geo_h","5_geo_l",
                  "10_env_h","10_env_l","10_geo_h","10_geo_l"),las=2,ylab="p")
  
  
}

fine_struct <- function(dataset,N=50){
  pops <- dataset
  all <- pops$all$output$off_means
  tests <- names(pops$pops_5_env$output)[grep("offset",names(pops$pops_5_env$output))]
  
  pops_5_env_high <- rowMeans(pops$pops_5_env$output[,head(tests,N)],na.rm = T)
  pops_5_geo_high <- rowMeans(pops$pops_5_geo$output[,head(tests,N)],na.rm = T)
  pops_5_env_low <- rowMeans(pops$pops_5_env$output[,tail(tests,N)],na.rm = T)
  pops_5_geo_low <- rowMeans(pops$pops_5_geo$output[,tail(tests,N)],na.rm = T)
  
  pops_10_env_high <- rowMeans(pops$pops_10_env$output[,head(tests,N)],na.rm = T)
  pops_10_geo_high <- rowMeans(pops$pops_10_geo$output[,head(tests,N)],na.rm = T)
  pops_10_env_low <- rowMeans(pops$pops_10_env$output[,tail(tests,N)],na.rm = T)
  pops_10_geo_low <- rowMeans(pops$pops_10_geo$output[,tail(tests,N)],na.rm = T)
  
  pops_5_env_high <- cor.test(x=all, y=pops_5_env_high, method = 'spearman')
  pops_5_geo_high <- cor.test(x=all, y=pops_5_geo_high, method = 'spearman')
  pops_5_env_low <- cor.test(x=all, y=pops_5_env_low, method = 'spearman')
  pops_5_geo_low <- cor.test(x=all, y=pops_5_geo_low, method = 'spearman')
  
  
  pops_10_env_high <- cor.test(x=all, y=pops_10_env_high, method = 'spearman')
  pops_10_geo_high <- cor.test(x=all, y=pops_10_geo_high, method = 'spearman')
  pops_10_env_low <- cor.test(x=all, y=pops_10_env_low, method = 'spearman')
  pops_10_geo_low <- cor.test(x=all, y=pops_10_geo_low, method = 'spearman')
  
  
  barplot(c(pops_5_env_high$estimate,pops_5_env_low$estimate,
            pops_5_geo_high$estimate,pops_5_geo_low$estimate,
            pops_10_env_high$estimate,pops_10_env_low$estimate,
            pops_10_geo_high$estimate,pops_10_geo_low$estimate),
          names=c("5_env_h","5_env_l","5_geo_h","5_geo_l",
                  "10_env_h","10_env_l","10_geo_h","10_geo_l"),las=2,ylab="Rho")
  
  
  barplot(c(pops_5_env_high$p.value,pops_5_env_low$p.value,
            pops_5_geo_high$p.value,pops_5_geo_low$p.value,
            pops_10_env_high$p.value,pops_10_env_low$p.value,
            pops_10_geo_high$p.value,pops_10_geo_low$p.value),
          names=c("5_env_h","5_env_l","5_geo_h","5_geo_l",
                  "10_env_h","10_env_l","10_geo_h","10_geo_l"),las=2,ylab="p")
  
  
}

fine_struct <- function(dataset,N=50){
  pops <- dataset
  all <- pops$all$output$off_means
  tests <- names(pops$pops_5_env$output)[grep("offset",names(pops$pops_5_env$output))]
  
  pops_5_env_high <- rowMeans(pops$pops_5_env$output[,head(tests,N)],na.rm = T)
  pops_5_geo_high <- rowMeans(pops$pops_5_geo$output[,head(tests,N)],na.rm = T)
  pops_5_env_low <- rowMeans(pops$pops_5_env$output[,tail(tests,N)],na.rm = T)
  pops_5_geo_low <- rowMeans(pops$pops_5_geo$output[,tail(tests,N)],na.rm = T)
  
  pops_10_env_high <- rowMeans(pops$pops_10_env$output[,head(tests,N)],na.rm = T)
  pops_10_geo_high <- rowMeans(pops$pops_10_geo$output[,head(tests,N)],na.rm = T)
  pops_10_env_low <- rowMeans(pops$pops_10_env$output[,tail(tests,N)],na.rm = T)
  pops_10_geo_low <- rowMeans(pops$pops_10_geo$output[,tail(tests,N)],na.rm = T)
  
  pops_5_env_high <- cor.test(x=all, y=pops_5_env_high, method = 'spearman')
  pops_5_geo_high <- cor.test(x=all, y=pops_5_geo_high, method = 'spearman')
  pops_5_env_low <- cor.test(x=all, y=pops_5_env_low, method = 'spearman')
  pops_5_geo_low <- cor.test(x=all, y=pops_5_geo_low, method = 'spearman')
  
  
  pops_10_env_high <- cor.test(x=all, y=pops_10_env_high, method = 'spearman')
  pops_10_geo_high <- cor.test(x=all, y=pops_10_geo_high, method = 'spearman')
  pops_10_env_low <- cor.test(x=all, y=pops_10_env_low, method = 'spearman')
  pops_10_geo_low <- cor.test(x=all, y=pops_10_geo_low, method = 'spearman')
  
  
  barplot(c(pops_5_env_high$estimate,pops_5_env_low$estimate,
            pops_5_geo_high$estimate,pops_5_geo_low$estimate,
            pops_10_env_high$estimate,pops_10_env_low$estimate,
            pops_10_geo_high$estimate,pops_10_geo_low$estimate),
          names=c("5_env_h","5_env_l","5_geo_h","5_geo_l",
                  "10_env_h","10_env_l","10_geo_h","10_geo_l"),las=2,ylab="Rho")
  
  
  barplot(c(pops_5_env_high$p.value,pops_5_env_low$p.value,
            pops_5_geo_high$p.value,pops_5_geo_low$p.value,
            pops_10_env_high$p.value,pops_10_env_low$p.value,
            pops_10_geo_high$p.value,pops_10_geo_low$p.value),
          names=c("5_env_h","5_env_l","5_geo_h","5_geo_l",
                  "10_env_h","10_env_l","10_geo_h","10_geo_l"),las=2,ylab="p")
  
  
 tests_rho <-  data.frame(test=c("pop_5_env_high","pop_5_env_low",
                    "pop_5_geo_high","pop_5_geo_low",
                    "pop_10_env_high","pop_10_env_low",
                    "pop_10_geo_high","pop_10_geo_low"),
             rho=c(pops_5_env_high$estimate,pops_5_env_low$estimate,
                   pops_5_geo_high$estimate,pops_5_geo_low$estimate,
                   pops_10_env_high$estimate,pops_10_env_low$estimate,
                   pops_10_geo_high$estimate,pops_10_geo_low$estimate),
             p=c(pops_5_env_high$p.value,pops_5_env_low$p.value,
                 pops_5_geo_high$p.value,pops_5_geo_low$p.value,
                 pops_10_env_high$p.value,pops_10_env_low$p.value,
                 pops_10_geo_high$p.value,pops_10_geo_low$p.value))
 return(tests_rho)
  
}


