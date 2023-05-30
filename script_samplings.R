rm(list=ls()); graphics.off()
reps <- 1000
library(lfmm);library(LEA);library(qvalue);library(vcfR);library(gplots);library(raster);library(fuzzySim);library(maps);library(BEDASSLE);library(qvalue);library(adegenet); library(gradientForest);library(magrittr);library(plot.matrix);library(MetBrewer)
source("~/Desktop/proyecto_bosques_niebla/muestreos/functions/functions_samplings.R")

setwd("/Users/jonasaaguirre/Desktop/proyecto_bosques_niebla/muestreos/analyses/populations/species/fagus/")
if(FALSE){
  load("all_SNPs/input_sp.R")
  frec <- input_sp$frec
  p <- frec[,seq(1,ncol(frec)-1,2)]
  q <- frec[,seq(2,ncol(frec),2)]
  tot <- p+q
  frec <- p/tot
  bios <- input_sp$env
  bios <- bios[,grep("bio",names(bios))]
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
  fdr <- qvalue(pv.lfmm,fdr.level = 0.1)
  cand <- colnames(temp_frec)[which(fdr$significant)]
  ref <- colnames(temp_frec)[which(!fdr$significant)]
  
  frec <- input_sp$frec
  loc_cand <- NULL
  for(i in 1:length(cand)){
    print(length(cand)-i)
    tp <- grep(cand[i],colnames(frec))
    if(length(tp)>2){
      print(i)
      next
    }
    loc_cand <- c(loc_cand,tp)
  }
  ref <- seq(1,ncol(frec)-1,2)
  ref <- sample(ref,size = 100000,replace = F)
  ref <- c(ref,ref+1)
  ref <- ref[order(ref)]
  loc <- c(loc_cand,ref)
  loc <- loc[order(loc)]
  frec <- frec[,loc]
  input_sp$frec <- frec
  save(input_sp,file = "inputs/input_sp.R")
  
}

run_analyes <- function(reps,species){
  setwd(paste("~/Desktop/proyecto_bosques_niebla/muestreos/analyses/populations/species/",species,sep = ""))
    
  # create stats
  if(TRUE){
    create_input(input = "inputs/input_sp.R",stacks = "~/Desktop/proyecto_bosques_niebla/capas_climaticas/",N_loc = 7000)
  }
  load("results/stats.R")

  frec <- output$frec
  stats <- output$stats
  lay <- output$lay
  offsets <- output$offsets %>% .[grep(lay,names(.))]
  if(TRUE){
    pdf("results/summaries.pdf")
    
    if(TRUE){
      # run random samples
      print("all subsamples")
      
      
      if(nrow(frec)<25){
        Npop <- c(5,10,15,20)
      }
      if(nrow(frec)>=25 & nrow(frec)<30){
        Npop <- c(5,10,15,20,25)
      }
      if(nrow(frec)>=30 ){
        Npop <- c(5,10,15,20,25,30)
      }
      
      pops <- vector("list",length(Npop))
      for(p in 1:length(Npop)){
        temp<- populations_subsample(frec=frec,fut=lay,Npop=Npop[p],stats=stats,reps=reps,plot.pdf = FALSE,names=paste("pop",pops[p],sep = ""),offsets = offsets)
        pops[[p]] <- temp
        names(pops)[p] <- paste("pop_",Npop[p],sep = "")
      }
      all<- populations_subsample(frec=frec,fut=lay,Npop="all",stats=stats,reps=2,plot.pdf = FALSE,names="all",offsets = offsets)
      pops$all <- all
      save(pops,file="results/pops.R")
      
      
    }
    
    load("results/pops.R")
    load("results/stats.R")
    fut <- output$lay
    cols <- colorRampPalette(c("lightgrey","lightblue","blue"))(length(pops))
    plot.boxplots(dataset = pops,outlier = F,cols = cols )
    
    summaries_pops <- plot.summaries(dataset = pops,cols=cols,struct = TRUE,fut=fut)
    save(summaries_pops,file = "results/summaries.R")
    
    
    load("results/stats.R")
    offsets <- output$offsets %>% .[,grep("offset",names(.))]
    boxplot(t(offsets),outline=F,las=2,ylab="Offsets",cex.axis=0.5,main="climatic projections")
    
    temp <- offsets
    cvs <- data.frame(matrix(NA,nrow = nrow(temp),2))
    rownames(cvs) <- rownames(temp)
    names(cvs) <- c("temp","proj")
    CV <-function(x){
      cv <- sd(x,na.rm = T)/abs(mean(x,na.rm = T))*100
      return(cv)
    }
    for(j in 1:nrow(cvs)){
      cv <- CV(as.vector(t(temp[j,])))
      cvs[j,2] <- cv
    }
    barplot(cvs$proj,ylab="CV",las=2,names=rownames(cvs),cex.names=0.7,main="climatic projections")
    
    
    dev.off()
    
  }
  
  
  
  
  if(TRUE){
    # run random samples
    load("results/stats.R")
    frec <- output$frec
    stats <- output$stats
    N <- 100
    Npop <- c(4,5,6,7,8,9,10,15,20)
    pops <- vector(mode = "list",length(Npop))
    for(i in 1:length(Npop)){
      temp <- min_pops(frec=frec,Npop = Npop[i],stats = stats,reps = N)
      pops[[i]] <- temp
      names(pops)[i] <- paste("pop",Npop[i],sep = "")
    }
  
    save(pops,file="results/pops5_10.R") 
  }
  
}


run_analyes_large <- function(reps,species){
  setwd(paste("~/Desktop/proyecto_bosques_niebla/muestreos/analyses/populations/species/",species,sep = ""))
  
  # create stats
  if(TRUE){
    create_input(input = "inputs/input_sp.R",stacks = "~/Desktop/proyecto_bosques_niebla/capas_climaticas/",N_loc = 15000,simplify = T,N_loc_sim = 200000)
  }
  load("results/stats.R")
  
  frec <- output$frec
  stats <- output$stats
  lay <- output$lay
  offsets <- output$offsets %>% .[grep(lay,names(.))]
  if(TRUE){
    pdf("results/summaries.pdf")
    
    if(TRUE){
      # run random samples
      print("all subsamples")
      
      
      if(nrow(frec)<25){
        Npop <- c(5,10,15,20)
      }
      if(nrow(frec)>=25 & nrow(frec)<30){
        Npop <- c(5,10,15,20,25)
      }
      if(nrow(frec)>=30 ){
        Npop <- c(5,10,15,20,25,30)
      }
      
      pops <- vector("list",length(Npop))
      for(p in 1:length(Npop)){
        temp<- populations_subsample(frec=frec,fut=lay,Npop=Npop[p],stats=stats,reps=reps,plot.pdf = FALSE,names=paste("pop",pops[p],sep = ""),offsets = offsets)
        pops[[p]] <- temp
        names(pops)[p] <- paste("pop_",Npop[p],sep = "")
      }
      all<- populations_subsample(frec=frec,fut=lay,Npop="all",stats=stats,reps=2,plot.pdf = FALSE,names="all",offsets = offsets)
      pops$all <- all
      save(pops,file="results/pops.R")
      
      
    }
    
    load("results/pops.R")
    load("results/stats.R")
    fut <- output$lay
    cols <- colorRampPalette(c("lightgrey","lightblue","blue"))(length(pops))
    plot.boxplots(dataset = pops,outlier = F,cols = cols )
    
    summaries_pops <- plot.summaries(dataset = pops,cols=cols,struct = TRUE,fut=fut)
    save(summaries_pops,file = "results/summaries.R")
    
    
    load("results/stats.R")
    offsets <- output$offsets %>% .[,grep("offset",names(.))]
    boxplot(t(offsets),outline=F,las=2,ylab="Offsets",cex.axis=0.5,main="climatic projections")
    
    temp <- offsets
    cvs <- data.frame(matrix(NA,nrow = nrow(temp),2))
    rownames(cvs) <- rownames(temp)
    names(cvs) <- c("temp","proj")
    CV <-function(x){
      cv <- sd(x,na.rm = T)/abs(mean(x,na.rm = T))*100
      return(cv)
    }
    for(j in 1:nrow(cvs)){
      cv <- CV(as.vector(t(temp[j,])))
      cvs[j,2] <- cv
    }
    barplot(cvs$proj,ylab="CV",las=2,names=rownames(cvs),cex.names=0.7,main="climatic projections")
    
    
    dev.off()
    
  }
  
  
  if(TRUE){
    # run random samples
    load("results/stats.R")
    frec <- output$frec
    stats <- output$stats
    N <- 100
    Npop <- c(4,5,6,7,8,9,10,15,20)
    pops <- vector(mode = "list",length(Npop))
    for(i in 1:length(Npop)){
      temp <- min_pops(frec=frec,Npop = Npop[i],stats = stats,reps = N)
      pops[[i]] <- temp
      names(pops)[i] <- paste("pop",Npop[i],sep = "")
    }
    
    save(pops,file="results/pops5_10.R") 
  }
  
}


if(TRUE){
  #run_analyes(reps = 1000,species = "mexicana")
  #run_analyes(reps = 1000,species = "parviglumis")
  #run_analyes_large(reps = 100,species = "lyrata")
  #run_analyes(reps = 1000,species = "fagus")
  #run_analyes(reps=1000,species = "empidonax")
  run_analyes(reps = 1000,species = "lyrata")
}

if(FALSE){
  #source("~/Desktop/proyecto_bosques_niebla/muestreos/functions/figures_final.R")
  rm(list=ls()); graphics.off()
  
  source("~/Desktop/proyecto_bosques_niebla/muestreos/functions/functions_samplings.R")
  
  
  #species <- c("fagus","lyrata","parviglumis","mexicana")
  species <- c("fagus")
  for(sp in species){
    print(sp)
    
    load(paste(sp,"/inputs/input_sp.R",sep = ""))
    frec <- input_sp$frec
    load(paste(sp,"/results/pops.R",sep = ""))
    smp <- names(pops)
    for(i in 1:length(smp)){
      
      print(smp[i])
      output <- pops[[smp[i]]]$output
      if(sp=="lyrata"){
        loc <- colnames(frec)[seq(1,ncol(frec)-1,2)]
        loc <- sample(1:length(loc),size = 20000,replace = F) %>% .[order(.)]
        
        print(loc)
        temp <- calculate_FST(frec = frec,output = output,loc = loc,rdmloc = T)
        
      }else{
        temp <- calculate_FST(frec = frec,output = output,loc = NA,rdmloc = F)
      }
      
      pops[[smp[i]]]$fst <- temp
      
    }
    
    
    save(pops,file = paste(sp,"/results/pops.R",sep = ""))
    
    
    
    
    
  }
  
  
}
