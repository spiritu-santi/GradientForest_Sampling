rm(list=ls()); graphics.off()
library(lfmm);library(LEA);library(adegenet);library(qvalue);library(vcfR);library(gplots);library(raster);library(fuzzySim);library(maps);library(BEDASSLE);library(qvalue);library(adegenet); library(gradientForest);library(magrittr);library(plot.matrix);library(MetBrewer);library(ggplot2);library(ggthemes);library(Rmisc)
source("~/Desktop/proyecto_bosques_niebla/muestreos/functions/functions_samplings.R")

# create input inds
if(FALSE){
  setwd("~/Desktop/proyecto_bosques_niebla/muestreos/analyses/individuals/species/mexicana/inputs/")
  create_input_ind(input = "input_genind.R")  
  setwd("~/Desktop/proyecto_bosques_niebla/muestreos/analyses/individuals/species/parviglumis/inputs/")
  create_input_ind(input = "input_genind.R")  
  setwd("~/Desktop/proyecto_bosques_niebla/muestreos/analyses/individuals/species/fagus/inputs/")
  create_input_ind(input = "input_genind.R")  
  setwd("~/Desktop/proyecto_bosques_niebla/muestreos/analyses/individuals/species/empidonax/inputs/")
  create_input_ind(input = "input_genind.R")  
  
}

setwd("~/Desktop/proyecto_bosques_niebla/muestreos/analyses/individuals/species/") 


if(FALSE){
  #species <- list.dirs() %>% strsplit(.,"/") %>% sapply(.,"[",2) %>% unique() %>% .[!is.na(.)]
  species <- "fagus"
  nrep <- 1000
  for(sp in species){
    print(sp)
    temp_input <- paste(sp,"/inputs/input_sp.R",sep = "")
    load(temp_input)
    output_inds <- individual_sampling(input = temp_input,Ntot = c(1:10),N_pop = 20,nrep = nrep)
    save(output_inds,file = paste(sp,"/results/output.R",sep = ""))  
    
  }
    
}


rm(list=ls()); graphics.off()
setwd("~/Desktop/proyecto_bosques_niebla/muestreos/analyses/individuals/species/") 
species <- list.dirs() %>% strsplit(.,"/") %>% sapply(.,"[",2) %>% unique() %>% .[!is.na(.)]
order_na <- NULL
N_loc <- NULL

mat_list <- vector(mode = "list",length = length(species))
names(mat_list) <- species

pdf("results_new.pdf")
data <- data.frame(sp=NA,N_ind=NA,rho=NA)
nombres <- vector(mode = "list",length = length(species))
for(i in 1:length(species)){
  load(paste(species[i],"/inputs/input_genind.R",sep = ""))
  tab_t <- input_genind$tab
  tab_t <- tab_t[,seq(1,ncol(tab_t)-1,2)]
  na <- which(is.na(tab_t)) %>% length(.)
  lna <- which(!is.na(tab_t)) %>% length(.)
  na <- na/(na+lna)
  order_na <- c(order_na,na); names(order_na)[i] <- species[i]
  N_loc <- c(N_loc,length(input_genind$loc.n.all)); names(N_loc)[i] <- species[i]
  
  load(paste(species[i],"/results/output.R",sep = ""))
  
  lims <- range(output_inds[,-c(1:2)],na.rm = T)
  cols <- heat.colors(11)[11:1]
  
  boxplot(output_inds[,c(4:ncol(output_inds),3)],col=c(rep("grey",ncol(output_inds[,-c(1:2)])-1),"red"),border=F,las=2,cex.axis=0.5,ylab="Offsets",main=species[i])
  names(output_inds) <- gsub("random_rep","rdm_rep",names(output_inds))
  names(output_inds) <- gsub("maxN","Nmax",names(output_inds))
  
  N <- names(output_inds)[grep("rep",names(output_inds))]  %>% strsplit(.,"_rep") %>% sapply(.,"[",1) %>% unique(.)
  nombres[[i]] <- N
  names(nombres)[i] <- paste("N_",length(N),sep = "")
  ranks <- vector("list",length(N))
  
  for( j in 1:length(N)){
    temp <- output_inds[grep(paste(N[j],"_",sep = ""),names(output_inds))]
    temp_rank <- NULL
    for(k in 1:ncol(temp)){
      corr <- cor.test(x=output_inds$all, y=temp[,k], method = 'spearman')
      temp_rank <- c(temp_rank,corr$estimate)
    }
    
    if(length(temp_rank)>1000){
      print(i)
      stop()
    }
    ranks[[j]]<- temp_rank
    names(ranks)[j]<- N[j]
    
    data <- rbind(data,
                  data.frame(sp=species[i],N_ind=N[j],rho=temp_rank))
    
  }
  
  #temp <- output_inds[grep(paste("N_random",sep = ""),names(output_inds))]
  #rands <- names(output_inds)[grep("random",names(output_inds))] %>% strsplit(.,"_rep") %>% sapply(.,"[",1) %>% unique(.)
  
  if(length(ranks$N_1)!=length(ranks$N_2)){
    Ntemp <- rep(NA,1000)
    Ntemp[1:length(ranks$N_1)] <- ranks$N_1
    ranks$N_1 <- Ntemp
  }
  boxplot(as.data.frame(ranks),las=2,cex.axis=0.5,ylab="Rank correlation",main=species[i])
  mat_list[[species[i]]] <- apply(as.data.frame(ranks),2,FUN = median,na.rm=T) 
  temp <- as.data.frame(ranks)[order(mat_list[[species[i]]])]
  
  boxplot(temp,las=2,cex.axis=0.5,ylab="Rank correlation",main=species[i])
  abline(v=grep("plus",names(temp)),col="grey",lty=2,lwd=2)
  
}
nombres <- names(nombres) %>% sub("N_","",.) %>% as.numeric() %>% max() %>% paste("N_",.,sep = "") %>% nombres[[.]] 

N <- length(nombres)
mat <- matrix(NA,nrow = N,ncol=length(species)) %>% data.frame(.)
rownames(mat) <- nombres 
names(mat) <- species 

for(i in 1:length(mat_list)){
  mat[names(mat_list[[i]]),names(mat_list)[i]] <- mat_list[[i]]
}


files_NA <- names(order_na[order(order_na)]) 
files_Nloc <- names(N_loc[order(N_loc,decreasing = T)]) 

library(plot.matrix)

colors = met.brewer(name="Hiroshige")
mat$temp <- NA; mat$temp <- rep(0:1,nrow(mat))[1:nrow(mat)]
plot(as.matrix(mat[,c(files_NA,"temp")]),ylab="",col=colors,las=2,xlab="",cex.axis=0.5,main="Ordered >NA")
plot(as.matrix(mat[,c(files_Nloc,"temp")]),ylab="",col=colors,las=2,xlab="",cex.axis=0.5,main="Ordered <Nloc")

dev.off()

data <- data[-1,]

data_N<- data_Nplus <- data

data_N <- data_N[-grep("N_plus",data_N$N_ind),]
data_Nplus <- data_Nplus[c(grep("N_10",data_Nplus$N_ind),grep("N_plus",data_Nplus$N_ind)),]


data_N$sp <- factor(data_N$sp,levels=files_Nloc)
data_Nplus$sp <- factor(data_Nplus$sp,levels=files_Nloc)

data_N$N_ind <- factor(data_N$N_ind,levels = unique(data_N$N_ind))
data_Nplus$N_ind <- factor(data_Nplus$N_ind,levels = unique(data_Nplus$N_ind))

colors = met.brewer(name="Hiroshige")

colors <- colors[c(10:6,5:1)]
ggplot(data_N, aes(x=sp, y=rho, fill=N_ind)) + 
  geom_boxplot(width = 1,outlier.size = 0.75)+
  ylab("Rank correlation") + 
  xlab("-NA->") +
  ggtitle("Individual samples")+
  theme_few() + scale_fill_manual(values=colors) -> p7

colors = met.brewer(name="Hiroshige")

colors <- colors[c(1,5:9)]

ggplot(data_Nplus, aes(x=sp, y=rho, fill=N_ind)) + 
  geom_boxplot(width = 1,outlier.size = 0.75)+
  ylab("Rank correlation") + 
  xlab("<-loc-") +
  ggtitle("Individual samples")+
  theme_few() + scale_fill_manual(values=colors) -> p8



ggsave("~/Desktop/p7.pdf" , plot = p7, 
       device = "pdf", width = 10, 
       height = 5, units= "in")
ggsave("~/Desktop/p8.pdf" , plot = p8, 
       device = "pdf", width = 10, 
       height = 5, units= "in")


individuals <- list(data_N=data_N,data_Nplus=data_Nplus)
save(individuals,file = "~/Desktop/proyecto_bosques_niebla/muestreos/manuscrito/correcciones/figures/lista_tabla_individuos.R")
