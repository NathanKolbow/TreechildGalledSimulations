# PURPOSE: generate the networks to be used in simulations
# OUTPUT: files "raw-n<ntaxa>.netfile" - newicks of generated networks before
#         various operations in "network_validator.jl" are applied to them, hence "raw"

library(SiPhyNetwork)
library(tidyverse)

lambdas <- c(0.5,0.6)
nus <- c(0.0125,0.02,0.03,0.05,0.07,0.08)
ntips<-c(20)#,20,30)

nreps<- 2000
set.seed(1234)
iter<-expand.grid(1:nreps,ntips,nus,lambdas)

empty_col<-rep(Inf,nrow(iter))
ret_data <- data.frame(ntaxa = empty_col,
                       nu    = empty_col,
                       lambda= empty_col,
                       rets  = empty_col,
                       level = empty_col,
                       blobs = empty_col,
                       TC    = empty_col)

i<-1
while(i <= nrow(iter)){

  nu <- iter[i,3]
  ntaxa <-iter[i,2]
  lambda <-iter[i,4]
  
  if((i %% 1000) == 0){
    print(paste(ntaxa,nu,lambda,i))
  }
  
  net<-sim.bdh.taxa.ssa(ntaxa,1,
                   lambda,0,nu,
                   hybprops = c(1,1,1),
                   hyb.inher.fxn = make.beta.draw(10,10))[[1]]
  nret<- nrow(net$reticulation)
  if(nret==0){
    next 
  }

  nrets <- nrow(net$reticulation)

  ##Compute the biconnected components for the level and number of blobs
  rt <- as.integer(length(net$tip.labels) + 1)
  edges <- rbind(net$edge, net$reticulation)
  hyb_nds <- net$reticulation[, 2]
  mode(edges) <- "integer"
  nNode <- length(net$tip.label) + net$Nnode
  blobs <- biconnectedComponents(edges, rt, nNode)
  net_level <- 1
  for (blob in blobs) {
    blob <- blob + 1
    blob_nds <- unique(as.vector(edges[blob, ]))
    blob_nds <- blob_nds[blob_nds %in% hyb_nds]
    net_level <- max(c(net_level, length(blob_nds)))
  }
  
  
  
 
  ret_data[i,]<- c(ntaxa,nu,lambda,nret,net_level,length(blobs),isTreeChild(net))
  i<-i+1
}



summed_dat <- ret_data %>% group_by(ntaxa,lambda,nu) %>%
  summarise(ave_ret_density = mean(rets)/mean(ntaxa),
            ave_ret = mean(rets),
            var_ret = var(rets)/n(),
            ave_level = mean(level),
            prop_level1 = sum(rets==1)/n(),
            prop_TC = sum(TC)/n())

view(summed_dat)


# CHOICES: (lambda, nu)
# n10, low density: (0.5, 0.03)
# n10, high density: (0.5, 0.08)
#
# n20, low density: (0.5, 0.02)
# n20, high density: (0.6, 0.05)
#
# n30, low density: (0.5, 0.0125)
# n30, high density: (0.6, 0.03)