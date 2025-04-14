library(SiPhyNetwork)
library(tidyverse)

# CHOICES: (lambda, nu)
# n10, low density: (0.5, 0.03)
# n10, high density: (0.5, 0.08)
#
# n20, low density: (0.5, 0.02)
# n20, high density: (0.6, 0.05)
#
# n30, low density: (0.5, 0.0125)
# n30, high density: (0.6, 0.03)






###########################################
# THESE NETWORKS ARE TREECHILD AND GALLED #
###########################################

set.seed(1234)
params <- data.frame(
    ntaxa = rep(c(10, 20, 30), each=2),
    lambda = c(0.5, 0.5, 0.5, 0.6, 0.5, 0.6),
    nu = c(0.03, 0.08, 0.02, 0.05, 0.0125, 0.03),
    density = rep(c("low", "high"), 3)
) %>%
    mutate(filename = paste0("data/networks/treechild-galled/raw-n", ntaxa, "-", density, "-dens.netfile"))
nreps <- 25

for(i in 1:nrow(params)) {
    nets <- sim.bdh.taxa.ssa(   # PLACEHOLDER LIST
            params$ntaxa[i], 50, params$lambda[i], 0, params$nu[1],
            hybprops = c(1, 1, 1),
            hyb.inher.fxn = make.beta.draw(10, 10)
        )
    neti <- 1

    while(length(nets) < 10*nreps) { # 10*nreps in case some networks aren't galled
        net <- sim.bdh.taxa.ssa(
            params$ntaxa[i], 1, params$lambda[i], 0, params$nu[1],
            hybprops = c(1, 1, 1),
            hyb.inher.fxn = make.beta.draw(10, 10)
        )[[1]]
        if(nrow(net$reticulation) == 0 || !isTreeChild(net)) next;

        nets[[neti]] <- net
        neti <- neti + 1
    }
    
    write.net(nets, params$filename[i])
}
