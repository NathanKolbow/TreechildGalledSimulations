library(SiPhyNetwork)
library(tidyverse)

# CHOICES: (lambda, nu)
# n10, low density: (0.5, 0.03)
# n10, high density: NA / REMOVED
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
    ntaxa = c(10, 20, 20, 30, 30),
    lambda = c(0.5, 0.5, 0.5, 0.5, 0.5),
    nu = c(0.04, 0.05, 0.125, 0.05, 0.05),   # values further tweaked so that desire density is achieved AFTER the network_valiator.jl script
    nsim = c(50, 50, 250, 100, 400),
    density = c("low", "low", "high", "low", "high")
) %>%
    mutate(filename = paste0("treechild-galled/raw-n", ntaxa, "-", density, "-dens.netfile"))
nreps <- 25

#for(i in 1:nrow(params)) {
for(i in c(3)) {
    nets <- sim.bdh.taxa.ssa(   # PLACEHOLDER LIST
            params$ntaxa[i], 50, params$lambda[i], 0, params$nu[1],
            hybprops = c(1, 1, 1),
            hyb.inher.fxn = make.beta.draw(10, 10)
        )
    neti <- 1
    ntarget <- params$nsim[i]

    while(length(nets) < ntarget) { # 10*nreps in case some networks aren't galled
        net <- sim.bdh.taxa.ssa(
            params$ntaxa[i], 1, params$lambda[i], 0, params$nu[i],
            hybprops = c(1, 1, 1),
            hyb.inher.fxn = make.beta.draw(10, 10)
        )[[1]]
        if(nrow(net$reticulation) == 0 || !isTreeChild(net)) next;

        nets[[neti]] <- net
        neti <- neti + 1
        cat("\r", params$ntaxa[i], " ", params$density[i], ": ", neti-1, "/", ntarget)
    }
    cat("\n")
    
    write.net(nets, params$filename[i])
}
