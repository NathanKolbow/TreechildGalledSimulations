using RCall, PhyloNetworks, StatsBase, DataFrames

const NREPS = 25
@rimport SiPhyNetwork
params = DataFrame(
    ntaxa = [10, 20, 20, 30, 30],
    lambda = [0.5, 0.5, 0.5, 0.5, 0.5],
    nu = [0.04, 0.08, 0.175, 0.05, 0.05],
    density = ["low", "low", "high", "low", "high"]
)
@rput params
R"set.seed(1234)"
R"library(SiPhyNetwork)"
R"library(ape)"


function validate_network(newick::String)::Union{HybridNetwork, Nothing}
    net = readnewick(newick)
    while shrink2cycles!(net) || shrink3cycles!(net) continue end
    if net.numhybrids == 0 return nothing end
    if !isgalled(net) return nothing end
    if !istreechild(net)[1] return nothing end
    return net
end



for (i, prow) in enumerate(eachrow(params))
    if i == 1 || i == 2 || i == 4 continue end # these are done and good!
    @rput i
    nets = []
    R"set.seed(1234)"

    # Need to do this so that the n10 networks match what is
    # currently in the file...
    R"""
        nets <- sim.bdh.taxa.ssa(
            params$ntaxa[i], 50, params$lambda[i], 0, params$nu[1],
            hybprops = c(1, 1, 1),
            hyb.inher.fxn = make.beta.draw(10, 10)
        )
        net <- 0
    """

    while length(nets) < NREPS
        print("\r$(prow["ntaxa"]) $(prow["density"]): $(length(nets))/$(NREPS)")
        R"""
            while(is.numeric(net) || nrow(net$reticulation) == 0 || !isTreeChild(net)) {
                net <- sim.bdh.taxa.ssa(
                    params$ntaxa[i], 1, params$lambda[i], 0, params$nu[i],
                    hybprops = c(1, 1, 1), hyb.inher.fxn = make.beta.draw(10, 10)
                )[[1]]
            }
            net_newick <- write.net(net)
            net <- 0
        """
        @rget net_newick

        if net_newick == "" continue end
        net = validate_network(net_newick)
        if net !== nothing
            push!(nets, net)
        end
        print("\r$(prow["ntaxa"]) $(prow["density"]): $(length(nets))/$(NREPS)")
    end
    println("\n\tdensity: $(mean(n.numhybrids for n in nets) / prow["ntaxa"])")

    writemultinewick(Vector{HybridNetwork}(nets), "treechild-galled/processed-n$(prow["ntaxa"])-$(prow["density"])-dens.netfile")
end

