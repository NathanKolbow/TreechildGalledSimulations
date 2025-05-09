# PURPOSE:
# 1. remove 2- and 3-cycles from all generated networks
# 2. remove galled networks
# 3. ensure that the final topologies are all fully identifiable

using PhyloNetworks, StatsBase


nreps = 25
ntaxas = [10, 20, 30]
net_dirs = [d for d in readdir(@__DIR__, join=true) if isdir(d)]

for net_dir in net_dirs
    length(readdir(net_dir)) > 0 || continue
    for density in ["low", "high"]
        net_files = [joinpath(net_dir, "raw-n$(j)-$(density)-dens.netfile") for j in ntaxas]

        for (ntaxa, file) in zip(ntaxas, net_files)
            if ntaxa == 10 && density == "high" continue end
            nets = readmultinewick(file)
            keepnets = Vector{HybridNetwork}()
            for (neti, net) in enumerate(nets)
                # 1. remove 2- and 3-cycles from the network
                while shrink2cycles!(net) || shrink3cycles!(net) continue end

                # 2. remove non-galled, non-TC networks
                if !isgalled(net) continue end
                if !istreechild(net)[1] continue end
                
                # 3. ensure that the final topologies are all fully identifiable
                net.numhybrids > 0 || continue
                push!(keepnets, net)
                length(keepnets) >= nreps && break
            end

            # 4. write out the result
            length(keepnets) == nreps || error("Only have $(length(keepnets)) / $(nreps) nets.")
            
            filename = replace(file, "raw-" => "processed-")
            @info "$(filename): $(mean(n.numhybrids for n in keepnets) / ntaxa)"
            writemultinewick(keepnets, filename)
        end
    end
end
