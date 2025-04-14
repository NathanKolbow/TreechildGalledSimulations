# PURPOSE:
# 1. remove 2- and 3-cycles from all generated networks
# 2. remove galled networks
# 3. ensure that the final topologies are all fully identifiable

using PhyloNetworks


nreps = 25
ntaxas = [10, 20, 30]
net_dirs = [d for d in readdir(@__DIR__, join=true) if isdir(d)]

for net_dir in net_dirs
    length(readdir(net_dir)) > 0 || continue
    for density in ["low", "high"]
        net_files = [joinpath(net_dir, "raw-n$(j)-$(density)-dens.netfile") for j in ntaxas]

        for (ntaxa, file) in zip(ntaxas, net_files)
            nets = readmultinewick(file)
            keepnets = Vector{HybridNetwork}()
            for (neti, net) in enumerate(nets)
                # 1. remove 2- and 3-cycles from the network
                while shrink2cycles!(net) || shrink3cycles!(net) continue end

                # 2. remove galled networks
                if !isgalled(net) continue end
                
                # 3. ensure that the final topologies are all fully identifiable
                net.numhybrids > 0 || continue
                push!(keepnets, net)
                length(keepnets) >= nreps && break
            end

            # 4. write out the result
            length(keepnets) == nreps || error("Only have $(length(keepnets)) / $(nreps) nets.")
            writemultinewick(keepnets, replace(file, "raw-" => "processed-"))
        end
    end
end
