
global ntot::Int = 0
open("tcg.table", "w+") do f
    for ntaxa in [10, 20, 30]
        d = joinpath(@__DIR__, "..", "err", string(ntaxa))
        if !isdir(d) mkdir(d) end
        if !isdir(replace(d, "err" => "out")) mkdir(replace(d, "err" => "out")) end
        if !isdir(replace(d, "err" => "log")) mkdir(replace(d, "err" => "log")) end

        for repid = 1:25
            d = joinpath(@__DIR__, "..", "err", string(ntaxa), string(repid))
            if !isdir(d) mkdir(d) end
            if !isdir(replace(d, "err" => "out")) mkdir(replace(d, "err" => "out")) end
            if !isdir(replace(d, "err" => "log")) mkdir(replace(d, "err" => "log")) end

            for ils in ["low", "high"]
                d = joinpath(@__DIR__, "..", "err", string(ntaxa), string(repid), ils)
                if !isdir(d) mkdir(d) end
                if !isdir(replace(d, "err" => "out")) mkdir(replace(d, "err" => "out")) end
                if !isdir(replace(d, "err" => "log")) mkdir(replace(d, "err" => "log")) end

                for ngt in [100, 1_000, 10_000]
                    d = joinpath(@__DIR__, "..", "err", string(ntaxa), string(repid), ils, string(ngt))
                    if !isdir(d) mkdir(d) end
                    if !isdir(replace(d, "err" => "out")) mkdir(replace(d, "err" => "out")) end
                    if !isdir(replace(d, "err" => "log")) mkdir(replace(d, "err" => "log")) end

                    for L in [500, 1000]
                        d = joinpath(@__DIR__, "..", "err", string(ntaxa), string(repid), ils, string(ngt), string(L))
                        if !isdir(d) mkdir(d) end
                        if !isdir(replace(d, "err" => "out")) mkdir(replace(d, "err" => "out")) end
                        if !isdir(replace(d, "err" => "log")) mkdir(replace(d, "err" => "log")) end

                        for folder in ["treechild-galled"]
                            d = joinpath(@__DIR__, "..", "err", string(ntaxa), string(repid), ils, string(ngt), string(L), folder)
                        if !isdir(d) mkdir(d) end
                        if !isdir(replace(d, "err" => "out")) mkdir(replace(d, "err" => "out")) end
                        if !isdir(replace(d, "err" => "log")) mkdir(replace(d, "err" => "log")) end

                            for density in ["low", "high"]
                                d = joinpath(@__DIR__, "..", "err", string(ntaxa), string(repid), ils, string(ngt), string(L), folder, density)
                                if !isdir(d) mkdir(d) end
                                if !isdir(replace(d, "err" => "out")) mkdir(replace(d, "err" => "out")) end
                                if !isdir(replace(d, "err" => "log")) mkdir(replace(d, "err" => "log")) end

                                for _ = 1:10
                                    global ntot
                                    write(f, "$(ntaxa),$(repid),$(ils),$(ngt),$(L),treechild-galled,$(density)\n")
                                    ntot += 1
                                    ntot >= 50 && error("Quitting early after 50 queue.")
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end