length(ARGS) == 10 || error("Format: julia infer_network.jl <ntaxa> <repid> <ils> <ngt> <L> <treechild-galled/nottreechild-galled/...> <density> <seqgen path> <iqtree path> <ASTRAL-IV path>")
ntaxa = parse(Int, ARGS[1])
repid = parse(Int, ARGS[2])
ils = ARGS[3]
ngt = parse(Int, ARGS[4])
L = parse(Int, ARGS[5])
folder = ARGS[6]
density = ARGS[7]
seqgen = abspath(ARGS[8])
iqtree = abspath(ARGS[9])
astral = abspath(ARGS[10])

isfile(seqgen) || error("seqgen path \"$(seqgen)\" does not exist or is not a file.")
isfile(iqtree) || error("iqtree path \"$(iqtree)\" does not exist or is not a file.")
isfile(astral) || error("astral path \"$(astral)\" does not exist or is not a file.")

ntaxa in [6, 10, 20, 30] || error("Invalid ntaxa.")
1 ≤ repid ≤ 25 || error("Invalid repid.")
ils in ["low", "high"] || error("Invalid ILS.")
ngt in [100, 1000, 10000] || error("Invalid ngt.")
L in [0, 500, 1000] || error("Invalid L.")
density in ["low", "high"] || error("Invalid density.")
folder in ["treechild-galled", "nottreechild-galled", "treechild-notgalled", "nottreechild-notgalled"] || error("Invalid folder.")
(L != 0 || ils == "low") || error("If ngt == 0, ils must be low.")
(L != 0 || density == "low") || error("If ngt == 0, density must be low.")

using Distributed

# Need to run these lines locally, but on Condor we call the script with
# `--project=.`, which is preferable on Condor. 
# @everywhere using Pkg
# @everywhere Pkg.activate(joinpath(@__DIR__, ".."))

@everywhere using Random, PhyloNetworks, SNaQ
using CSV, DataFrames, PhyloCoalSimulations, StatsBase


################################
# GENERATE INPUT DATA FOR SNAQ #
################################

q::Array{Float64} = Array{Float64}(undef, binomial(ntaxa, 4), 3)    # THIS IS THE `q` THAT WILL BE USED AS INPUT
global t0::HybridNetwork    # THIS IS THE STARTING TOPOLOGY THAT WILL BE USED
truenet = readnewick(readlines(joinpath(@__DIR__, "processed-n$(ntaxa)-$(density)-dens.netfile"))[repid])
nhyb = truenet.numhybrids

@info "Setting ILS."
global targetmean::Float64 = ils == "low" ? 2.0 : ils == "high" ? 0.5 : error("ILS must be low or high.")
global sumE::Float64 = 0.0
global nE::Int = 0
for E in truenet.edge
    global nE, sumE
    getchild(E).leaf && continue
    sumE += E.length == -1. ? 0.0 : E.length
    nE += E.length == -1. ? 0 : 1
end
global refactor::Float64 = targetmean / (sumE / nE)
for E in truenet.edge
    global refactor
    E.length = E.length == -1. ? -1. : E.length * refactor
end

@info "Generating true gene trees."
gts = simulatecoalescent(truenet, ngt, 1);
estgts = Array{HybridNetwork}(undef, ngt);

@info "Generating true gt --> sequences --> est gts"
tdir = abspath("./temp_dir-$(ntaxa)-$(repid)-$(ils)-$(ngt)-$(L)-$(density)-$(folder)/")
if isdir(tdir) rm(tdir, recursive=true) end
mkdir(tdir);

if L == 0
    # No simulated gene trees, so set `estgts` to true gts as a quick hack
    estgts .= gts
else
    # Here, we conduct the process in batches, with 1 gene tree going through the pipeline
    # at a time per processor allocated (say there are `np` processors allocated).
    estgts = @distributed (vcat) for j = 1:ngt
        # 1. write gt
        gtfile = joinpath(tdir, "truegts-$(j).tre")
        seqfile = joinpath(tdir, "seq-$(j).seq")
        writenewick(gts[j], gtfile)
        
        # 2. simulate sequence
        run(pipeline(`$(seqgen) -q -s0.036 -n1 -f0.3,0.2,0.2,0.3 -mHKY -l$(L) "$(gtfile)" -z$(j)`, seqfile))
        rm(gtfile)

        # 3. infer gt
        iqtreeprefix = joinpath(tdir, "iqtree-$(j)")
        run(`$(iqtree) -s "$(seqfile)" -pre "$(iqtreeprefix)" -seed $(j) -quiet`)
        rm(seqfile)

        # 4. read inferred gt
        estgt = readnewick("$(iqtreeprefix).treefile")
        for suffix in [".bionj", ".ckp.gz", ".log", ".model.gz", ".iqtree", ".mldist", ".treefile"]
            rm("$(iqtreeprefix)$(suffix)")
        end
        estgt
    end
end


@info "Computing GTEE."
gtee = mean([hardwiredclusterdistance(tnet, enet, false) / (2*tnet.numtaxa-6) for (tnet, enet) in zip(gts, estgts)])
@info "gtee = $(gtee) - $(gtee * (2*truenet.numtaxa-6))"


@info "Counting qCFs."
qstatic, _ = countquartetsintrees(estgts)


@info "Converting `qstatic` to a Float64 array."
for j = 1:length(qstatic)
    for k = 1:3
        q[j, k] = qstatic[j].data[k]
    end
end


@info "Inferring a starting tree topology with ASTRAL-IV."
estgtfile = joinpath(tdir, "estgts.tre")
writemultinewick(estgts, estgtfile)
t0file = joinpath(tdir, "t0.tre")
run(`$(astral) -v 0 -i "$(estgtfile)" -o $(t0file)`)
t0 = readnewick(t0file)
for node in t0.node
    if !node.leaf && !node.hybrid
        node.name = ""
    else
        try
            parse(Float64, node.name)
            node.name = ""
        catch e
        end
    end
end
for E in t0.edge if E.length == -1. E.length = 0.0 end end
rm(estgtfile)
@info "ASTRAL-IV starting topology: $(writenewick(t0, round=true))"


@info "Collecting garbage everywhere."
@everywhere GC.gc();


############
# RUN SNAQ #
############

@info "Running `multi_search`."
rt = @elapsed search_net, _ = multi_search(
    t0, q, nhyb; runs=20, restrictions=SNaQ.restriction_set(; require_strongly_tree_child=true),
    prehybprob=0.0, preopt=false, maxequivPLs=1500
)
println("\tTook $(round(rt / 60, digits=2)) minutes.")

@info "Computing PLs."
SNaQ.semidirect_network!(truenet)
truePL = SNaQ.compute_loss(truenet, q)
trueoptPL = optimize_bls!(truenet, q)
estPL = SNaQ.compute_loss(search_net, q)


@info "Saving data."
outfile = "n$(ntaxa)-r$(repid)-$(density)-dens.out"
outdir = "/mnt/ws/home/nkolbow/TreechildGalledSimulations/data/out/"
outdir = joinpath(outdir, "$(folder)")
if !isdir(outdir) mkdir(outdir) end
outdir = joinpath(outdir, "ils$(ils)")
if !isdir(outdir) mkdir(outdir) end
outdir = joinpath(outdir, "ngt$(ngt)")
if !isdir(outdir) mkdir(outdir) end
outdir = joinpath(outdir, "L$(L)")
if !isdir(outdir) mkdir(outdir) end
outfile = joinpath(outdir, outfile)
if !isfile(outfile) touch(outfile) end


CSV.write(
    outfile,
    DataFrame(
        truenet=writenewick(truenet),
        estnet=writenewick(search_net),
        truePL=truePL,
        trueoptPL=trueoptPL,
        estPL=estPL,
        gtee=gtee,
        rt=rt,
        ntaxa=ntaxa,
        repid=repid
    ),
    append=true
)

