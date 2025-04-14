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
using CSV, DataFrames, PhyloCoalSimulations


################################
# GENERATE INPUT DATA FOR SNAQ #
################################

q::Array{Float64} = Array{Float64}(undef, binomial(ntaxa, 4), 3)    # THIS IS THE `q` THAT WILL BE USED AS INPUT
global t0::HybridNetwork    # THIS IS THE STARTING TOPOLOGY THAT WILL BE USED
truenet = readnewick(readlines(joinpath(@__DIR__, "processed-n$(ntaxa)-$(density)-dens.netfile"))[repid])
nhyb = truenet.numhybrids

@info "Generating true gene trees."
gts = simulatecoalescent(truenet, ngt, 1);
estgts = Array{HybridNetwork}(undef, ngt);

@info "Generating true gt --> sequences --> est gts"
tdir = joinpath("temp_dir-$(ntaxa)-$(repid)-$(ils)-$(ngt)-$(L)-$(density)-$(folder)")
if isdir(tdir) rm(tdir, recursive=true) end
mkdir(tdir);

if L == 0
    # No simulated gene trees, so set `estgts` to true gts as a quick hack
    estgts .= gts
else
    Threads.@threads for j = 1:length(gts)
        gt = gts[j]
        print("\r$(j) / $(ngt) inferred gene trees")

        # 1. write gt
        gtfile = joinpath(tdir, "truegts-$(j).tre")
        seqfile = joinpath(tdir, "seq-$(j).seq")
        writenewick(gts[j], gtfile)
        
        # 2. simulate sequence
        run(pipeline(`$(seqgen) -q -s0.036 -n1 -f0.3,0.2,0.2,0.3 -mHKY -l$(L) "$(gtfile)" -z$(j)`, seqfile))
        rm(gtfile)

        # 3. infer gt
        iqtreeprefix = joinpath(tdir, "iqtree-$(j)")
        run(`$(iqtree) -nt $(Threads.nthreads()) -s "$(seqfile)" -pre "$(iqtreeprefix)" -seed $(j) -quiet`)
        rm(seqfile)

        # 4. read inferred gt
        estgts[j] = readnewick("$(iqtreeprefix).treefile")
        for suffix in [".bionj", ".ckp.gz", ".log", ".model.gz", ".iqtree", ".mldist", ".treefile"]
            rm("$(iqtreeprefix)$(suffix)")
        end
    end
    println()
end

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
run(`$(astral) -v 0 -i "$(estgtfile)" -o $(t0file) -t $(Threads.nthreads())`)
t0 = readnewick(t0file)
for E in t0.edge if E.length == -1. E.length = 0.0 end end
rm(estgtfile)


@info "Collecting garbage everywhere."
@everywhere GC.gc();


############
# RUN SNAQ #
############

@info "Running `multi_search`."
rt = @elapsed search_net, _ = multi_search(
    t0, q, nhyb; runs=20, restrictions=SNaQ.restriction_set(; require_strongly_tree_child=true)
)
println("\tTook $(round(rt / 60, digits=2)) minutes.")

@info "Computing PLs."
SNaQ.semidirect_network!(truenet)
truePL = SNaQ.compute_loss(truenet, q)
trueoptPL = optimize_bls!(truenet, q)
estPL = SNaQ.compute_loss(search_net, q)


@info "Saving data."
outfile = "n$(ntaxa)-r$(repid)-$(density)-dens.out"
outdir = joinpath(@__DIR__, "../data/out/")
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
        rt=rt,
        ntaxa=ntaxa,
        repid=repid
    ),
    append=true
)

