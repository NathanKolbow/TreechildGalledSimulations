#!/bin/bash

# Need this for astral-IV
export PATH
. /etc/profile.d/modules.sh
module load gcc-13.2.0

# unzip the Julia tarball
echo "Unzipping tarball."
if [ ! -f "./julia-1.11.4/bin/julia" ]; then
    tar -xzf julia-1.11.4-linux-x86_64.tar.gz
fi
jl=./julia-1.11.4/bin/julia

# compile the Julia project from `Project.toml`
export JULIA_DEPOT_PATH=$PWD
echo "Compiling project."
JULIA_DEPOT_PATH=$PWD ${jl} --project=. ./compile_project.jl

# run the inference script
echo "Running inference script."
JULIA_DEPOT_PATH=$PWD ${jl} --project=. -p6 -t6 ./infer_TCG_network.jl "$@"
