#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8
#SBATCH --time=0:10:00

source /usr/local/sbin/modules.sh
module load julia
julia script.jl 20 1.0 1.0
