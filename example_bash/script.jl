# include the files with functions to be used in the script
include("../6-dmrg/dmrg.jl")
@time begin  # block for timing the script
    dmrgHeisenberg(N = parse(Int, ARGS[1]), delta = parse(Float64, ARGS[2]), J = parse(Float64, ARGS[3]))
end