# include the files with functions to be used in the script
include("dmrg.jl")
@time begin  # block for timing the script
    dmrgHeisenberg(N = 20, delta = 1.0, J = 1.0)
end

# For parsing the arguments from shell use:
#=
include("dmrg.jl")
@time begin  # block for timing the script
    dmrgHeisenberg(N = parse(Int, ARGS[1]), delta = parse(Float64, ARGS[2]), J = parse(Float64, ARGS[3]))
end
=#