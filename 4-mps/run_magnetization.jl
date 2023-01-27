include("magnetization.jl") #include the files with functions to be used in the script
@time begin  # block for timing the script
    expectSzAtSites(h = 0.0, J = 4.0, N = 20)
end

# For parsing the arguments from shell use:
#=
include("magnetization.jl")
@time begin  # block for timing the script
    expectSzAtSites(h = parse(Float64, ARGS[1]), J = parse(Float64, ARGS[2]), N = parse(Int, ARGS[3]))    
end 
=#