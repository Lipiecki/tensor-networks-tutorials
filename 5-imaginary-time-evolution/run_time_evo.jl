include("time_evo.jl") #include the files with functions to be used in the script
@time begin  # block for timing the script
    imTimeEvo(N = 20, delta = 1.0, J = 1.0, tau = 0.1, ttotal = 25.0, cutoff = 1E-8)
end

# For calling the script and passing the arguments from the shell
#=
include("time_evo.jl")
@time begin  # block for timing the script
    imTimeEvo(N = parse(Int, ARGS[1]), delta = parse(Float64, ARGS[2]), J = parse(Float64, ARGS[3]), tau = parse(Float64, ARGS[4]), ttotal = parse(Float64, ARGS[5]), cutoff = parse(Float64, ARGS[6]))
end 
=#
