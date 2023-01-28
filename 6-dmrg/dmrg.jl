using ITensors
using KrylovKit: eigsolve # import only eigsolve function as other functions might conflict with ITensors

# A. Run this script in REPL and call dmrgHeisenberg()
# or
# B. Run run_dmrg.jl file from shell (julia run_dmrg.jl)

# For repeated calls avoid running (julia run_dmrg.jl) multiple times,
# because it creates a new instance of Julia and will compile the function
# each time. Either run the script in REPL and call the functions multiple times or
# modify the run_dmrg.jl file to perform multiple calls.

"""
sweeporder(n) returns a list of tuples that contains 
the site indices and the direction of the sweep (:left or :right) for a
DMRG sweep of a chain of length n starting from left to right, 
then from right to left.
"""
function sweeporder(n)
    result = []
    for i in 1:n-1
        push!(result, (i, :right))
    end
    for i in n-1:-1:1
        push!(result, (i, :left))
    end
    return result
end

"""
dmrgHeisenberg(; N = 20, delta = 1.0, J = 1.0)

Run DMRG on the Heisenberg model with N sites
"""
function dmrgHeisenberg(; N::Integer = 20, delta::Real = 1.0, J::Real = 1.0)
    println("Running DMRG on the Heisenberg model with N = $N, delta = $delta, J = $J")
    sites = siteinds("S=1/2", N; conserve_qns=false)  
    
    os = OpSum()

    for j in 1:(N - 1)
        os += delta,"Sz",j,"Sz",j+1
        os += J/2,"S+",j,"S-",j+1
        os += J/2,"S-",j,"S+",j+1
    end

    H = MPO(os, sites)
    psi = randomMPS(sites)
  
    sweeps = Sweeps(5)
    maxdim!(sweeps, 10, 20, 50, 100, 200)
    cutoff!(sweeps, 1E-10)

    PH = ProjMPO(H)
    for sw in 1:sweeps.nsweep
        for (pos, direction) in sweeporder(N)    
            
            orthogonalize!(psi, pos)
            position!(PH, psi, pos)

            phi = psi[pos]*psi[pos+1]

            vals, vecs = eigsolve(PH, phi, 1, :LM, ishermitian=true) # Krylov subspace method of finding the ground state

            #energy = vals[1] 
            phi::ITensor = vecs[1]

            U, S, V = svd(phi, inds(psi[pos]), cutoff=sweeps.cutoff[sw], maxdim=sweeps.maxdim[sw])
            
            if direction == :right
                psi[pos] = U
                psi[pos+1] = S*V
            else
                psi[pos] = U*S
                psi[pos+1] = V
            end

            Epsi = apply(H, psi)
            E = inner(psi, Epsi)
            println("Sweep $sw, bond ($pos, $(pos+1)), direction $direction, energy = $E")
        end
    end
end