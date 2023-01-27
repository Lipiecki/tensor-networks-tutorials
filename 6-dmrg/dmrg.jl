using ITensors 
using KrylovKit: eigsolve # import only eigsolve function as other functions might conflict with ITensors

# An example of how to call dmrgHeisenberg() is provided in the run_dmrg.jl file

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

    sites = siteinds("S=1/2", N; conserve_qns=false)  
    
    os = OpSum()

    for j in 1:(N - 1)
        os += delta,"Sz",j,"Sz",j+1
        os += J/2,"S+",j,"S-",j+1
        os += J/2,"S-",j,"S+",j+1
    end

    H = MPO(os, sites)
    ψ = randomMPS(sites)
  
    sweeps = Sweeps(5)
    maxdim!(sweeps, 10, 20, 50, 100, 200)
    cutoff!(sweeps, 1E-10)

    PH = ProjMPO(H)
    for sw in 1:sweeps.nsweep
        for (pos, direction) in sweeporder(N)    
            
            orthogonalize!(ψ, pos)
            position!(PH, ψ, pos)

            ϕ = ψ[pos]*ψ[pos+1]

            vals, vecs = eigsolve(PH, ϕ, 1, ishermitian=true) # Krylov subspace method of finding the ground state

            #energy = vals[1] 
            ϕ::ITensor = vecs[1]

            U, S, V = svd(ϕ, inds(ψ[pos]), cutoff=sweeps.cutoff[sw], maxdim=sweeps.maxdim[sw])
            
            if direction == :right
                ψ[pos] = U
                ψ[pos+1] = S*V
            else
                ψ[pos] = U*S
                ψ[pos+1] = V
            end

            Epsi = apply(H, ψ)
            E = inner(ψ, Epsi)
            println("Sweep $sw, bond ($pos, $(pos+1)), direction $direction, energy = $E")
        end
    end
end