using ITensors

# This file contains two functions for calculating expected value of Sz operator at each site
# in the chain of spin 1/2 particles with transverse-field Ising Hamiltonian

# A. Run this script in REPL and call functions expectSzAtSites() and expectSz()
# or
# B. Call these functions through the run_magnetization.jl file

# A. will be faster for repeated calls, as the functions will be compiled only once

"""
expectSzAtSites(h, J, N)

    prints the expectation value of Sz operator at each site, 
    i.e. ⟨ψⱼ|Szⱼ|ψⱼ⟩ for ∈ [1, N], for a random state |ψ⟩,
    where:
    - h is the strength of the transverse magnetic field,
    - J is the strength of the nearest-neighbor interaction,
    - N is the number of sites in the chain

    > calling expectSzAtSites() uses default values of h = 0.0, J = 4.0, N = 50.
"""
function expectSzAtSites(; h::Real = 0.0, J::Real = 4.0, N::Integer = 50)

    sites = siteinds("S=1/2", N, conserve_qns=false)
    
    # in Julia instead of ampo = AutoMPO(sites) we use OpSum() and after adding the operators 
    # we construct the hamiltonian as H = MPO(ampo, sites)
    
    os = OpSum()
    for j in 1:N-1
        os += -J,"Sz",j,"Sz",j+1
    end
    for j in 1:N
        os += -h*2.0,"Sx",j
    end

    H = MPO(os, sites)          # Hamiltonian

    psi0 = randomMPS(sites)     # Initial Matrix Product State
    sweeps = Sweeps(5)          # number of sweeps to use in the dmrg procedure
    setmaxdim!(sweeps, 20)      # maximum bond dimension to use in the dmrg procedure
    setcutoff!(sweeps, 1E-10)   # cutoff for truncation in the dmrg procedure

    E, psi = dmrg(H, psi0, sweeps, outputlevel = 0) # dmrg procedure
    println("Ground state energy = $E")
    
    # Now we want to calculate the expectation value of Sz at each site
    for j in 1:N
        
        Szj = op(sites, "Sz", j)            # Sz operator acting on site j ≡ Szⱼ
        orthogonalize!(psi, j)              # move the center of orthogonality to site j
        bra_psi = dag(prime(psi[j]))        # bra of the MPS psi at site j, i.e. ⟨ψⱼ|
        Sz_expect = bra_psi*Szj*psi[j]      # ⟨Szⱼ⟩ = ⟨ψⱼ|Szⱼ|ψⱼ⟩
        
        println("Sz at site ", j, " = ", Sz_expect[1]) # [1] is needed because Sz_expect is a single-element Tensor
    end
end

"""
expectSz(h, J, N)

    prints the expectation value of Sz operator on all sites, divided by the number of sites,
    i.e. ⟨ψ|Sz|ψ⟩/N, for a random state |ψ⟩,
    where:
    - h is the strength of the transverse magnetic field,
    - J is the strength of the nearest-neighbor interaction,
    - N is the number of sites in the chain

    > calling expectSz(h) uses default values of h = 0.0, J = 4.0, N = 50.
"""
function expectSz(; h::Real = 0.0, J::Real = 4.0, N::Integer = 50)
    
    sites = siteinds("S=1/2", N, conserve_qns=false)

    os = OpSum()
    for j in 1:N-1
        os += -J,"Sz",j,"Sz",j+1
    end
    for j in 1:N
        os += -h*2.0,"Sx",j
    end

    H = MPO(os, sites)          # Hamiltonian

    psi0 = randomMPS(sites)     # Initial Matrix Product State
    sweeps = Sweeps(5)          # number of sweeps to use in the dmrg procedure
    setmaxdim!(sweeps, 20)      # maximum bond dimension to use in the dmrg procedure
    setcutoff!(sweeps, 1E-10)   # cutoff for truncation in the dmrg procedure

    E, psi = dmrg(H, psi0, sweeps, outputlevel = 0) # dmrg procedure
    println("Ground state energy = $E")
    
    os = OpSum()
    
    for j in 1:N
        os += "Sz", j
    end
    
    Sz = MPO(os, sites)
    Szpsi = apply(Sz, psi)          # Apply the MPO Sz to the MPS psi, i.e. Sz|ψ⟩
    Sz_expect = inner(psi, Szpsi)   # Apply the inner product to the MPS |ψ> and the MPS Sz|ψ⟩, i.e. ⟨ψ|Sz|ψ⟩
    println("<Sz>/N = ", Sz_expect/N)

end
