using ITensors, LinearAlgebra

let
    # Initial product state
    s1 = Index(2, "s1")
    s2 = Index(2, "s2")
    psi = ITensor(s1, s2)
    psi[s1(1), s2(2)] = 1.0

    # S⁺ on first site
    Sp1 = ITensor(s1, prime(s1))
    Sp1[s1(1), prime(s1)(2)] = 1.0
    Sp1[s1(2), prime(s1)(1)] = 0.0 

    # S⁻ on first site
    Sm1 = ITensor(s1, prime(s1))
    Sm1[s1(1), prime(s1)(2)] = 0.0
    Sm1[s1(2), prime(s1)(1)] = 1.0

    # S⁺ on second site
    Sp2 = ITensor(s2, prime(s2))
    Sp2[s2(1), prime(s2)(2)] = 1.0
    Sp2[s2(2), prime(s2)(1)] = 0.0

    # S⁻ on second site
    Sm2 = ITensor(s2, prime(s2))
    Sm2[s2(1), prime(s2)(2)] = 0.0
    Sm2[s2(2), prime(s2)(1)] = 1.0

    # Sz on first site
    Sz1 = ITensor(s1, prime(s1))
    Sz1[s1(1), prime(s1)(1)] = 0.5
    Sz1[s1(2), prime(s1)(2)] = -0.5

    # Sz on second site
    Sz2 = ITensor(s2, prime(s2))
    Sz2[s2(1), prime(s2)(1)] = 0.5
    Sz2[s2(2), prime(s2)(2)] = -0.5

    # Two-site Heisenberg Hamiltonian
    H = Sz1*Sz2 + 0.5*Sp1*Sm2 + 0.5*Sm1*Sp2


    β = 100.0
    expH = exp(-β*H, ishermitian=true)

    psibeta = expH*psi
    psibeta = noprime(psibeta)
    normalize!(psibeta)
    
    println("psibeta = ", psibeta)
    println("En = ", (dag(prime(psibeta))*H*psibeta)[1]) # [1] is used to extract the value from 0-dimensional ITensor

    # SVD
    U, S, V = svd(psibeta, s1, maxdim=1)
    newpsi = U*S*V
    overlap = dot(psibeta, newpsi)
    println("overlap for maxdim = 1: ", overlap)

    U, S, V = svd(psibeta, s1, maxdim=2)
    newpsi = U*S*V
    overlap = dot(psibeta, newpsi)
    println("overlap for maxdim = 2: ", overlap)
end