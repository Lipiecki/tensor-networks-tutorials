using ITensors, LinearAlgebra

# Example of using SVD to compute von Neumann entropy of mixed state
let
    s1 = Index(2, "s1")
    s2 = Index(2, "s2")

    sing = ITensor(s1, s2) # singlet state
    prod = ITensor(s1, s2) # product state

    # set singlet state
    sing[s1(1), s2(2)] = 1.0/sqrt(2)
    sing[s1(2), s2(1)] = -1.0/sqrt(2)

    # set product state
    prod[s1(1), s2(2)] = 1.0

    psi = ITensor(s1, s2) # mixed state
    
    for mix in 0.0:0.1:1.0
        psi = (1-mix)*prod + mix*sing
        U, d, V = svd(psi, s1)
        if d[2, 2] ≈ 0.0
            entropy = -(d[1, 1]^2)*(log(ℯ, d[1, 1]^2)) # ℯ is typed as \euler
        else
            entropy = -((d[1, 1]^2)*(log(ℯ, d[1, 1]^2)) + (d[2, 2]^2)*(log(ℯ, d[2, 2]^2)))
        end
        println("mix = ", mix, ",  von Neumann entropy = ", entropy)
    end
end