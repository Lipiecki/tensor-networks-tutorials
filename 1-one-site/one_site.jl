using ITensors, LinearAlgebra

let
    # Deifne our Index
    s = Index(2, "s")

    # Define Operator Sx
    Sx = ITensor(s, prime(s))
    Sx[s[1], prime(s)[2]] = 0.5
    Sx[s[2], prime(s)[1]] = 0.5

    # Define single-site wavefunction
    psi = ITensor(s)
    psi[s[1]] = 1.0
    psi[s[2]] = 1.0

    # Normalize psi
    normalize!(psi) # or psi = psi/norm(psi)
    println("<ψ|ψ>: $(dot(psi, psi))")

    phi = Sx*psi
    normalize!(phi)
    println("<psi|phi> = $(dot(psi, phi))")
end