using ITensors, DelimitedFiles

# Run this script in REPL and call imTimeEvo()
# or
# run run_time_evo.jl file from shell `julia run_time_evo.jl`

# For repeated calls avoid running `julia run_time_evo.jl` multiple times,
# because it creates a new instance of Julia and will compile the function
# each time. Either run the script in REPL and call the functions multiple times or
# modify the run_time_evo.jl file to perform multiple calls.

function imTimeEvo(; N::Integer = 20, delta::Real = 1.0, J::Real = 1.0, tau::Real = 0.1, ttotal::Real = 25.0, cutoff::Real = 1E-8)
  # Number of time steps
  steps = Int(ttotal/tau)
  # Array to store <Sz> at each time step
  Sz_n_t = zeros(N, steps)

  # Make an array of 'site' indices
  s = siteinds("S=1/2", N; conserve_qns=true)

  # Make gates (1,2),(2,3),(3,4),...
  gates = ITensor[]  
  os = OpSum()

  for j in 1:(N - 1)
    s1 = s[j]
    s2 = s[j + 1]

    hj =
      delta*op("Sz", s1) * op("Sz", s2) +
      J/2* op("S+", s1) * op("S-", s2) +
      J/2* op("S-", s1) * op("S+", s2)
    
    Gj = exp(-im * tau/2 * hj)
    push!(gates, Gj)

    # os is used for creating Hamiltonian
    os += delta,"Sz",j,"Sz",j+1
    os += J/2,"S+",j,"S-",j+1
    os += J/2,"S-",j,"S+",j+1
  end

  H = MPO(os, s)

  # Include gates in reverse order too
  # (N,N-1),(N-1,N-2),...
  append!(gates, reverse(gates))

  # Initialize psi to be a product state (alternating up and down)
  c = div(N, 2) # center site

  # domain-wall initial state
  psi = productMPS(s, n -> n < c ? "Up" : "Dn")

  # Compute and print <Sz> at each time step
  # then apply the gates to go to the next time step
  for step in 1:steps
    t = tau*step
    Epsi = apply(H, psi)
    E = inner(psi, Epsi)
    println("$t $E")
    for site in 1:N
      Sz_n_t[site, step] = expect(psi, "Sz"; sites = site)
    end

    t≈ttotal && break
    psi = apply(gates, psi; cutoff)
    normalize!(psi)
  end
  # save <Sz> at each time step
  writedlm("Sz_n_t.dat", Sz_n_t)
end