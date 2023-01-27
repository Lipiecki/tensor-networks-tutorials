using ITensors, LinearAlgebra

# Simple example of SVD and the impact of truncation size on the error 
let
    M = 
      [ 0.435839 0.223707 0.10;
        0.435839 0.223707 -0.10;
        0.223707 0.435839 0.10;
        0.223707 0.435839 -0.10 ]
    
    maxdim = minimum(size(M))
    U, d, V = svd(M)
    
    # Truncate to keep_vals singular values
    keep_vals = 2
    Dtrunc = zeros((maxdim, maxdim))
    
    for i in 1:keep_vals
        Dtrunc[i, i] = d[i]
    end

    Mtrunc = U*Dtrunc*(V')
    diff = norm((M - Mtrunc))
    println("|M - Mtrunc|^2 = ", diff^2)
    println("d = ", d)
    println("Kept ", keep_vals, " singular values")
end