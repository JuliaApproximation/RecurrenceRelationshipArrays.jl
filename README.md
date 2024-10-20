# RecurrenceRelationshipArrays.jl
A Julia package for caching solutions to recurrence relationships


This package implements arrays associated with recurrence relationships
as implemented in [RecurrenceRelationships.jl](https://github.com/JuliaApproximation/RecurrenceRelationships.jl). 

# Matrix-valued Clenshaw

We can
construct multiplication matrices associated with orthogonal polynomials using `Clenshaw`. The follow represents multiplication by $U_0(x) + 2U_1(x) + 3U_2(x)$ acting on a Chebyshev-U expansion:
```julia
julia> using RecurrenceRelationshipArrays, FillArrays, InfiniteArrays

julia> X = SymTridiagonal(Fill(0.0,∞), Fill(1/2,∞)) # Jacobi matrix associated with Chebyshev U
ℵ₀×ℵ₀ SymTridiagonal{Float64, Fill{Float64, 1, Tuple{InfiniteArrays.OneToInf{Int64}}}} with indices OneToInf()×OneToInf():
 0.0  0.5   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   …  
 0.5  0.0  0.5   ⋅    ⋅    ⋅    ⋅    ⋅      
  ⋅   0.5  0.0  0.5   ⋅    ⋅    ⋅    ⋅      
  ⋅    ⋅   0.5  0.0  0.5   ⋅    ⋅    ⋅      
  ⋅    ⋅    ⋅   0.5  0.0  0.5   ⋅    ⋅      
  ⋅    ⋅    ⋅    ⋅   0.5  0.0  0.5   ⋅   …  
  ⋅    ⋅    ⋅    ⋅    ⋅   0.5  0.0  0.5     
 ⋮                        ⋮              ⋱  

julia> rec_U = Fill(2,∞), Zeros{Int}(∞), Ones{Int}(∞); # recurrence coefficients for Chebyshev U

julia> Clenshaw([1,2,3], rec_U..., X)
ℵ₀×ℵ₀ Clenshaw{Float64} with 3 degree polynomial:
 1.0  2.0  3.0   ⋅    ⋅    ⋅    ⋅    ⋅   …  
 2.0  4.0  2.0  3.0   ⋅    ⋅    ⋅    ⋅      
 3.0  2.0  4.0  2.0  3.0   ⋅    ⋅    ⋅      
  ⋅   3.0  2.0  4.0  2.0  3.0   ⋅    ⋅      
  ⋅    ⋅   3.0  2.0  4.0  2.0  3.0   ⋅      
  ⋅    ⋅    ⋅   3.0  2.0  4.0  2.0  3.0  …  
  ⋅    ⋅    ⋅    ⋅   3.0  2.0  4.0  2.0     
 ⋮                        ⋮              ⋱  
```


# Minimal solutions to 3-term recurrences

The type `RecurrenceArray` allows us to represent a minimal solution
of a recurrence relationship, i.e., the Stieltjes transform of the orthogonal polynomials. 
This automatically combines forward and backward recurrence, and so will also will work on the support of the measure. Here is a simple example:
```julia
julia> z = 1.0001; # point close to the support

julia> ξ = inv(z + sign(z)sqrt(z^2-1)); # exact formula for Stieltjes transform of sqrt(1-x^2). The next term is ξ^2.

julia> r = RecurrenceArray(z, rec_U, [ξ,ξ^2])
ℵ₀-element RecurrenceArray{Float64, 1, Float64, Fill{Int64, 1, Tuple{InfiniteArrays.OneToInf{Int64}}}, Zeros{Int64, 1, Tuple{InfiniteArrays.OneToInf{Int64}}}, Ones{Int64, 1, Tuple{InfiniteArrays.OneToInf{Int64}}}} with indices OneToInf():
 0.9859575108273005
 0.9721122131567664
 0.9584613379288637
 0.9450021549685467
 0.9317319724392233
 0.9186481363043878
 0.9057480297968131
 ⋮

julia> z = [0.1+0im, 1.0001, 10.0]; # can evaluate at multiple points

julia> ξ = @. inv(z + sign(z)sqrt(z^2-1));

julia> RecurrenceArray(z, rec_U, [ξ'; ξ'.^2])
ℵ₀×3 RecurrenceArray{ComplexF64, 2, Vector{ComplexF64}, Fill{Int64, 1, Tuple{InfiniteArrays.OneToInf{Int64}}}, Zeros{Int64, 1, Tuple{InfiniteArrays.OneToInf{Int64}}}, Ones{Int64, 1, Tuple{InfiniteArrays.OneToInf{Int64}}}} with indices OneToInf()×Base.OneTo(3):
       0.1+0.994987im  0.985958+0.0im    0.0501256+0.0im
     -0.98+0.198997im  0.972112+0.0im   0.00251258+0.0im
    -0.296-0.955188im  0.958461+0.0im  0.000125945-0.0im
    0.9208-0.390035im  0.945002+0.0im   6.31305e-6-0.0im
   0.48016+0.877181im  0.931732+0.0im   3.16446e-7+0.0im
 -0.824768+0.565471im  0.918648+0.0im    1.5862e-8-0.0im
 -0.645114-0.764087im  0.905748+0.0im  7.95095e-10-0.0im
          ⋮                            
```
