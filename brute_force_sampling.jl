include("supplement_functions.jl")
using FileIO, JLD2

numruns = 10000
paramlen = 3*num_eigenvectors+1
all_x = Array{Float64}(undef, paramlen, numruns)
all_f = Float64[]

for i = 1:numruns
  x = randn(paramlen)
  all_x[:,i] = x
  push!(all_f, f(x))
end

fsort = sort(all_f, rev = true);
sortind = sortperm(all_f, rev = true);
all_x_sorted = all_x[:,sortind]

FileIO.save("all_f.jld2","all_f",all_f)
FileIO.save("all_x.jld2","all_x",all_x)
