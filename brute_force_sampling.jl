include("supplement_functions.jl")
using FileIO, JLD2

numruns = 10000
paramlen = 3*num_eigenvectors+1
all_x = Array{Float64}(undef, paramlen, numruns)
all_h = Array{Float64}(undef, nfreenode, numruns)
all_f = Float64[]

for i = 1:numruns
  x = randn(paramlen)
  all_x[:,i] = x
  hfree = solveforhfree(x2logKs(x))
  all_h[:,i] = hfree
  push!(all_f, full(hfree)[critical_point_node]*9.807*997*1e-6)
end


FileIO.save("all_f.jld2","all_f",all_f)
FileIO.save("all_x.jld2","all_x",all_x)
FileIO.save("all_h.jld2","all_h",all_h)
