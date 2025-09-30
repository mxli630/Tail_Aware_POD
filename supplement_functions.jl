# ================================================================
# supplemental_functions.jl
#
# Purpose:
#   Utility functions and wrappers for seeding RNG, data handling,
#   and calling the forward solver.
#
# Notes:
#   - Contains helper routines used across sampling and ROM scripts.
#   - QoI and its gradient are defined in `f(x)` and `grad_f(x)`.
# ================================================================
include(joinpath(@__DIR__, "external", "weeds_fracture.jl"))
using Statistics
using Plots
using SciPy
import KernelDensity
import Random
Random.seed!(1)
using FileIO, JLD2

isfreenode, nodei2freenodei, freenodei2nodei = DPFEHM.getfreenodes(length(Qs), dirichletnodes)
ns_reverse = [201,101]
nfreenode = sum(isfreenode)


f(x) = solveforh(x2logKs(x))[critical_point_node]*9.807*997*1e-6 #convert from head (meters) to pressure (MPa) (for water 25 C)
grad_f(x) = Zygote.gradient(f, x)[1]



function solveforhfree(logKs)
    Ks_neighbors = logKs2Ks_neighbors(logKs)
    isfreenode, nodei2freenodei, freenodei2nodei = DPFEHM.getfreenodes(length(Qs), dirichletnodes)
    args = (zeros(sum(isfreenode)), Ks_neighbors, neighbors, areasoverlengths, dirichletnodes, dirichleths, Qs, ones(length(Qs)), ones(length(Qs)))
    b = -DPFEHM.groundwater_residuals(args...)
    A = DPFEHM.groundwater_h(args...)
    return DPFEHM.amg_solver(A,b)
end

function full(h)
    h = map(i->isfreenode[i] ? h[nodei2freenodei[i]] : dirichleths[i], 1:length(Qs))
    return h
end

function sampler(nsamp,nfreenode = 19899)
    p = Array{Float64}(undef, nsamp, nfreenode)
    for i = 1:nsamp
        x = randn(3*num_eigenvectors+2);
        p[i,:] = transpose(solveforhfree(x2logKs(x)))
    end
    return p
end

function mode_gen(P)
    # P is a nsamp * nfreenode matrix 
    # each row of P is a solution
    nsamp = size(P,1);
    P = P./(sqrt(nsamp-1));
    F = svd(P);
    mode = F.Vt;
    lambda = (F.S).^2;
    return transpose(mode),lambda
end

function modegen_weight(P, weight)
    # P is a nsamp * nfreenode matrix 
    # each row of P is a solution
    # P contains unscaled samples
    nsamp = size(P,1)
    P = weight * P;
    P = P./(sqrt(nsamp-1))
    F = svd(P)
    Vtilde = transpose(F.Vt)
    eig = (F.S).^2
    return Vtilde, eig
end

function l2_norm(v)
    N = length(v);
    return norm(v,2)/sqrt(N);
end

function global_l2_error(h_red, h_truth)
    # assume that both h_red and h_truth are of length 19899
    # other entries are filled with 0 anyway
    err = h_truth - h_red;
    return l2_norm(err);
end

function sample_around_optimizer(x_res, nsamp = 1000)
    x_around = Array{Float64}(undef, 46, nsamp)
    h_around = Array{Float64}(undef, nfreenode, nsamp)
    isweight = Float64[]; all_f = Float64[];
    for i = 1:nsamp
        x = x_res + randn(46)
        h = solveforhfree(x2logKs(x))
        w_num = exp(-0.5*transpose(x)*x)
        w_denom = exp(-0.5*transpose(x-x_res)*(x-x_res))
        append!(isweight, w_num/w_denom)
        x_around[:,i] = x
        append!(all_f, full(h)[critical_point_node]*9.807*997*1e-6)
        h_around[:,i] = h
    end
    return x_around, h_around, all_f, isweight
end


function reduced_model_error(param,sample,nmode, mode)
    # param is a collection of nsamps of single parameters, 
    # should be (length of single parameter) x nsamps;
    # sample should be nfreenode x nsamps
    m = length(nmode);
    nsamp = size(param,2);
    gerror = Array{Float64}(undef,nsamp, m);
    perror = Array{Float64}(undef,nsamp, m);
    fp = Array{Float64}(undef, nsamp, m);
    for i = 1:nsamp
        truth = sample[:,i];
        truthl2 = l2_norm(truth);
        f_truth = f(param[:,i]);
        x_now = param[:,i];
        for j = 1:m
            hp = reduced_sln(x_now, mode[:,1:nmode[j]])
            fp[i,j] = full(hp)[critical_point_node]*9.807*997*1e-6;
            gerror[i,j] = global_l2_error(hp,truth)/truthl2;
            perror[i,j] = abs(fp[i,j] - f_truth)/f_truth;
        end
    end
    g_median = Float64[]; p_median = Float64[]
    for i = 1:m
        append!(g_median, median(gerror[:,i]))
        append!(p_median, median(perror[:,i]))
    end
    return gerror, perror, fp, g_median, p_median
end

