# ================================================================
# generate_tail_sample.jl
#
# Purpose:
#   Run LDT-guided sampling to generate "tail" parameter realizations
#   and evaluate the Quantity of Interest (QoI).
#
# Notes:
#   - The ALM function combining the rate function and the constraint $f(x)>=z$
#     is `penalty_func_lambda(x,mu,lambda)`
#   - `mu` and `lambda` are updated .
#   - `z` is the value for where the tail of QoI starts. In our 
#     example, the tail is chosen as the top 1% of QoI.
#   - The rate function depends on constants associated with the
#     log permeability field, QoI, definition of tail. If you
#     changed anything related, you must update the constants  
#      as well.
#
# Outputs:
#   The `sampled_tails.jld2` file that contains the sampled 
#   parameters, snapshots, QoIs, and the importance sampling weight.
# ================================================================

include("supplement_functions.jl")


z = 0.0011

penalty_func_lambda(x, mu, lambda, z) = 0.05*(transpose(x)*x) - lambda*(f(x) - z) + (mu/2)*((f(x) - z)^2);
penalty_grad_lambda(x, mu, lambda, z) = 0.1.*x + (mu*(f(x)-z)-lambda).*grad_f(x);

lambda = 1024; mu = 1e8;
x = randn(3*num_eigenvectors+1)
x = randn(3*num_eigenvectors+1)

# ALM method
res = SciPy.optimize.fmin_bfgs(penalty_func_lambda,x,fprime = penalty_grad_lambda, args = (mu,lambda,z), maxiter = 40, gtol = 1e-5, full_output = 1, retall = true);
x_res = res[1];

x_tail, h_tail, f_tail, isweight_tail = sample_around_optimizer(x_res)

@save "sampled_tails.jld2" x_tail h_tail f_tail isweight_tail;
