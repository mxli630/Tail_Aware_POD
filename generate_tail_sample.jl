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
#   - `mu` and `lambda` are updated.
#   - `z` is the value for where the tail of QoI starts. In our 
#     example, the tail is chosen as the top 1% of QoI.
#
# Outputs:
#   The `sampled_tails.jld2` file that contains the sampled 
#   parameters, snapshots, QoIs, and the importance sampling weight.
# ================================================================

include("supplement_functions.jl")


z = 0.0008

penalty_func_lambda(x, mu, lambda, z) = 0.05*(transpose(x)*x) - lambda*(f(x) - z) + (mu/2)*((f(x) - z)^2);
penalty_grad_lambda(x, mu, lambda, z) = 0.1.*x + (mu*(f(x)-z)-lambda).*grad_f(x);

# ALM method to update the constants mu and lambda
lambda_now = 100; mu_now = 1e4;
x = randn(parameter_len)
for i = 1:8
    res = SciPy.optimize.fmin_bfgs(penalty_func_lambda,x_now,fprime = penalty_grad_lambda, args = (mu_now,lambda_now,z), maxiter = 40, gtol = 1e-5, full_output = 1, retall = true);
    lambda_now = lambda_now - mu_now*(f(res[1])-z);
    mu_now = mu_now*2;
    x_now = res[1] + 0.5*randn(parameter_len);
end
lambda = lambda_now; mu = mu_now

# find the LDT optimizer
x = randn(parameter_len)
res = SciPy.optimize.fmin_bfgs(penalty_func_lambda,x,fprime = penalty_grad_lambda, args = (mu,lambda,z), maxiter = 40, gtol = 1e-5, full_output = 1, retall = true);
x_res = res[1];

# sample around the optimzer
x_tail, h_tail, f_tail, isweight_tail = sample_around_optimizer(x_res)

@save "sampled_tails.jld2" x_tail h_tail f_tail isweight_tail;
