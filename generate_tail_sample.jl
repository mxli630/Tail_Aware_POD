include("supplement_functions.jl")

penalty_func_lambda(x, mu, lambda, z) = 0.05*(transpose(x)*x) - lambda*(f(x) - z) + (mu/2)*((f(x) - z)^2);
penalty_grad_lambda(x, mu, lambda, z) = 0.1.*x + (mu*(f(x)-z)-lambda).*grad_f(x);

lambda = 1024; mu = 1e8; z = 0.0011
x = randn(3*num_eigenvectors+1)
x = randn(3*num_eigenvectors+1)
res = SciPy.optimize.fmin_bfgs(penalty_func_lambda,x,fprime = penalty_grad_lambda, args = (mu,lambda,z), maxiter = 40, gtol = 1e-5, full_output = 1, retall = true);
x_res = res[1];

x_tail, h_tail, f_tail, isweight_tail = sample_around_optimizer(x_res)

@save "sampled_tails.jld2" x_tail h_tail f_tail isweight_tail;
