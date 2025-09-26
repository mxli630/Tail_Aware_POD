include("supplement_functions.jl")

penalty_func_lambda(x, mu, lambda, z) = 0.05*(transpose(x)*x) - lambda*(f(x) - z) + (mu/2)*((f(x) - z)^2);
penalty_grad_lambda(x, mu, lambda, z) = 0.1.*x + (mu*(f(x)-z)-lambda).*grad_f(x);

lambda = 1024; mu = 1e8; z = 0.0011
x = randn(3*num_eigenvectors+1)
res = SciPy.optimize.fmin_bfgs(penalty_func_lambda,x,fprime = penalty_grad_lambda, args = (mu,lambda,z), maxiter = 40, gtol = 1e-5, full_output = 1, retall = true);
x_res = res[1];

x_around, h_around, f_around, isweight = sample_around_optimizer(x_res)

RleftIS = Diagonal(isweight.*(f_around.^2))
V_aroundIS, eig_aroundIS = modegen_weight(transpose(h_around), RleftIS)

nmode = [10, 15, 20, 25, 30, 35, 40, 45, 50, 75, 100, 200];
all_x_test = FileIO.load("all_x_train.jld2","all_x");
all_h_test = FileIO.load("all_h_train.jld2","all_h");
all_f_test = FileIO.load("all_f.jld2","all_f");

sortind = sortperm(all_f_test, rev = true);
all_x_sorted = all_x_test[:,sortind];
all_h_sorted = all_h_test[:,sortind];

tail_x = all_x_sorted[:,1:100]
tail_h = all_h_sorted[:,1:100]

ge_aroundIS, pe_aroundIS, f_aroundIS, gm_aroundIS, pm_aroundIS = reduced_model_error(tail_x, tail_h, nmode, V_aroundIS)
