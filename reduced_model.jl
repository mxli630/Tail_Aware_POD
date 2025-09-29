include("supplement_functions.jl")

@load "sampled_tails.jld2" x_tail h_tail f_tail isweight_tail;

RleftIS = Diagonal(isweight_tail.*(f_tail.^2))
V_tailIS, eig_tailIS = modegen_weight(transpose(h_tail), RleftIS)

nmode = [10, 15, 20, 25, 30, 35, 40, 45, 50, 75, 100, 200];
all_x_test = FileIO.load("all_x_train.jld2","all_x");
all_h_test = FileIO.load("all_h_train.jld2","all_h");
all_f_test = FileIO.load("all_f.jld2","all_f");

sortind = sortperm(all_f_test, rev = true);
all_x_sorted = all_x_test[:,sortind];
all_h_sorted = all_h_test[:,sortind];

tail_x = all_x_sorted[:,1:100]
tail_h = all_h_sorted[:,1:100]

ge_aroundIS, pe_aroundIS, f_aroundIS, gm_aroundIS, pm_aroundIS = reduced_model_error(tail_x, tail_h, nmode, V_tailIS)
