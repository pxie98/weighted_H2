function [fhist, x, f, alg, func, iter_hist] = dfo1_4test(max_iter, func, x_0, clist)
  if (~exist('dfoipv', 'dir'))
    warning('dfoipv seems unavailable on your computer. Forget about testing it.');
    x = x_0;
    f = fun(x_0);
  else
    options = struct("max_iter", max_iter, "maxfev", 100*length(x_0), "init_delta", 1, "tol_delta", 1e-6, "tol_f", 1e-6, "tol_norm_g", 1e-6, "sample_gen", "auto", "verbosity", 0);
    alg = "DFO";
    cd 'dfoipv'
    % [res, ~, fhist] = bb_optimize(func, x_0, alg, [1/3, 1/3, 1/3], options);
    [res, ~, fhist, iter_hist] = bb_optimize(func, x_0, alg, clist, options);
    cd '..'
    x = res.x;
    f = func(x);
  end
end
