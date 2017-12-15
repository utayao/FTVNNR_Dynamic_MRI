function xhat = Solve_TVLR(A, B, pars, F_gt)

pars.max_iter = 200;
pars.tol = 1e-4;
pars.verbose = 1;  % It prints out evaluation each step, may set to 0 for faster speed
pars.L =1; 
pars.t1 = 4;
pars.debug_output = 0;

xhat = TVLR_opt(A, B, pars, F_gt);