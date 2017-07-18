function [sigma, dsigma] = fun_sigma(P)

sigma = [P.sigma_neg*ones(P.n_neg+1,1)*(1-P.epsl_neg)^1.5; zeros(P.n_sep-1,1); ...
    P.sigma_pos*ones(P.n_pos+1,1)*(1-P.epsl_pos)^1.5];
dsigma = zeros(P.nj,1);