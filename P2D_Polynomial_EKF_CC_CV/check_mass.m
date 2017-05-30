% check mass conservation in solid
    c =  sol(P.idx_c:P.ns:P.ns*P.nj, k);
    css =  x(P.idx_css:P.nx:P.nx*P.nj, k);
    nLi_neg = P.epss_neg*sum(c(1:P.bnd_sep_neg-1)+c(2:P.bnd_sep_neg))/2*P.L_neg/P.n_neg;
    nLi_pos = P.epss_pos*sum(c(P.bnd_pos_sep:P.nj-1)+c(P.bnd_pos_sep+1:P.nj))/2*P.L_pos/P.n_pos;
    
    nLi(k-1) = nLi_neg+nLi_pos;
    
    q =  sol(P.idx_q:P.ns:P.ns*P.nj, k);
    qLi_neg = P.epss_neg*sum(q(1:P.bnd_sep_neg-1)+q(2:P.bnd_sep_neg))/2*P.L_neg/P.n_neg;
    qLi_pos = P.epss_pos*sum(q(P.bnd_pos_sep:P.nj-1)+q(P.bnd_pos_sep+1:P.nj))/2*P.L_pos/P.n_pos;
    qLi(k-1) = qLi_neg+qLi_pos;
        
    c =  sol_est(P.idx_c:P.ns:P.ns*P.nj, k);
    css =  x_est(P.idx_css:P.nx:P.nx*P.nj, k);
    nLi_neg = P.epss_neg*sum(c(1:P.bnd_sep_neg-1)+c(2:P.bnd_sep_neg))/2*P.L_neg/P.n_neg;
    nLi_pos = P.epss_pos*sum(c(P.bnd_pos_sep:P.nj-1)+c(P.bnd_pos_sep+1:P.nj))/2*P.L_pos/P.n_pos;
    
    nLi_est(k-1) = nLi_neg+nLi_pos;
    
    q =  sol_est(P.idx_q:P.ns:P.ns*P.nj, k);
    qLi_neg = P.epss_neg*sum(q(1:P.bnd_sep_neg-1)+q(2:P.bnd_sep_neg))/2*P.L_neg/P.n_neg;
    qLi_pos = P.epss_pos*sum(q(P.bnd_pos_sep:P.nj-1)+q(P.bnd_pos_sep+1:P.nj))/2*P.L_pos/P.n_pos;
    qLi_est(k-1) = qLi_neg+qLi_pos;
    
    %disp(nLi)
    
%     css =  x1(P.idx_css:P.nx:P.nx*P.nj, k);
%     nLi_neg = P.epss_neg*sum(css(1:P.bnd_sep_neg-1)+css(2:P.bnd_sep_neg))/2*P.L_neg/P.n_neg;
%     nLi_pos = P.epss_pos*sum(css(P.bnd_pos_sep:P.nj-1)+css(P.bnd_pos_sep+1:P.nj))/2*P.L_pos/P.n_pos;
%     
%     nLi1 = nLi_neg+nLi_pos;
    %disp(nLi1/nLi_est)    
    
    cl =  x_est(P.idx_cl:P.nx:P.nx*P.nj, k);
    nLi_neg = P.epsl_neg*sum(cl(1:P.bnd_sep_neg-1)+cl(2:P.bnd_sep_neg))/2*P.L_neg/P.n_neg;
    nLi_pos = P.epsl_pos*sum(cl(P.bnd_pos_sep:P.nj-1)+cl(P.bnd_pos_sep+1:P.nj))/2*P.L_pos/P.n_pos;
    nLi_sep = P.epsl_sep*sum(cl(P.bnd_sep_neg:P.bnd_pos_sep-1)+cl(P.bnd_sep_neg+1:P.bnd_pos_sep))/2*P.L_sep/P.n_sep;
    
    nl_Li_est(k-1) = nLi_neg+nLi_sep+nLi_pos;
    %disp(nl_Li)
    
    cl =  x(P.idx_cl:P.nx:P.nx*P.nj, k);
    nLi_neg = P.epsl_neg*sum(cl(1:P.bnd_sep_neg-1)+cl(2:P.bnd_sep_neg))/2*P.L_neg/P.n_neg;
    nLi_pos = P.epsl_pos*sum(cl(P.bnd_pos_sep:P.nj-1)+cl(P.bnd_pos_sep+1:P.nj))/2*P.L_pos/P.n_pos;
    nLi_sep = P.epsl_sep*sum(cl(P.bnd_sep_neg:P.bnd_pos_sep-1)+cl(P.bnd_sep_neg+1:P.bnd_pos_sep))/2*P.L_sep/P.n_sep;
    nl_Li(k-1) = nLi_neg+nLi_sep+nLi_pos;
    %disp(nl_Li)