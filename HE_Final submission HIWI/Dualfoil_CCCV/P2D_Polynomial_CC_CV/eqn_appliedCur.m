function [g, Jrow] = eqn_appliedCur(k,P,x)

% jn = x(P.idx_jn:P.nx:P.nx*P.nj, :);
% hx = P.L_neg/P.n_neg;
% area = 3*P.epss_neg/P.Rp_neg;

Jrow = zeros(1,P.nx*P.nj+2);
if isequal(P.mode, 'CC')
    g = 0;
    Jrow(P.nj*P.nx+1) = 1;
elseif isequal(P.mode, 'CV')
    curold = 0;
    vold = x(P.nj*P.nx+2,1);

    Rint = -(4.1-vold)/(x(P.nj*P.nx+1,k)-curold)+P.RG;
    OCP = x(P.nj*P.nx+2,k) + Rint*x(P.nj*P.nx+1,k);
    g = x(P.nj*P.nx+1,k) - (OCP-4.1)/Rint;

    Jrow(P.nj*P.nx+2) = 1/Rint;
    Jrow(P.nj*P.nx+1) = -1 + (Rint-P.RG)/Rint;
    
    
    
    
    % 2nd attempt
%     il = x(P.idx_il:P.nx:P.nx*P.nj, :); 
%     phis = x(P.idx_phis:P.nx:P.nx*P.nj, :);
%     hx = P.L_pos/P.n_pos;
%     g_pos = - P.L_pos*x(P.nj*P.nx+1,k)/P.sigma_pos + ...
%         hx*sum(il(P.bnd_pos_sep:P.nj-1,k)+il(P.bnd_pos_sep+1:P.nj,k))/2/P.sigma_pos;
%     
% %     Jrow([(P.nj-1)*P.nx+P.idx_il, (P.bnd_pos_sep-1)*P.nx+P.idx_il]) = -hx*P.sigma_pos/2;
% %     Jrow((P.bnd_pos_sep:P.nj-2)*P.nx+P.idx_il) = -hx*P.sigma_pos;
%     
%     hx = P.L_neg/P.n_neg;
%     g_neg = - P.L_neg*x(P.nj*P.nx+1,k)/P.sigma_neg + ...
%         hx*sum(il(1:P.bnd_sep_neg-1,k)+il(2:P.bnd_sep_neg,k))/2/P.sigma_neg;
%     
% %     Jrow([P.idx_il, (P.bnd_sep_neg-1)*P.nx+P.idx_il]) = -hx*P.sigma_neg/2;
% %     Jrow((1:P.bnd_sep_neg-2)*P.nx+P.idx_il) = -hx*P.sigma_neg;
% 
%     eta = overpot(x, k, P);
%     css = x(P.idx_css:P.nx:P.nx*P.nj, :);
%     jn = x(P.idx_jn:P.nx:P.nx*P.nj, :);
%     phil = x(P.idx_phil:P.nx:P.nx*P.nj, :);
%     [Eeq, dEeq_dcss] = Eeq_neg(css(P.bnd_sep_neg,k)/P.csmax_neg, P);
%     
% %     g = -0*phis(end,k)+0*phis(1,k)+ g_pos+g_neg + 0*phis(P.bnd_pos_sep,k)- 0.005;%-P.RG*x(P.nj*P.nx+1,k);
% %     Jrow(P.nj*P.nx+1) = P.L_pos/P.sigma_pos+P.L_neg/P.sigma_neg;%+P.RG;
% %     Jrow(P.nj*P.nx+1) = 1;
%     
%     g = -x(P.nj*P.nx+2,k) + g_pos+g_neg + phis(P.bnd_pos_sep,k)- ...
%         phil(P.bnd_sep_neg,k)-0.175193184028335-P.Rf_neg*P.F*jn(P.bnd_sep_neg,k)-eta(P.bnd_sep_neg)-P.RG*x(P.nj*P.nx+1,k);
%     Jrow(P.nj*P.nx+1) = -P.L_pos/P.sigma_pos-P.L_neg/P.sigma_neg-P.RG;
%     Jrow(P.nj*P.nx+2) = 1;
%     Jrow((P.bnd_pos_sep-1)*P.nx+P.idx_phis) = 1;
%     Jrow((P.bnd_sep_neg-1)*P.nx+P.idx_phil) = -1;
%     Jrow((P.bnd_sep_neg-1)*P.nx+P.idx_jn) = -P.Rf_neg*P.F;
%     Jrow((P.bnd_sep_neg-1)*P.nx+P.idx_css) = dEeq_dcss;
%     Jrow(P.idx_phis) = 1;
%     Jrow(P.idx_phis*P.nj) = -1;
    
    
end

