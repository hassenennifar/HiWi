function [g, J] = approxBattery(k,P,x,cur,ts)

J = zeros(P.nx*P.nj+2); % Jacobian Matrix

g = zeros(P.nx*P.nj+2,1);

j = 1; % node 1
[g(j:j*P.nx), J(j:j*P.nx,j:(j+1)*P.nx)] = jac(j,k,P,x,cur,ts);
% [sigma, ~] = fun_sigma(P);
% if j < P.bnd_sep_neg % negative electrode
%     hx = P.L_neg/P.n_neg;
%     J(j*P.nx, P.nx*P.nj+1) = hx/sigma(j);
% end

for j = 2:P.nj-1
    [g((j-1)*P.nx+1:j*P.nx), J((j-1)*P.nx+1:j*P.nx,(j-2)*P.nx+1:(j+1)*P.nx)] ...
        = jac(j,k,P,x,cur,ts);
%     if j <= P.bnd_pos_sep && j >= P.bnd_sep_neg % separator 
%         J(j*P.nx-1, P.nx*P.nj+1) = -1;    
%     end
%     
%     [sigma, ~] = fun_sigma(P);
%     if j < P.bnd_sep_neg % negative electrode
%         hx = P.L_neg/P.n_neg;
%         J(j*P.nx, P.nx*P.nj+1) = hx/sigma(j);
%     elseif j == P.bnd_sep_neg % set constant current at separtor il=cur
%         hx = P.L_neg/P.n_neg;
%         J(j*P.nx, P.nx*P.nj+1) = hx/sigma(j);
%     elseif j < P.nj && j >= P.bnd_pos_sep % positive electrode
%         hx = P.L_pos/P.n_pos;
%         J(j*P.nx, P.nx*P.nj+1) = hx/sigma(j);
%     end
end

j = P.nj; % node P.nj
[g((j-1)*P.nx+1:j*P.nx), J((j-1)*P.nx+1:j*P.nx, (j-2)*P.nx+1:j*P.nx)] ...
    = jac(j,k,P,x,cur,ts);
% hx = P.L_pos/P.n_pos;
% J(j*P.nx, P.nx*P.nj+1) = hx/sigma(j);




[g(P.nj*P.nx+1), J(P.nj*P.nx+1,:)] = eqn_appliedCur(k,P,x);
[g(P.nj*P.nx+2), J(P.nj*P.nx+2,:)] = eqn_outVoltage(k,P,x);