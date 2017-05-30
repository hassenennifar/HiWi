function [g, J] = approxBattery(k,P,x,cur,ts,sol)

J = zeros(P.nx*P.nj); % Jacobian Matrix

g = zeros(P.nx*P.nj,1);

j = 1; % node 1
[g(j:j*P.nx), J(j:j*P.nx,j:(j+1)*P.nx)] = jac(j,k,P,x(:,k-1:k),cur,ts(k-1:k),sol(:,k-1:k));

for j = 2:P.nj-1
    [g((j-1)*P.nx+1:j*P.nx), J((j-1)*P.nx+1:j*P.nx,(j-2)*P.nx+1:(j+1)*P.nx)] ...
        = jac(j,k,P,x(:,k-1:k),cur,ts(k-1:k),sol(:,k-1:k));
end

j = P.nj; % node P.nj
[g((j-1)*P.nx+1:j*P.nx), J((j-1)*P.nx+1:j*P.nx, (j-2)*P.nx+1:j*P.nx)] ...
    = jac(j,k,P,x(:,k-1:k),cur,ts(k-1:k),sol(:,k-1:k));