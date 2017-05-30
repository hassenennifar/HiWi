function [g, J] = approxBattery_EKF(k,P,x,cur,ts,sol)

J = zeros((P.nx-2)*P.nj); % Jacobian Matrix

g = zeros((P.nx-2)*P.nj,1);

j = 1; % node 1
[g(j:j*(P.nx-2)), J(j:j*(P.nx-2),j:(j+1)*(P.nx-2))] = jac_EKF(j,k,P,x(:,k-1:k),cur,ts,sol(:,k-1:k));

for j = 2:P.nj-1
    [g((j-1)*(P.nx-2)+1:j*(P.nx-2)), J((j-1)*(P.nx-2)+1:j*(P.nx-2),(j-2)*(P.nx-2)+1:(j+1)*(P.nx-2))] ...
        = jac_EKF(j,k,P,x(:,k-1:k),cur,ts,sol(:,k-1:k));
end

j = P.nj; % node P.nj
[g((j-1)*(P.nx-2)+1:j*(P.nx-2)), J((j-1)*(P.nx-2)+1:j*(P.nx-2), (j-2)*(P.nx-2)+1:j*(P.nx-2))] ...
    = jac_EKF(j,k,P,x(:,k-1:k),cur,ts,sol(:,k-1:k));