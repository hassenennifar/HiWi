function [gj, Jj] = jac(j,k,P,x,cur,ts,sol)

gj = zeros(P.nx,1);

if j == 1 || j == P.nj
    Jj = zeros(P.nx, P.nx*2);
else
    Jj = zeros(P.nx, P.nx*3);
end

[gj(1), Jj(1,:)] = eqn_matBalLiquid(j,2,P,x);

[gj(2), Jj(2,:)] = eqn_ohmLawLiquid(j,2,P,x);

[gj(3), Jj(3,:)] = eqn_BVK_Eigen(j,2,P,x,sol,ts);

[gj(4), Jj(4,:)] = eqn_matBalSolid_Eigen(j,2,P,x,sol,ts);

[gj(5), Jj(5,:)] = eqn_curDivLiquid(j,2,P,x,cur);

[gj(6), Jj(6,:)] = eqn_ohmLawSolid(j,2,P,x,cur);






