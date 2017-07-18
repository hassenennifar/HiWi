function [gz, Jz] = jac_EKF(j,k,P,x,cur,ts,sol)

gz = zeros(P.nx-2,1);

if j == 1 || j == P.nj
    Jz = zeros(P.nx-2, (P.nx-2)*2);
else
    Jz = zeros(P.nx-2, (P.nx-2)*3);
end

%[gj(1), Jj(1,:)] = eqn_matBalLiquid(j,2,P,x);

[gz(1), Jz(1,:)] = eqn_ohmLawLiquid_EKF(j,2,P,x);

[gz(2), Jz(2,:)] = eqn_BVK_EKF(j,2,P,x,sol);

%[gz(3), Jz(3,:)] = eqn_matBalSolid_EKF(j,2,P,x,sol);

[gz(3), Jz(3,:)] = eqn_curDivLiquid_EKF(j,2,P,x,cur);

[gz(4), Jz(4,:)] = eqn_ohmLawSolid_EKF(j,2,P,x,cur);






