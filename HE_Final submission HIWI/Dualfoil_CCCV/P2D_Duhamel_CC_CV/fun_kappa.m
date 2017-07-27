function [kappa, dkappa] = fun_kappa(cl, P)

tdkap = exp((P.Ebarkap)*(P.T-298.15)/(P.T*298.15));

kappa = tdkap*(0.0911+1.9101*cl./1000 - ...
        1.052*((cl./1000).^2) + 0.1554*((cl./1000).^3));
    
dkappa = tdkap*(1.9101/1000 - 2*1.052*cl./1000./1000 ...
        + 0.1554*3*((cl./1000).^2)./1000);

    