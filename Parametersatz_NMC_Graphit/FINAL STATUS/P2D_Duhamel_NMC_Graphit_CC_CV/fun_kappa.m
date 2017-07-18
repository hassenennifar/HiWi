function [kappa, dkappa] = fun_kappa(cl, P)
%% Dualfoil5.f
% tdkap = exp((P.Ebarkap)*(P.T-298.15)/(P.T*298.15));
% 
% kappa = tdkap*(0.0911+1.9101*cl./1000 - ...
%         1.052*((cl./1000).^2) + 0.1554*((cl./1000).^3));
%     
% dkappa = tdkap*(1.9101/1000 - 2*1.052*cl./1000./1000 ...
%         + 0.1554*3*((cl./1000).^2)./1000);
    
%% Julius collocation method
% eff = 1.5;
% T = P.T;
% cl = cl*(1e3/1e9); % convert [mol/m3] to [mmol/mm3]
% factor = -10.5+6.68e2*cl+0.074*T-1.78e1*T*cl-6.96e-5*T.^2+2.8e-2*cl*T.^2+4.94e5*cl.^2-8.86e2*cl.^2*T;
% kappa =  0.1*cl.*factor.^2*eff; % electrolyte kappa [S/mm]
% kappa = kappa/1e-3; % adapt unit [S/mm] to [S/m]
% 
% dfactor_dcl = (6.68e2-1.78e1*T+2.8e-2*T.^2)+(4.94e5-8.86e2*T)*2*cl;
% dkappa = 0.1*eff*(factor.^2+2*cl.*factor.*dfactor_dcl);  
% dkappa = dkappa*(1/1e-3)/(1e-3/1e-9); % adapt unit from [(S/mm)/(mmol/mm3)] to [(S/m)/(mol/m3)]


%% disable
% kappa = ones(size(cl));
% dkappa = zeros(size(cl));

%% Daten von COMSOL Modell NMC
T = P.T;
cl = cl/1e3;

kappa = cl.*(-10.5+0.074*T-6.96e-5*T^2+0.668*cl-0.0178*cl*T ...
    +2.8e-5*cl*T^2+0.494*cl.^2-8.86e-4*cl.^2*T).^2*0.1;

dkappa = (-10.5+0.074*T-6.96e-5*T^2+0.668*cl-0.0178*cl*T ...
    +2.8e-5*cl*T^2+0.494*cl.^2-8.86e-4*cl.^2*T).^2*0.1 ...
    + 2*cl.*(-10.5+0.074*T-6.96e-5*T^2+0.668*cl-0.0178*cl*T ...
    + 2.8e-5*cl*T^2+0.494*cl.^2-8.86e-4*cl.^2*T)*0.1 ...
    .*(0.668-0.0178*T+ 2.8e-5*T^2+2*0.494*cl-2*8.86e-4*cl*T);
dkappa = dkappa/1e3;

