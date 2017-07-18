function [k,dk_dcl,dk_dT]=kappa(cl,T,eff,mode)

factor = -10.5+6.68e2*cl+0.074*T-1.78e1*T*cl-6.96e-5*T.^2+2.8e-2*cl*T.^2+4.94e5*cl.^2-8.86e2*cl.^2*T;
k =  0.1*cl.*factor.^2*eff; % electrolyte kappa [S/mm]
%COMSOL: 0.1*cl*(-10.5+6.68e-1*cl+0.074*T-1.78e-2*T*cl-6.96e-5*T^2+2.8e-5*cl*T^2+4.94e-1*cl^2-8.86e-4*cl^2*T)^2
dk_dcl = [];
dk_dT = [];

if nargin == 4 && strcmp(mode,'partials')
    dfactor_dcl = (6.68e2-1.78e1*T+2.8e-2*T.^2)+(4.94e5-8.86e2*T)*2*cl;
    dfactor_dT = (0.074-1.78e1*cl)-(6.96e-5+2.8e-2*cl-8.86e2*cl.^2)*2*T;
    dk_dcl = 0.1*eff*(factor.^2+2*cl.*factor.*dfactor_dcl);  
    dk_dT = 0.1*eff*cl.*factor.*dfactor_dT; 
end

end