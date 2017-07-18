function [Dl, dDl] = fun_Dl(cl, P)
% cl in mol/m3

% tdd = exp((P.EbarD)*(P.T-298.15)/(P.T*298.15));
% Dl = 5.34e-10*exp(-0.65*cl./1000)*tdd;
% dDl = -0.65*Dl./1000*tdd;

cl = cl/1e3;
Dl = 10.^(-4.43-54./(P.T-(229+5*cl))-0.22*cl)*1e-4;
dDl = log(10)*Dl.*(-270./(P.T-5*cl-229).^2-0.22);
dDl = dDl/1e3;