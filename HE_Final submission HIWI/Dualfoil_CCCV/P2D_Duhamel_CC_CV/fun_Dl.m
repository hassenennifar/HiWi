function [Dl, dDl] = fun_Dl(cl, P)
% cl in mol/m3

tdd = exp((P.EbarD)*(P.T-298.15)/(P.T*298.15));
Dl = 5.34e-10*exp(-0.65*cl./1000)*tdd;
dDl = -0.65*Dl./1000*tdd;