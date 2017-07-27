function [Dl, dDl] = fun_Dl(cl, P)
% cl in mol/m3

% tdd = exp((P.EbarD)*(P.T-298.15)/(P.T*298.15));
% Dl = 5.34e-10*exp(-0.65*cl./1000)*tdd;
% dDl = -0.65*Dl./1000*tdd;

% cl = cl/1e3;
% Dl = 10.^(-4.43-54./(P.T-(229+5*cl))-0.22*cl)*1e-4;
% dDl = log(10)*Dl.*(-270./(P.T-5*cl-229).^2-0.22);
% dDl = dDl/1e3;

T = P.T;
eff = 1.5;
cl = cl*(1e3/1e9); % convert mol/m3 to mmol/mm3 
Dl = exp(-5.6-124.3./(T-229-5e3*cl)-506*cl)*eff; % electrolyte diffusion [mm^2/s]
Dl = Dl*1e-6; % convert from [mm2/s] to [m2/s]

dDl = 1e6*Dl.*(5e3*124.3./(T-229-5e3*cl).^2-506);
dDl =  dDl*(1e-6)/(1e-3/1e-9);% convert from [(mm2/s)/(mmol/mm3)]