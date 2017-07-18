function [D,dD_dcl,dD_dT]=Dl(cl,T,eff,mode)

D = exp(-5.6-124.3./(T-229-5e3*cl)-506*cl)*eff; % electrolyte diffusion [mm^2/s]
dD_dcl = [];
dD_dT = [];

if nargin == 4 && strcmp(mode,'partials')
    dD_dcl = D.*(5e3*124.3./(T-229-5e3*cl).^2-506);  
    dD_dT = D*124.3./(T-229-5e3*cl).^2;
end

end