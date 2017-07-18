function [Difl,dDifl_dcl,dDifl_dT ]= diffusionCoef(cl,T,eff,mode)

Difl = exp(-5.6-124.3./(T-229-5e3*cl)-506*cl)*eff; % electrolyte diffusion [mm^2/s]
dDifl_dcl = [];
dDifl_dT = [];

if nargin == 4 && strcmp(mode,'partials')
    dDifl_dcl = Difl.*(5e3*124.3./(T-229-5e3*cl).^2-506);  
    dDifl_dT = Difl*124.3./(T-229-5e3*cl).^2;
end

end