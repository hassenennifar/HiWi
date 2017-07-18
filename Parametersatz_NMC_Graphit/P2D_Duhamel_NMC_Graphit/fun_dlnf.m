function [dlnf_dc, d2lnf_dc] = fun_dlnf(cl, P)

cl = cl/1e3;
dlnf_dlnc = (0.601-0.24*cl.^0.5+0.982*(1-0.0052*(P.T-298.15))*cl.^1.5)/(1-0.38)-1; %combined "activity" factor
dlnf_dc = dlnf_dlnc./(cl*1e3);
d2lnf_dlnc =(-0.5*0.24*cl.^(-0.5)+0.982*(1-0.0052*(P.T-298.15))*1.5*cl.^0.5)/(1-0.38);
d2lnf_dc  = -dlnf_dlnc./cl.^2 + d2lnf_dlnc./cl;
d2lnf_dc = d2lnf_dc/1e6;

% dlnf_dc = zeros(size(cl));
% d2lnf_dc = zeros(size(cl));