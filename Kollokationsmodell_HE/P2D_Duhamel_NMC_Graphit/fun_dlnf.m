function [dlnf_dc, d2lnf_dc, dlnf_dlnc, d2lnf_dlnc] = fun_dlnf(cl, P)

% cl = cl/1e3;
% dlnf_dlnc = (0.601-0.24*cl.^0.5+0.982*(1-0.0052*(P.T-298.15))*cl.^1.5)/(1-0.38)-1; %combined "activity" factor
% dlnf_dc = dlnf_dlnc./(cl*1e3);
% d2lnf_dlnc =(-0.5*0.24*cl.^(-0.5)+0.982*(1-0.0052*(P.T-298.15))*1.5*cl.^0.5)/(1-0.38);
% d2lnf_dc  = -dlnf_dlnc./cl.^2 + d2lnf_dlnc./cl;
% d2lnf_dc = d2lnf_dc/1e6;

% dlnf_dc = zeros(size(cl));
% d2lnf_dc = zeros(size(cl));

%% Julius collocation method
cl = cl*(1e3/1e9); % convert [mol/m3] to [mmol/mm3]

%act = (1+dlnf/dcl)(1-t_plus)*2R/F;
T = P.T;
F = 96.4853;               % Faraday constant [As/mmol]
R = 8.31446e-3;             % Gas constant [J/(K mmol)]
act = 0.601-7.5895*cl.^0.5+0.982*(1-0.0052*(T-298.15))*3.1623e4*cl.^1.5;
dlnf_dlnc = act; % convert to dlnf/dlnc
dlnf_dc = dlnf_dlnc./(cl*1e-3/1e-9);

d2lnf_dlnc = 0.5*-7.5895*cl.^-0.5+0.982*(1-0.0052*(T-298.15))*1.5*3.1623e4*cl.^0.5;
d2lnf_dc  = -dlnf_dlnc./cl.^2 + d2lnf_dlnc./cl;
d2lnf_dc = d2lnf_dc/(1e-3/1e-9)^2; % convert [(.)/(mmol/mm3)^2] to [(mol/m3)^2]