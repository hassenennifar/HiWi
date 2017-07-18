function [act,dact_dcl,dact_dT]=activity(cl,T,mode)
%act = (1+dlnf/dcl)(1-t_plus)*2R/F;
act = 0.601-7.5895*cl.^0.5+0.982*(1-0.0052*(T-298.15))*3.1623e4*cl.^1.5; %combined "activity" factor
%COMSOL: dlnf/dcl=(0.601-0.24*cl^0.5+0.982*(1-0.0052*(T-298.15))*cl^1.5)/(1-t_plus)-1
dact_dcl = [];
dact_dT = [];
if nargin == 3 && strcmp(mode,'partials')
    dact_dcl = 0.5*-7.5895*cl.^-0.5+0.982*(1-0.0052*(T-298.15))*1.5*3.1623e4*cl.^0.5;
    dact_dT = 0.982*-0.0052*T*3.1623e4*cl.^1.5;
end

end