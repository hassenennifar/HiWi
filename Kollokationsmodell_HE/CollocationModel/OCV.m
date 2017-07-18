function [U0,U,dU0dcs,dUdcs,dU0dTdcs] = OCV(cssurf,T,dom,p,mode)

U0 = ppval(p.Eeq{dom},cssurf/p.cs_max(dom));
dU0dT = ppval(p.dEeqdT{dom},cssurf/p.cs_max(dom));
U = U0+T*dU0dT;

if nargin ==5 && strcmp(mode,'partials')
    dU0dcs = ppval(p.dEeqdx{dom},cssurf/p.cs_max(dom))/p.cs_max(dom);
    dU0dTdcs = ppval(p.d2EeqdTdx{dom},cssurf/p.cs_max(dom))/p.cs_max(dom);
    dUdcs = dU0dcs + T*dU0dTdcs;
end

end