function [y,yp] = initialConditions(t,y,p)

ind = p.precomp.indices;
n = p.precomp.n; m = p.precomp.m;
neg = 1:n(1); pos = n(1)+1:n(1)+n(3);

if isempty(y)
    y=zeros(ind.T,1); yp = y;
    cs = p.x0.*p.cs_max; r = p.precomp.r;
    rcs_neg = repmat(r(:,1)*cs(1),n(1),1); rcs_pos = repmat(r(:,3)*cs(3),n(3),1);
    y(ind.rcs) = [rcs_neg; rcs_pos];                                      
    y(ind.cl) = p.cl_0;
    y(ind.T) = p.T0; 
    y(ind.phis(neg)) = 0;
    y(ind.phis(pos)) =  OCV(cs(3),p.T0,3,p)-OCV(cs(1),p.T0,1,p);  
    
end 
    y(ind.ir(neg)) = p.i(t)/p.d(1);
    y(ind.ir(pos)) = -p.i(t)/p.d(3);
    
    
    y(ind.il(neg)) = p.precomp.x_neg* p.i(t)/p.d(1);
    y(ind.il(pos)) = p.i(t)*(1-(p.precomp.x_pos-p.d(1)-p.d(2))/p.d(3));
    
    T = y(ind.T);
    alpha = p.F/(2*p.R*T); 
    cs = y(ind.rcs(m))/p.rp(1); cl = y(ind.cl);
    i0_neg = p.kprimed(1)*sqrt(cs*(p.cs_max(1)-cs)*cl(1));
    eta_neg = asinh(p.i(t)/p.d(1)/i0_neg)/alpha;
    phil_neg =  - OCV(cs,T,1,p) - eta_neg;
    y(ind.phil) = phil_neg;
    cs = y(ind.rcs(end))/p.rp(3); 
    i0_pos = p.kprimed(3)*sqrt(cs*(p.cs_max(3)-cs)*cl(end));
    eta_pos = asinh(p.i(t)/p.d(1)/i0_pos)/alpha;
    y(ind.phis(pos)) = eta_pos + phil_neg+ OCV(cs,T,3,p);
    yp = NewmanDAE(t,y,p);

end