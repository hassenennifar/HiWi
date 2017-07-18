function f = NewmanDAE(t,y,p)
%this function implements the space discretized Newman equations
n = p.precomp.n; m = p.precomp.m; ind = p.precomp.indices;
f = zeros(size(y));

T = y(ind.T);
alpha = p.F/(2*p.R*T);

q=0;

%loop through domains: anode, separator, cathode
left = 1; leftEl = 1;
for k = 1:3

% indices and components of state vector defined on all domains
right = left+n(k)-1; dom = left:right; left = right;

ind_cl = ind.cl(dom); cl = y(ind_cl); 
ind_phil = ind.phil(dom); phil= y(ind_phil);

% first derivative matrix and Clenshaw-Curtis weights
DM = p.precomp.DM{k}; w = p.precomp.w{k};

if k~=2
    % indices and components of state vector defined on the electrodes only
    rightEl = leftEl+n(k)-1; dom = leftEl:rightEl; 
    ind_ir = ind.ir(dom); ir = y(ind_ir); 
    ind_il = ind.il(dom); il = y(ind_il); 
    ind_phis = ind.phis(dom); phis = y(ind_phis);
    dom = (leftEl-1)*m+1:m*rightEl; leftEl = rightEl + 1;
    ind_rcs = ind.rcs(dom); rcs = y(ind_rcs);

    % reaction current   
    cssurf = rcs(m:m:end)/p.rp(k); %surface concentration
    i0 = p.kprimed(k)*sqrt(cssurf.*(p.cs_max(k)-cssurf).*cl);
    [U0,U] = OCV(cssurf,T,k,p);
    eta = phis - phil - U; 
    f(ind_ir)  = ir - i0.*sinh(alpha*eta); %Butler-Volmer equation
    q = q-w*(ir.*U0);
    
    % solid-phase concentration times radial distance: rcs,t = Ds(T)/R^2*rcs,rr, Ds(T)/R^2(rcs(r=1) -rcs) + ir = 0; 
    f(ind_rcs) = p.Ds(k)*p.precomp.As{k}*rcs + p.precomp.Bs{k}*ir;

    % electrolyte current
    if k == 1
        outer = 1; inner = 2:n(1); interface = n(1); 
        f(ind_phil(1)) = phis(1);   %electric ground
    else
        outer = n(3); inner = 1:n(3)-1; interface = 1; 
    end  
    f(ind_il(outer)) = il(outer);
    f(ind_il(inner)) = DM(inner,:)*il-ir(inner);
    f(ind_cl(outer)) = DM(outer,:)*cl;
    f(ind_cl(interface)) = f(ind_cl(interface)) + p.effl(k)*DM(interface,:)*cl;
    % solid phase potential
    f(ind_phis(1:end-1)) = p.sigma(k)*(DM(1:end-1,:)*phis)+p.i(t)-il(1:end-1);
    f(ind_phis(end)) = il(interface)-p.i(t);
   
else 
    il = ones(n(k),1)*p.i(t); ir = zeros(n(k),1);
    b = [1 n(k)];
    for j = 1:2
        f(ind_cl(b(j))) = f(ind_cl(b(j))) - p.effl(k)*DM(b(j),:)*cl;
    end
end

% electrolyte potential
f(ind_phil(2:end)) = kappa(cl(2:end),T,p.effl(k)).*(DM(2:end,:)*phil-activity(cl(2:end),T)/alpha.*(DM(2:end,:)*log(cl)))+il(2:end);

% electrolyte diffusion
f(ind_cl(2:end-1)) = (DM(2:end-1,:)*(Dl(cl,T,p.effl(k)).*(DM*cl))+(1-p.t_plus)/p.F*ir(2:end-1))/p.epsl(k);

end

% bulk temperature model

f(ind.T) = ((-phis(end)*p.i(t)+q)/p.l-p.hA*(T-p.Tamb))/p.cap;

end