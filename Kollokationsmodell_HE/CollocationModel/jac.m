function dfdy = jac(t,y,p)
% analytical compution of Jacobian for Newman DAE
n = p.precomp.n; m = p.precomp.m; ind = p.precomp.indices;
dfdy = zeros(p.precomp.nvars);

T = y(ind.T);

    
alpha = p.F/(2*p.R*T);

%loop through domains
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
    ind_ir = ind.ir(dom); ir =  y(ind_ir);
    ind_il = ind.il(dom);
    ind_phis = ind.phis(dom); phis = y(ind_phis);
    dom = (leftEl-1)*m+1:m*rightEl; leftEl = rightEl + 1;
    ind_rcs = ind.rcs(dom); rcs = y(ind_rcs);

    % reaction current   
    cssurf = rcs(m:m:end)/p.rp(k); %surface concentration
    i0 = p.kprimed(k)*sqrt(cssurf.*(p.cs_max(k)-cssurf).*cl);
    [U0,U,dU0dcs,dUdcs,dUdTdcs] = OCV(cssurf,T,k,p,'partials');  %[U0,U] = OCV(cssurf,T,k,p); q = q-w*(ir.*U0);
    
    dfdy(ind.T,ind_ir) = -w.*U0'/(p.l*p.cap);
    dfdy(ind.T,ind_rcs(m:m:end)) = -w.*(ir.*dU0dcs)'/(p.l*p.cap*p.rp(k));
    
    eta = phis - phil - U; 
    ir = i0.*sinh(alpha*eta);
    dirdeta = diag(i0.*cosh(alpha*eta)*alpha);
    
    dfdy(ind_ir,ind_ir) = eye(n(k));
    dfdy(ind_ir,ind_phis) = - dirdeta;
    dfdy(ind_ir,ind_phil) = dirdeta;
    dfdy(ind_ir,ind_cl) = -diag(ir./(2*cl)); %-diag(ir./(2*cl));%
    dfdy(ind_ir,ind_rcs(m:m:end)) = (-diag(ir.*(p.cs_max(k)-2*cssurf)./(2*(cssurf.*(p.cs_max(k)-cssurf))))+dirdeta*diag(dUdcs))/p.rp(k);
  

    % solid-phase concentration times radial distance: rcs,t = Ds(T)/R^2*rcs,rr, Ds(T)/R^2(rcs(r=1) -rcs) + ir = 0; 
    dfdy(ind_rcs,ind_rcs) = p.Ds(k)*p.precomp.As{k}; dfdy(ind_rcs,ind_ir) = p.precomp.Bs{k};
    
    % electrolyte current
    if k == 1
        outer = 1; inner = 2:n(1); interface = n(1); 
         dfdy(ind_phil(1),ind_phis(1)) = 1;   %electric ground
    else
        outer = n(3); inner = 1:n(3)-1; interface = 1; 
    end  
    % electrolyte current: il,x = ir; il=0 at x=0,l;
    dfdy(ind_il(outer),ind_il(outer)) = 1;
    dfdy(ind_il(inner),ind_il) = DM(inner,:);
    dfdy(ind_il(inner),ind_ir(inner)) = -eye(n(k)-1);

    % electrolyte concentration: (De_eff(ce,T)ce,x),x+(1-t_+)/F*ir=0; ce,x=0 at
    % x=0,l; De_eff*ce,x continuous at electrode-separator boundary
  
    dfdy(ind_cl(outer),ind_cl) = DM(outer,:);

    dfdy(ind_cl(interface),ind_cl) = dfdy(ind_cl(interface),ind_cl) + p.effl(k)*DM(interface,:); 
    
    %J(ind_cl(2:end-1),ind_T) = Dx.DM{1}(2:end-1,:,1)*(p.eff(1)*p.DeDT(cl,T).*(Dx.DM{1}(:,:,1)*cl));
    
    % solid phase potential
    dfdy(ind_phis(1:end-1),ind_phis) = p.sigma(k)*DM(1:end-1,:);
    dfdy(ind_phis(1:end-1),ind_il(1:end-1)) = -eye(n(k)-1);
    dfdy(ind_phis(end),ind_il(interface)) = 1;
    dfdy(ind_phil(2:end),ind_il(2:end)) = eye(n(k)-1);
    dfdy(ind_cl(2:end-1),ind_ir(2:end-1))=(1-p.t_plus)/p.F*eye(n(k)-2)/p.epsl(k);
   
else 
    b = [1 n(k)];
    for j = 1:2
       dfdy(ind_cl(b(j)),ind_cl) = dfdy(ind_cl(b(j)),ind_cl) - p.effl(k)*DM(b(j),:);
    end
end

% electrolyte potential

[kap,dkap_dcl,~]=kappa(cl(2:end),T,p.effl(k),'partials'); [act,dact_dcl,~]=activity(cl(2:end),T,'partials');
%kap=kappa(cl(2:end),T,p.effl(k)); act=activity(cl(2:end),T);

dfdy(ind_phil(2:end),ind_phil) = diag(kap)*DM(2:end,:);
dfdy(ind_phil(2:end),ind_cl) =  - diag(kap.*act)*(DM(2:end,:)*diag(1./cl))/alpha + [zeros(n(k)-1,1) diag( dkap_dcl.*(DM(2:end,:)*phil-act.*(DM(2:end,:)*log(cl))/alpha) - kap.*dact_dcl.*(DM(2:end,:)*log(cl))/alpha)]; 
                         
% electrolyte diffusion  
[D,dD_dcl,~]=Dl(cl,T,p.effl(k),'partials');

dfdy(ind_cl(2:end-1),ind_cl) = DM(2:end-1,:)*(diag(D)*DM+diag(dD_dcl.*(DM(:,:)*cl)))/p.epsl(k);

end

% bulk temperature
dfdy(ind.T,ind.phis(end)) = -p.i(t)/(p.l*p.cap);
dfdy(ind.T,ind.T) = -p.hA/p.cap;

dfdy = sparse(dfdy);

end