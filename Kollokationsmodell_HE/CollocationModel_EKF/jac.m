function dfdy = jac(t,y,p)
% JAC  computes the Jacobian matrix dfdy at time t and state y(t), where 
% f is the right-hand side of the DAE My'=f(t,y) that is computed in 
% NewmanDAE.m Accordingly any change in NewmanDAE.m makes a corresponding 
% change in jac.m necessary. While the computation of the Jacobian is not
% necessary to use ode15s, it is essential for speed.
%
% dfdy - Jacobian matrix of f defined by NewmanDAE.m
%
% Copyright (c) 2017 Julius Zwirner <mailto:julius.zwirner@posteo.de>  and
% TU München. See license.txt for further information.
% July 2017.

n = p.precomp.n; m = p.precomp.m; ind = p.precomp.indices;
dfdy = zeros(p.precomp.d);

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
D = p.precomp.D{k}; w = p.precomp.w{k};

if k~=2
    % indices and components of state vector defined on the electrodes only
    rightEl = leftEl+n(k)-1; dom = leftEl:rightEl; 
    ind_ir = ind.ir(dom); ir =  y(ind_ir);
    ind_il = ind.il(dom);
    ind_phis = ind.phis(dom); phis = y(ind_phis);
    dom = (leftEl-1)*m+1:m*rightEl; leftEl = rightEl + 1;
    ind_cs = ind.cs(dom); rcs = y(ind_cs);

    % derivatives of reaction current   
    cssurf = rcs(m:m:end)/p.rp(k); %surface concentration
    i0 = p.kprimed(k)*sqrt(cssurf.*(p.cs_max(k)-cssurf).*cl);
    [U0,U,dU0dcs,dUdcs,dUdTdcs] = OCV(cssurf,T,k,p,'partials'); 
    %[U0,U] = OCV(cssurf,T,k,p); q = q-w*(ir.*U0); too small can ignore
    
    dfdy(ind.T,ind_ir) = -w.*U0'/(p.l*p.cap);
    dfdy(ind.T,ind_cs(m:m:end)) = -w.*(ir.*dU0dcs)'/(p.l*p.cap*p.rp(k));
    
    eta = phis - phil - U; 
    ir = i0.*sinh(alpha*eta);
    dirdeta = diag(i0.*cosh(alpha*eta)*alpha);
    
    dfdy(ind_ir,ind_ir) = eye(n(k));
    dfdy(ind_ir,ind_phis) = - dirdeta;
    dfdy(ind_ir,ind_phil) = dirdeta;
    dfdy(ind_ir,ind_cl) = -diag(ir./(2*cl)); 
    dfdy(ind_ir,ind_cs(m:m:end)) = (-diag(ir.*(p.cs_max(k)-2*cssurf)./(2*(cssurf.*(p.cs_max(k)-cssurf))))...
        +dirdeta*diag(dUdcs))/p.rp(k);
  

    % derivatives of transformed particle diffusion equation (3.155)
    dfdy(ind_cs,ind_cs) = p.Ds(k)*p.precomp.A{k}; 
    dfdy(ind_cs,ind_ir) = p.precomp.b{k};
    
    % electrolyte current
    if k == 1
        outer = 1; inner = 2:n(1); interface = n(1); 
         dfdy(ind_phil(1),ind_phis(1)) = 1;   %electric ground
    else
        outer = n(3); inner = 1:n(3)-1; interface = 1; 
    end  
    % electrolyte current: il,x = ir; il=0 at x=0,l;
    dfdy(ind_il(outer),ind_il(outer)) = 1;
    dfdy(ind_il(inner),ind_il) = D(inner,:);
    dfdy(ind_il(inner),ind_ir(inner)) = -eye(n(k)-1);

    % electrolyte concentration: (De_eff(ce,T)ce,x),x+(1-t_+)/F*ir=0; ce,x=0 at
    % x=0,l; De_eff*ce,x continuous at electrode-separator boundary
  
    dfdy(ind_cl(outer),ind_cl) = D(outer,:);

    dfdy(ind_cl(interface),ind_cl) = dfdy(ind_cl(interface),ind_cl) + p.effl(k)*D(interface,:); 
    
    %J(ind_cl(2:end-1),ind_T) = Dx.DM{1}(2:end-1,:,1)*(p.eff(1)*p.DeDT(cl,T).*(Dx.DM{1}(:,:,1)*cl));
    
    % solid phase potential
    dfdy(ind_phis(1:end-1),ind_phis) = p.sigma(k)*D(1:end-1,:);
    dfdy(ind_phis(1:end-1),ind_il(1:end-1)) = -eye(n(k)-1);
    dfdy(ind_phis(end),ind_il(interface)) = 1;
    dfdy(ind_phil(2:end),ind_il(2:end)) = eye(n(k)-1);
    dfdy(ind_cl(2:end-1),ind_ir(2:end-1))=(1-p.t_plus)/p.F*eye(n(k)-2)/p.epsl(k);
   
else 
    b = [1 n(k)];
    for j = 1:2
       dfdy(ind_cl(b(j)),ind_cl) = dfdy(ind_cl(b(j)),ind_cl) - p.effl(k)*D(b(j),:);
    end
end

% electrolyte potential

[kap,dkap_dcl,~]=kappa(cl(2:end),T,p.effl(k),'partials'); [act,dact_dcl,~]=activityCoef(cl(2:end),T,'partials');
%kap=kappa(cl(2:end),T,p.effl(k)); act=activity(cl(2:end),T);

dfdy(ind_phil(2:end),ind_phil) = diag(kap)*D(2:end,:);
dfdy(ind_phil(2:end),ind_cl) =  - diag(kap.*act)*(D(2:end,:)*diag(1./cl))/alpha + [zeros(n(k)-1,1) diag( dkap_dcl.*(D(2:end,:)*phil-act.*(D(2:end,:)*log(cl))/alpha) - kap.*dact_dcl.*(D(2:end,:)*log(cl))/alpha)]; 
                         
% electrolyte diffusion  
[Difl,dDifl_dcl,~]= diffusionCoef(cl,T,p.effl(k),'partials');

dfdy(ind_cl(2:end-1),ind_cl) = D(2:end-1,:)*(diag(Difl)*D+diag(dDifl_dcl.*(D(:,:)*cl)))/p.epsl(k);

end

% bulk temperature
dfdy(ind.T,ind.phis(end)) = -13.6e-5/(p.l*p.cap);
dfdy(ind.T,ind.T) = -p.hA/p.cap;

dfdy = sparse(dfdy); % crucial for speed

end