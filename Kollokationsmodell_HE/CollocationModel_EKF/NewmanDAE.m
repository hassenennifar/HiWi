function f = NewmanDAE(t,y,p)
% NEWMANDAE  calculates the right-hand side f of the spatially discretized
% Newman model in the form My'=f(t,y) at time t and state y(t). 
% 
% f - right-hand side of spatially discretized Newman model
%
% Copyright (c) 2017 Julius Zwirner <mailto:julius.zwirner@posteo.de>  and
% TU München. See license.txt for further information.
% July 2017.

n = p.precomp.n; m = p.precomp.m; ind = p.precomp.indices;
f = zeros(size(y));

T = y(ind.T);
alpha = p.F/(2*p.R*T);

q=0;

%loop through domains
left = 1; leftEl = 1;
for k = 1:3

% indices and components of state vector defined on all domains
right = left+n(k)-1; dom = left:right; left = right;

ind_cl = ind.cl(dom); cl = y(ind_cl); 
ind_phil = ind.phil(dom); phil= y(ind_phil);

% first derivative matrix and Clenshaw-Curtis weights
D = p.precomp.D{k}; w = p.precomp.w{k};

if k~=2 % electrodes
    
    % indices and components of state vector defined in electrodes only
    
    rightEl = leftEl+n(k)-1; dom = leftEl:rightEl; 
    ind_ir = ind.ir(dom); ir = y(ind_ir); 
    ind_il = ind.il(dom); il = y(ind_il); 
    ind_phis = ind.phis(dom); phis = y(ind_phis);
    dom = (leftEl-1)*m+1:m*rightEl; leftEl = rightEl + 1;
    ind_cs = ind.cs(dom); rcs = y(ind_cs); % work with r*cs instead of cs

    % reaction current  (Gleichung 2.5 in paper) 
    
    cssurf = rcs(m:m:end)/p.rp(k); %surface concentration
    i0 = p.kprimed(k)*sqrt(cssurf.*(p.cs_max(k)-cssurf).*cl);
    [U0,U] = OCV(cssurf,T,k,p);
    eta = phis - phil - U;                      % (2.6)
    f(ind_ir)  = ir - i0.*sinh(alpha*eta);      % (2.5)
    
    q = q-w*(ir.*U0); % heat term 
    
    % transformed particle diffusion equation (3.155)
    f(ind_cs) = p.Ds(k)*p.precomp.A{k}*rcs + p.precomp.b{k}*ir;

    % electrolyte current
    
    if k == 1
        outer = 1; inner = 2:n(1); interface = n(1); 
        f(ind_phil(1)) = phis(1);                   % electric ground BC 
    else
        outer = n(3); inner = 1:n(3)-1; interface = 1; 
    end  
    
    f(ind_il(outer)) = il(outer);                     % no il current BC
    
    f(ind_il(inner)) = D(inner,:)*il-ir(inner); % electrode-separator BC
    
    f(ind_cl(outer)) = D(outer,:)*cl;   % no concentration flow outer BC  
    
    f(ind_cl(interface)) = f(ind_cl(interface)) + ...
         p.effl(k)*D(interface,:)*cl;  % continued electrode-separator BC 
    
    % ohms law (2.8)
    f(ind_phis(1:end-1)) = p.sigma(k)*(D(1:end-1,:)*phis)+13.6e-5-il(1:end-1);
    f(ind_phis(end)) = il(interface)-13.6e-5;        % is = 0 in separator
   
else %separator
    
    il = ones(n(k),1)*13.6e-5; ir = zeros(n(k),1);
    
    % electrode-separator BC for
    b = [1 n(k)];
    for j = 1:2
        f(ind_cl(b(j))) = f(ind_cl(b(j))) - p.effl(k)*D(b(j),:)*cl;
    end
end

% electrolyte potential
f(ind_phil(2:end)) = kappa(cl(2:end),T,p.effl(k)).*(D(2:end,:)*phil - ... 
    activityCoef(cl(2:end),T)/alpha.*(D(2:end,:)*log(cl)))+il(2:end);

% electrolyte diffusion (2.16)
f(ind_cl(2:end-1)) = (D(2:end-1,:)*(diffusionCoef(cl,T,p.effl(k)).*(D*cl))+ ...
    (1-p.t_plus)/p.F*ir(2:end-1))/p.epsl(k);

end

% bulk temperature model (2.17)

f(ind.T) = ((-phis(end)*13.6e-5+q)/p.l-p.hA*(T-p.Tamb))/p.cap;

end