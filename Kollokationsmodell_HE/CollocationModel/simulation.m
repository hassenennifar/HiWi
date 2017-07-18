% Numerical solution of the Newman model for a Li-Ion cell
function result = simulation(n,m,profile,plots)

% set model parameters
p = parameters;

% precompute spatial discretization and other quantities that need to be calculated only once: 
p.precomp = precompute(n,m,p);

if strcmp(profile,'10C on/off')
    %input current: 10C dischrage/charge for 1s each
    p.i = @(t)10*p.C*((1<t & t<=2)-(6<t & t<=7));         
    p.tend = 20;
elseif strcmp(profile,'1C charge')
    %input current: constant 1C charge
    p.x0 = [0.7 NaN 0.5];
    p.i = @(t)-p.C*(0<=t);          
    p.tend = 100; 
elseif strcmp(profile,'10C discharge')
    %input current: constant 10C discharge for 70s
    p.i = @(t)10*p.C;          
    p.tend = 70;
else
    %input current: constant 1C discharge for 3600s
    p.i = @(t)p.C*(t<3600);          
    p.tend = 5000;
end

% initialize result
tout = []; yout = []; t = 0; y = [];

% set integrator options
options = odeset('Mass',p.precomp.M,'MassSingular','yes','MStateDependence','none','RelTol',1e-9,...
    'Events',@(t,y)myEventsFcn(t,y,p),'Jacobian',@(t,y)jac(t,y,p),'Stats','on');

while t < p.tend
y = initialConditions(t,y,p);
% tic;
[t,y] = ode15s(@(t,y)NewmanDAE(t,y,p),[t,p.tend],y,options);
toc
tout = [tout;t]; yout = [yout;y]; t = tout(end)+1e-3; y = yout(end,:)'; disp(' ');%#ok<AGROW>
end

result = interprete(p,tout,yout);
if ~strcmp(plots,'none')
    batteryplot(result,plots);
end
end