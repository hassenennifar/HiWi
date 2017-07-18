function result = simulation(n,m,profile,plotting)
% SIMULATION  calculates the Newman model for a Li-Ion cell. n,m describe
% the spatial discretization, profile specifies the input current i(t) and
% plotting allows for a selection of plot options.
% 
% result - structure containing all the calculated quantities
%
% Copyright (c) 2017 Julius Zwirner <mailto:julius.zwirner@posteo.de>  and
% TU München. See license.txt for further information.
% July 2017.


% set model parameters

p = parameters;

% precompute quantities that need to be calculated only once 

p.precomp = precompute(n,m,p);

% choose input current profile

p = current(profile,p);



% set integrator options

options = odeset('Mass',p.precomp.M,...
                 'MassSingular','yes',...
                  'MStateDependence','none',...
                  'AbsTol',1e-4,'RelTol',1e-5,...
                  'Events',@(t,y)myEventsFcn(t,y,p),...
                  'Jacobian',@(t,y)jac(t,y,p),...
                  'Stats','off');

% initialize result

tout = []; yout = []; t = 0; y = []; tjump = [];  

tic;

while t < p.tf

% estimate initial conditions for current y and i(t), if y is still empty
% we estimate initial conditions for initial parameters 

y = initialConditions(t,y,p);

% we integrate the DAE with our options

[t,y] = ode15s(@(t,y)NewmanDAE(t,y,p),[t,p.tf],y,options);

% if the integration finishes before the end time this can be because of 
% two reasons given in my EventFunction. Either we have a very steep flank
% in the input current which ode15s can't handle. In this case we want to 
% append our solution and start the integration afresh

tout = [tout;t]; yout = [yout;y]; %#ok<AGROW>
t = tout(end)+1e-3; y = yout(end,:)';  %#ok<AGROW> 
tjump = [tjump; tout(end)];  

% Or the voltage leaves the safe operating range in which our model is 
% valid. In this case we end the integration. 

V = y(p.precomp.indices.phis(end));

if  V < p.Vmin || V > p.Vmax; p.tf= t; end  

%disp(' '); %if 'Stats' are on this aids readability

end

toc;

% utlility function to perform unit conversions etc. 
result = interprete(p,tout,tjump,yout);

% plot results according to 'plotting'
if ~strcmp(plotting,'none'), batteryplot(result,plotting);end

end