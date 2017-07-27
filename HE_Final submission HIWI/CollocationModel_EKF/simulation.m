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

tout = []; yout = []; t = 0; y = []; tjump = [];  dt=2;
tic;
% estimate initial conditions for current y and i(t), if y is still empty
% we estimate initial conditions for initial parameters 

y = initialConditions(t,y,p);

R = 1e-6;

mat = load('results_2s.mat');
yTrue = mat.result.y;
rng(1); % Fix the random number generator for reproducible results
yMeas = yTrue .*(1+sqrt(R)*randn(size(yTrue))*0); % sqrt(R): Standard deviation of noise

Nsteps = size(yMeas,1); % Number of time-steps
Nstates = size(yMeas,2);
yCorrected = zeros(Nsteps,Nstates); % Corrected state estimates
PCorrected = zeros(Nsteps,Nstates,Nstates); % Corrected state estimation error covariances
e = zeros(Nsteps,1); % Residuals (or innovations)


% Your initial state guess at time k, utilizing measurements up to time k-1: xhat[k|k-1]
initialStateGuess = y; % xhat[k|k-1]
dynamic_ind = [p.precomp.indices.cl, p.precomp.indices.cs];
algebr_ind = 1:Nstates; algebr_ind(dynamic_ind) = [];
xk = y(dynamic_ind);
vk = y(algebr_ind);
% Construct the filter
ekf = extendedKalmanFilter(...
    @vdpStateFcn,... % State transition function
    @vdpMeasurementFcn,... % Measurement function
    xk,...
    'HasAdditiveMeasurementNoise',false);
ekf.MeasurementNoise = R;
indices = [p.precomp.indices.cl, p.precomp.indices.cs];
M = p.precomp.M; Mt = full(M(indices, indices)); mt = diag(Mt);
% indices = find(mt);
Mt = Mt(indices, indices);
ekf.ProcessNoise = Mt*1e2;

% timevector = t:dt:p.tf;

for k=1:Nsteps
    % Let k denote the current time
    %
    % Residuals (or innovations): Measured output - Predicted output
    e(k) = yMeas(k,p.precomp.indices.phis(end)) - vdpMeasurementFcn(ekf.State(dynamic_ind),vk,p); % ekf.State is x[k|k-1] at this point
    % Incorporate the measurements at time k into the state estimates by
    % using the "correct" command. This updates the State and StateCovariance
    % properties of the filter to contain x[k|k] and P[k|k]. These values
    % are also produced as the output of the "correct" command.
    [yCorrected(k,:), PCorrected(k,:,:)] = correct(ekf,yMeas(k,p.precomp.indices.phis(end)));
    % Predict the states at next time step, k+1. This updates the State and
    % StateCovariance properties of the filter to contain x[k+1|k] and
    % P[k+1|k]. These will be utilized by the filter at the next time-step
    predict(ekf,p);
    
    %% 
    % fsolve for algebraic states (consistent solution)
end

% while t < p.tf
% 
% % estimate initial conditions for current y and i(t), if y is still empty
% % we estimate initial conditions for initial parameters 
% 
% y = initialConditions(t,y,p);
% 
% % we integrate the DAE with our options
% timevector = t:dt:p.tf;
% [t,y] = ode15s(@(t,y)NewmanDAE(t,y,p),timevector,y,options);
% 
% % if the integration finishes before the end time this can be because of 
% % two reasons given in my EventFunction. Either we have a very steep flank
% % in the input current which ode15s can't handle. In this case we want to 
% % append our solution and start the integration afresh
% 
% tout = [tout;t]; yout = [yout;y]; % #ok<AGROW>
% t = tout(end)+2; y = yout(end,:)';  % #ok<AGROW> 
% tjump = [tjump; tout(end)];  
% 
% % Or the voltage leaves the safe operating range in which our model is 
% % valid. In this case we end the integration. 
% 
% V = y(p.precomp.indices.phis(end));
% 
% if  V < p.Vmin || V > p.Vmax; p.tf= t; end  
% 
% %disp(' '); %if 'Stats' are on this aids readability
% 
% end

toc;

% utlility function to perform unit conversions etc. 
% result = interprete(p,tout,tjump,yout);
% 
% % plot results according to 'plotting'
% if ~strcmp(plotting,'none'), batteryplot(result,plotting);end
% 
% end