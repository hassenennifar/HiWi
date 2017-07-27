function p = parameters()

L_neg = 50e-3;      %[mm] Thickness of graphite electrode
L_sep = 20e-3;      %[mm] Thickness of separator
L_pos = 45e-3;      %[mm] Thickness of positive electrode (NMC)

Ds_neg = 3e-8;      %[mm^2/s] Solid phase Li-diffusivity graphite
Ds_pos = 1e-9;      %[mm^2/s] Solid phase Li-diffusivity NMC
rp_neg = 10e-3;      %[mm] Particle radius graphite
rp_pos = 5e-3;      %[mm] Particle radius NMC
epss_pos = 0.55;      %Solid phase vol-fraction positive electrode
epsl_pos = 0.35;      % Electrolyte phase vol-fraction positive electrode
epss_neg = 0.6;      % Solid phase vol-fraction negative electrode
epsl_neg = 0.35;      % Electrolyte phase vol-fraction negative electrode
epsl_sep = 0.5;      % Porosity separator
cs_max_neg = 28e-3;      %[mmol/mm^3] Max solid phase concentration negative electrode
cs_max_pos = 50e-3;      %[mmol/mm^3] Max solid phase concentration NMC
sigma_neg = 0.100;    %[S/mm] Solid phase conductivity graphite
sigma_pos = 1e-3;      %[S/mm]
k_neg = 2e-8;      %[mm/s] Reaction rate coefficient negative electrode
k_pos = 2e-8;      %[mm/s] Reaction rate coefficient positive electrode
p.t_plus = 0.38;      % Transport number electrolyte
cl_ref = 1e-6;          % [mmol/mm^3] reference concentration
p.C = 13.6e-6;             %[A/mm^2] 1C rate

%thermal parameters
h = 1e-5;                  %[W/(mm^2*K)] convective heat transport
rho = 1.626e-3;                 % [g/mm^3] density
cp = 0.75;                    % [J/gK] specific heat capacity
A2Vratio = 0.253;          % [1/mm] surface area to volume ratio for 18650 round cell
p.Tamb = 298.15;              %ambient temperature
p.cap = rho*cp;
p.hA = h*A2Vratio;

%voltage limits
p.Vmax = 4.2;
p.Vmin = 2.4;

%physical constants
p.R = 8.31446e-3;             % Gas constant [J/(K mmol)]
p.F = 96.4853;               % Faraday constant [As/mmol]

% store quantities for convenience in vectors 
p.d = [L_neg L_sep L_pos];
p.l = sum(p.d);
p.epss = [epss_neg NaN epss_pos];
p.epsl = [epsl_neg epsl_sep epsl_pos];
p.rp = [rp_neg NaN rp_pos];
p.as = 3*p.epss./p.rp;

p.cs_max = [cs_max_neg NaN cs_max_pos];

p.kprimed = 2*p.F*p.as*sqrt(k_neg*k_pos/cl_ref);
p.Ds = [Ds_neg NaN Ds_pos]./p.rp.^2;

%effective values after Bruggeman
brug = 1.5;         
p.effl = p.epsl.^brug;
effs = p.epss.^brug;
p.sigma = [sigma_neg NaN sigma_pos].*effs;

%initial conditions 
p.T0 = 298.15;                %[K] inital Temperature
p.cl_0 = 1e-3;             % [mmol/mm^3] Initial electrolyte salt concentration
x0_neg = 0.8;              %  Initial concentration graphite
x0_pos = 0.4;             % Initial concentration NMC
p.x0 = [x0_neg NaN x0_pos];

%functions
filename = 'Eeq_neg_Graphite.txt';
delimiter = ' ';
data = importdata(filename,delimiter);
x = data(:,1); y = data(:,2);
p.Eeq{1} = splinefit(x,y,12);
p.dEeqdx{1} = ppdiff(p.Eeq{1},1);

filename = 'Eeq_pos_NMC.txt';
data = importdata(filename,delimiter);
x = data(:,1); y = data(:,2);
breaks = [0 0.9  0.975 1];
p.Eeq{3}  = splinefit(x,y,breaks);
p.dEeqdx{3}  = ppdiff(p.Eeq{3},1);

filename = 'dEeqdT_neg_Graphite.txt';
data = importdata(filename,delimiter);
x = data(:,1); y = 1e-3*data(:,2); % units of original data are mV/K - we use V/K
breaks = [0 0.04 0.07 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.48 0.51 0.54 0.6 1];
p.dEeqdT{1}  = splinefit(x,y,breaks);
p.d2EeqdTdx{1}  = ppdiff(p.dEeqdT{1},1);

filename = 'dEeqdT_pos_NMC.txt';
data = importdata(filename,delimiter);
x = data(:,1); y = 1e-3*data(:,2); % units of original data are mV/K - we use V/K 
p.dEeqdT{3}  = splinefit(x,y,8);
p.d2EeqdTdx{3}  = ppdiff(p.dEeqdT{3},1);

end
