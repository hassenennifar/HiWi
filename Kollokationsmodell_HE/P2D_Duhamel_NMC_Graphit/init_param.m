function P = init_param()

P.F = 96485;  % Faraday constant (As/mol)
P.R = 8.314;  % gas constant (J/mol-K)

P.idx_cl = 1; % variable index liquid concentration 
P.idx_phil = 2; % variable index liquid potential
P.idx_css = 3; % variable index solid surface concentration
P.idx_il = 4; % variable index liquid current density
P.idx_jn = 5; % variable index pore-wall molar flux
P.idx_phis = 6; % variable index solid phase potential

P.nx = 6; % number of states
P.nc = 2; % number of control variables

P.lim = 1000;  % limit on number of iterations
P.tolrel = 1e-4; % relative tolerance
P.tolabs = 1e-10; % abolute tolerance

P.n_neg = 3; % number of nodes in negative electrode
P.n_sep = 1; % number of nodes in separator
P.n_pos = 4; % number of nodes in positive electrode
P.nj = P.n_neg+P.n_sep+P.n_pos+1; % total number of nodes

P.bnd_sep_neg = P.n_neg+1; % boundary node between separator and anode
P.bnd_pos_sep = P.n_neg+P.n_sep+1; % boundary node between cathode and separator

P.L_neg = 50e-6;  % thickness of negative electrode (m)
P.L_sep = 20e-6;  % thickness of separator (m)
P.L_pos = 45e-6;  % thickness of positive electrode (m)

P.hcn = 0e-06;  % thickness of negative electrode current collector (m)
P.hcp = 0e-06;  % thickness of positive electrode current collector (m)

P.cl0 = 1000;  % initial electrolyte concentration (mol/m3)
P.sto_neg = 0.8;  % initial stochiometric parameter for negative particle
P.sto_pos = 0.4;  % initial stochiometric parameter for positive particle
P.dtmax = 1.02;  % maximum time step size (s)

P.Ds_neg = 3e-14;  % diffusion coefficient in negative solid (m2/s)
P.Ds_pos = 1e-15;  % diffusion coefficient in positive solid (m2/s)

P.Ds_neg0 = P.Ds_neg; % save value at 25degC for Arrhenius coefficient
P.Ds_pos0 = P.Ds_pos; % save value at 25degC for Arrhenius coefficient

P.Rp_neg = 10e-6;  % radius of negative particles (m)
P.Rp_pos = 5e-6;  % radius of positive particles (m)

P.epsl_neg = 0.35;  %  volume fraction of electrolyte in negative electrode
P.epsl_sep = 0.5;  %  volume fraction of electrolyte in separator
P.epsl_pos = 0.35;  %  volume fraction of electrolyte in positive electrode

P.epf_neg = 0.05;  %  volume fraction of inert filler in negative electrode
P.epf_pos = 0.10;  %  volume fraction of inert filler in positive electrode

P.epss_neg = 1-P.epsl_neg-P.epf_neg; % volume fraction of solid phase in negative electrode
P.epss_pos = 1-P.epsl_pos-P.epf_pos; % volume fraction of solid phase in negative electrode

P.sigma_neg = 100;  %  conductivity of solid negative matrix (S/m)
P.sigma_pos = 1;  %  conductivity of solid positive matrix (S/m)

P.rka_neg = 2e-11;  %  reaction rate constant for negative insertion reaction
P.rka_pos = 2e-11;  %  reaction rate constant for positive insertion reaction

P.rka_neg0 = P.rka_neg; % save value at 25degC for Arrhenius coefficient
P.rka_pos0 = P.rka_pos; % save value at 25degC for Arrhenius coefficient

P.Rf_neg = 0;  %  anode film resistance (out of place) (ohm-m2) 
P.Rf_pos = 0;  %  cathode film resistance (out of place) (ohm-m2)

P.Rf_neg0 = P.Rf_neg; % save value at 25degC for Arrhenius coefficient
P.Rf_pos0 = P.Rf_pos; % save value at 25degC for Arrhenius coefficient

P.cot_neg = 372;  %  coulombic capacity of negative material (mAh/g)
P.cot_pos = 274;  %  coulombic capacity of positive material (mAh/g)

P.re = 1324;  %  density of electrolyte (kg/m3)
P.rs_neg = 1800;  %  density of negative insertion material (kg/m3)
P.rs_pos = 5010;  %  density of positive insertion material (kg/m3)
P.rf = 1800;  %  density of inert filler (kg/m3)
P.rpl = 1780;  %  density of polymer phase (kg/m3)
P.rc = 552;  %  density of separator material (kg/m3)
P.rcn = 8954;  %  density of negative current collector (kg/m3)
P.rcp = 2707;  %  density of positive current collector (kg/m3)

P.csmax_neg = 28e3;%3600*P.cot_neg*P.rs_neg/P.F; % max solid phase concentration negative electrode (mol/m3)
P.csmax_pos = 50e3;%3600*P.cot_pos*P.rs_pos/P.F; % max solid phase concentration positive electrode (mol/m3)
P.cs0_neg = P.sto_neg*P.csmax_neg; % init solid phase concentration negative electrode (mol/m3)
P.cs0_pos = P.sto_pos*P.csmax_pos; % init solid phase concentration positive electrode (mol/m3)

P.htc = .10;  % heat-transfer coefficient with external medium (W/m2K)
P.Cp = 750;  % heat capacity of cell (J/kgK)
P.Tam  = 298.15;  % ambient temperature (K)
P.T  = P.Tam;  % temperature (K)
P.A2Vratio = 253; % [1/m] surface area to volume ratio for 18650 round cell
P.rho = 1626;
P.cap = P.rho*P.Cp;
P.hA = P.htc*P.A2Vratio;

P.RG = 8e-4*0; % total resistance in foils, leads, and contacts, ohm-m2
P.RGext = P.RG/4; % resistance outside cell

P.nLi_solid_matlab=P.epss_neg*P.cs0_neg*P.L_neg+P.epss_pos*P.cs0_pos*P.L_pos;
P.totLiold = (P.epsl_neg*P.L_neg+P.epsl_sep*P.L_sep+P.epsl_pos*P.L_pos)*P.cl0;

P.EbarDs_neg = 0;  % activation energy, anode solid diffusion
P.EbarDs_pos = 0;  % activiation energy, cathode solid diffusion
P.Ebarkap = 4000;  % activation energy electrolyte conductivity
P.EbarD  = 4000;  % activation energy electrolyte diffusion
P.Ebarka_neg = 4000;  % activation energy negative kinetics
P.Ebarka_pos = 4000;  % activation energy positive kinetics
P.Ebarks1a = 4000;  % activation energy O2 side rxn
P.Ebarks1c = 4000;  % activation energy O2 side rxn
P.Ebarks2a = 4000;  % activation energy H2 side rxn
P.Ebarks2c = 4000;  % activation energy H2 side rxn
P.Ebarks3a = 4000;  % activation energy shuttle side rxn
P.Ebarks3c = 4000;  % activation energy shuttle side rxn
P.Ebarr_neg = 4000;  % activation energy, film resistance anode
P.Ebarr_pos = 4000;  % activation energy, film resistance cathode

