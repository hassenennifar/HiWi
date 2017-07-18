clear;

%discretization parameters
n = [4 2 5]; %number of nodes in anode, separator, cathode
m = 25; %number of nodes in particle dimension

%Choose input profile from following options: 
% '1C discharge', 
% '1C charge',
% '10C on/off'
% '10 discharge'

profile = '1C discharge';

%Choose which plots to display
% 'only time plots', 
% 'all plots',
% 'none'

plots = 'all plots';
tic
result = simulation(n,m,profile,plots);