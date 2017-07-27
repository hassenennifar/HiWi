% run simulation with a selection of current profiles and plotting options

clear;

% spatial discretization in anode, separator, cathode

n = [10 10 10];

% spatial discretization in particle dimension

m = 25;

% Choose input profile from following options: 
 
profile = '1C discharge' 
%profile = '1C charge'
%profile = '10C on/off'
%profile = '10C discharge'
%profile = 'volatile'
%profile = 'artemis'

% Choose output plots (all is slow) 

%plotting = 'none';
%plotting = 'timePlotsOnly';
plotting = 'all';

result=simulation(n,m,profile,plotting);

