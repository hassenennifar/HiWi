function [value,isterminal,direction] = myEventsFcn(t,y,p)
% MYEVENTFUNCTION  checks for two events: A steep flank in input current
% and the voltage safe voltage range in which model is valid. This function
% can be given as an input to the integrator ode15s and will cause it to
% terimate if either event occurs.
%
% value - event occurs if value is zero
%
% isterminal - should event cause intergrator to stop?
%
% direction - considers direction from which value approaces zero
%
% Copyright (c) 2017 Julius Zwirner <mailto:julius.zwirner@posteo.de>  and
% TU München. See license.txt for further information.
% July 2017.

V = y(p.precomp.indices.phis(end));

value = double((abs(p.i(t+1e-3)-p.i(t))<p.C) && p.Vmin < V && V < p.Vmax);
isterminal = 1; % yes
direction = 0; % doesn't matter
end
