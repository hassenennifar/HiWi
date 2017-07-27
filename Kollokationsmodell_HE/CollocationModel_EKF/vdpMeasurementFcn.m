function yk = vdpMeasurementFcn(xk,vk,p)
% vdpMeasurementNonAdditiveNoiseFcn Example measurement function for discrete
% time nonlinear state estimators with non-additive measurement noise.
%
% yk = vdpNonAdditiveMeasurementFcn(xk,vk)
%
% Inputs:
%    xk - x[k], states at time k
%    vk - v[k], measurement noise at time k
%
% Outputs:
%    yk - y[k], measurements at time k
%
% The measurement is the first state with multiplicative noise
%
% See also extendedKalmanFilter, unscentedKalmanFilter

%   Copyright 2016 The MathWorks, Inc.

%#codegen

% The tag %#codegen must be included if you wish to generate code with 
% MATLAB Coder.

n = p.precomp.n; m = p.precomp.m; ind = p.precomp.indices;
algebr_ind = length(xk);
T = vk(ind.T-algebr_ind);
alpha = p.F/(2*p.R*T);
V=zeros(2,1);
%loop through domains
left = 1; leftEl = 1;


for k = 1:3
    % indices and components of state vector defined on all domains
    right = left+n(k)-1; dom = left:right; left = right;

    ind_cl = ind.cl(dom); cl = xk(ind_cl); 
    ind_phil = ind.phil(dom)-algebr_ind; phil= vk(ind_phil);

    % first derivative matrix and Clenshaw-Curtis weights
    D = p.precomp.D{k}; w = p.precomp.w{k};

    if k~=2 % electrodes

        % indices and components of state vector defined in electrodes only

        rightEl = leftEl+n(k)-1; dom = leftEl:rightEl; 
        ind_ir = ind.ir(dom)-algebr_ind; ir = vk(ind_ir); 
        ind_il = ind.il(dom)-algebr_ind; il = vk(ind_il); 
        ind_phis = ind.phis(dom)-algebr_ind; phis = vk(ind_phis);
        dom = (leftEl-1)*m+1:m*rightEl; leftEl = rightEl + 1;
        ind_cs = ind.cs(dom); rcs = xk(ind_cs); % work with r*cs instead of cs
        
        cssurf = rcs(m:m:end)/p.rp(k); %surface concentration
        i0 = p.kprimed(k)*sqrt(cssurf.*(p.cs_max(k)-cssurf).*cl);
        [U0,U] = OCV(cssurf,T,k,p);
        eta = asinh(ir./i0)/alpha;
        PHIS = eta + phil + U0;
        if k == 1
            V(1) = PHIS(1);
        else
            V(2) = PHIS(end);
        end
    end
end
yk = V(2)-V(1);
end


