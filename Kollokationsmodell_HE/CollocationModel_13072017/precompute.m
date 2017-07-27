function precomp = precompute(n,m,p)
% PRECOMPUTE  performs all the calculations and in particular the spatial
% discretization that needs to be performed only once per simulation and
% stores it in the structure precomp
%
% precomp - data structure containing precalculated values
%
% Copyright (c) 2017 Julius Zwirner <mailto:julius.zwirner@posteo.de>  and
% TU München. See license.txt for further information.
% July 2017.

N = sum(n)-2;                       % nr of nodes along x-axis
Nel = n(1) + n(3);                  % nr of nodes discretizing electrodes
precomp.n = n; precomp.m = m;

% left/right domain boundaries with/out separator

indices.left = [1 n(1) n(1)+n(2)-1];              
indices.right = [n(1) n(1)+n(2)-1 N];             
indices.leftElectrode = [1 NaN n(1)+1];           
indices.rightElectrode = [n(1) NaN n(1)+n(3)];    

% specify components of state vector

indices.cl   = 1:N;                                % length: N
indices.cs   = N+1:N+m*Nel;                         % length: m*Nel;
indices.ir   = N+m*Nel+1:N+(1+m)*Nel;                 % length: Nel
indices.il   = N+(1+m)*Nel+1:N+(2+m)*Nel;           % length: Nel
indices.phis = N+(2+m)*Nel+1:N+(3+m)*Nel;           % length: Nel
indices.phil = N+(3+m)*Nel+1:2*N+(3+m)*Nel;         % length: N
indices.T    = 2*N+(3+m)*Nel+1;                     % length: 1

d = indices.T;                                %length of state vector
precomp.indices = indices; precomp.d = d;   %store result

% store total dxd sparse diagonal mass matrix (we have 1s for the equation 
% containing temporal derivatives, 0s for algebraic equations)

Mcl = [0; ones(n(1)-2,1); 0; ones(n(2)-2,1); 0; ones(n(3)-2,1); 0];
Mir = zeros(Nel,1); Mil = zeros(Nel,1); Mphis = zeros(Nel,1);
Mphil = zeros(N,1);
Mcs = repmat([ones(m-1,1); 0],Nel,1);
MT = 1;
precomp.M = spdiags([Mcl; Mcs; Mir; Mil; Mphis; Mphil; MT],0,d,d);

% compute Chebychev differentiation matrices and Clenshaw-Curtis weights
% for x discretization in the three domains

xi = cell(3,1); D = cell(3,1); w = cell(3,1);

for k = 1:3 
    [D{k},xi{k},w{k}] = cheb(n(k)); 
    D{k} = 2/p.d(k)*D{k}; w{k} = w{k}*p.d(k)/2; % rescale 
end

precomp.D = D; precomp.w = w; % store result

% map domain discretizations onto model geometry 

precomp.x_neg = p.d(1)*(xi{1}+1)/2; 
precomp.x_sep = p.d(1)+ p.d(2)*(xi{2}+1)/2; 
precomp.x_pos = p.d(1)+p.d(2)+ p.d(3)*(xi{3}+1)/2;  
precomp.x =  [precomp.x_neg; precomp.x_sep(2:end-1); precomp.x_pos];   
                         
% compute differentiation matrix for particle dimension and matrices 
% A, b for solid-state diffusion equation in the form 
% M du/dt = Ds/Rs^2*A*u + b*ir (see paper for details)

M = 2*m+1; [D,rho] = cheb(M); 
precomp.r = rho(m+2:M)*p.rp;                    % nodes in the r-dimension
D2 = D^2;                                       % second derivative
D = D(M,m+2:M)- D(M,m:-1:1);                    % use symmetry
D2 = D2(m+2:M-1,m+2:M) - D2(m+2:M-1,m:-1:1);    % use symmetry
A = [D2; (D - fliplr(eye(1,m)))];               % last row is due to BC
b = flipud(eye(m,1))/p.F;                       % [0  ...  0 1]' for BC

% diagonalise first m-1 rows and collumns of A. The eigenvectors should
% always be real, but I haven't proved this. If not we simply don't 
% diagonalise.

[V,~] = eig(D2(:,1:end-1)); %eigenvectors
V = blkdiag(V,1); % append 1 along diagonal
if ~isreal(V), V = eye(m); end
precomp.V = V;

A = V\A*V;          % is diagonal except for last row and column
A(abs(A)<1e-9)=0;   % clean rounding errors

         
% precomp.wr = wr(m+2:M+1)';

% currently not needed
%C = -(D(M+1,m+2:M+1)+D(m+1,m:-1:1))/D(M+1,m+1);
%D = - 1/D(M+1,m+1);

% finally replicate results for every node in electrodes and store result  

precomp.A{1} = kron(eye(n(1)),A); 
precomp.A{3} = kron(eye(n(3)),A);

bReps = repmat({b},1,max(n(1),n(3)));
precomp.b{1} = blkdiag(bReps{1:n(1)})/p.as(1);
precomp.b{3} = blkdiag(bReps{1:n(3)})/p.as(3);

end