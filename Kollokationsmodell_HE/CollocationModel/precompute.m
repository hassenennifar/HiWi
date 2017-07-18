function precomp = precompute(n,m,p)

N = sum(n)-2;                       % nr of nodes along x-axis
Nel = n(1) + n(3);                  % nr of nodes discretizing the electrodes
precomp.n = n; precomp.m = m;

indices.left = [1 n(1) n(1)+n(2)-1];
indices.right = [n(1) n(1)+n(2)-1 N];
indices.leftElectrode = [1 NaN n(1)+1];
indices.rightElectrode = [n(1) NaN n(1)+n(3)];

% specify components of the state vector

indices.cl   =  1:N;                                % length: N
indices.ir   =  N+1:N+Nel;                          % length: Nel
indices.il   = N+Nel+1:N+2*Nel;                     % length: Nel
indices.phis   = N+2*Nel+1:N+3*Nel;                 % length: Nel
indices.phil   = N+3*Nel+1:2*N+3*Nel;               % length: N
indices.rcs   = 2*N+3*Nel+1:2*N+(3+m)*Nel;          % length: m*Nel;
indices.T = 2*N+(3+m)*Nel+1;                        % length: 1
precomp.indices = indices;
nvars = indices.T;
precomp.nvars=nvars;


%total mass matrix: 
precomp.M = spdiags([0; ones(n(1)-2,1); 0; ones(n(2)-2,1); 0; ones(n(3)-2,1); 0; zeros(3*Nel+N,1); 
                    repmat([ones(m-1,1); 0],Nel,1); 1],0,nvars,nvars);

% compute differentiation matrices and Clenshaw-Curtis weights
xi = cell(3,1); DM = cell(3,1); w = cell(3,1);

for k = 1:3
[xi{k},DM{k}] = chebdif(n(k),1); 
DM{k} = 2/p.d(k)*DM{k}; %DM{k}(:,:,2) = 4/p.d(k)^2*DM{k}(:,:,2);
w{k} = clencurt(n(k)-1)*p.d(k)/2;
end

% nodes along model axis: we map the domain discretizations onto model geometry 
precomp.x_neg = p.d(1)*(xi{1}+1)/2; precomp.x_sep = p.d(1)+ p.d(2)*(xi{2}+1)/2; precomp.x_pos = p.d(1)+p.d(2)+ p.d(3)*(xi{3}+1)/2;  

precomp.x = [p.d(1)*(xi{1}(1:end-1)+1)/2; p.d(1)+ p.d(2)*(xi{2}(1:end)+1)/2; p.d(1)+p.d(2) + p.d(3)*(xi{3}(2:end)+1)/2]; 
       
% nodes along model domains
     
       
% store differentiation matrices and Clenshaw-Curtis weights   
precomp.DM = DM;
precomp.w = w;
                         
% compute differentiation matrix for particle dimension and use it to compute matrices
% Ms, As, Bs for solid-state diffusion equation in the form Ms du/dt =
% Ds/Rs^2*As*u + Bs*ir (see paper for details)

M = 2*m; [rho,DM] = chebdif(M+1,2); wr = clencurt(M);
DN = DM(M+1,m+2:M+1,1)- DM(M+1,m:-1:1,1); DN2 = DM(m+2:M,m+2:M+1,2) - DM(m+2:M,m:-1:1,2); 

% nodes in the particle or pseudodimension
precomp.r = rho(m+2:M+1)*p.rp;            
precomp.rho = rho(m+2:M+1);
precomp.wr = wr(m+2:M+1)';

As = [DN2; DN - fliplr(eye(1,m))];
Bs = flipud(eye(m,1))/p.F;
%Cs = -(DM(M+1,m+2:M+1,1)+DM(m+1,m:-1:1,1))/DM(M+1,m+1,1);
%Ds = - 1/DM(M+1,m+1,1);

%replicate for each node in x-Dimension.  
precomp.As{1} = kron(eye(n(1)),As);
precomp.As{3} = kron(eye(n(3)),As);
BsReps = repmat({Bs},1,max(n(1),n(3)));
precomp.Bs{1} = blkdiag(BsReps{1:n(1)})/p.as(1);
precomp.Bs{3} = blkdiag(BsReps{1:n(3)})/p.as(3);




end