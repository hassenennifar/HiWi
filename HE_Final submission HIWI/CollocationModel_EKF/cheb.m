function [D,x,w] = cheb(N)
% CHEB  calculates Chebyshev collocation points on [-1,1], corresponding
% differentiation matrix and the Clenshaw-Curtis quadrature weights. The
% code is a minor modification of the two files cheb.m and clencurt.m  
% from L.N. Trefethen's book Spectral Methods in MATLAB.
%
% [D,x,w] = cheb(N)
%
% x - Chebyshev collocation points in ascending order
%      [column-vector of length N]
%
% D - differentiation matrix
%      [square matrix, size N]
%
% w - Clenshaw-Curtis quadrature weights 
%     [row-vector of length N]

N = N-1; % we use this to get outputs of desired size

if N==0, D=0; x=1; w = -1; return, end

theta = pi*(0:N)'/N; 
x = -cos(theta);            % -sign to change order
c = [2 ones(1,N-1) 2] .* ((-1).^(0:N));  % c = [2 -1 1 -1 ... 1 -1 2]
X = repmat(x,1,N+1);
D_numer = c' * (1./c);
D_denom = (X - X') + eye(N+1);
D = D_numer./D_denom;
D = D - diag(sum(D,2));

w = zeros(1,N+1); ii = 2:N; v = ones(N-1,1);
  if mod(N,2)==0 
    w(1) = 1/(N^2-1); w(N+1) = w(1);
    for k=1:N/2-1, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
    v = v - cos(N*theta(ii))/(N^2-1);
  else
    w(1) = 1/N^2; w(N+1) = w(1);
    for k=1:(N-1)/2, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
  end
  w(ii) = 2*v/N;
end