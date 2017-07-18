clear; 
N = 30; 
 
%convergence plots
n0= [10 10 10]; m0 = 10;
convergence = zeros(N,7);
allResults = cell(N,1);

for k = 1:N
    k
    results = simulation(n0,m0+k);
    allResults{k}=results;
    convergence(k,1)=results.phis(end,end);
    convergence(k,2)=results.cssurf(end,1);
    convergence(k,3)=results.cssurf(end,end);
    convergence(k,4)=results.phil(end,1);
    convergence(k,5)=results.cl(end,1);
    convergence(k,6)=results.ir(end,end);       
end


close all;
figure('Name','Convergence-Plot')

subplot(3,2,1)
plot(n0(1)+(1:N),convergence(:,1))
legend('V')
subplot(3,2,2)
plot(n0(1)+(1:N),convergence(:,2))
legend('c_s left boundary')
subplot(3,2,3)
hold on
plot(n0(1)+(1:N),convergence(:,3))
legend('c_s right boundary')
subplot(3,2,4)
plot(n0(1)+(1:N),convergence(:,4))
legend('phil left boundary')
subplot(3,2,5)
plot(n0(1)+(1:N),convergence(:,5))
legend('cl left boundary')
subplot(3,2,6)
plot(n0(1)+(1:N),convergence(:,6))
legend('ir right boundary')








