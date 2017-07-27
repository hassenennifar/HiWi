filenames = {'Eeq_neg_Graphite.txt','Eeq_pos_NMC.txt','dEeqdT_neg_Graphite.txt','dEeqdT_pos_NMC.txt'};
ylabels = {'$\phi_0$ [V] (Graphit)','$\phi_0$ [V]  (NMC)','$\partial \phi_0/\partial T$ [mV/K] (Graphit)','$\partial \phi_0/\partial T$ [mV/K] (NMC)'};
delimiter = ' ';
close all

for n =1:2
    filename = filenames{n};
    data = importdata(filename,delimiter);
    x = sort(data(:,1)); y = sort(data(:,2));
    pp = spline(x,y);
    xx = x(1):0.001:x(end);
    yy = ppval(pp,xx);
    subplot(1,2,n)
    plot(xx,yy,x,y,'k.')
    xlabel('St\"ochiometrie $x$','Interpreter','latex','FontSize',14)
    ylabel(ylabels{n},'Interpreter','latex','FontSize',14)
end

figure
m = [25 10];

for n =1:2
    filename = filenames{n+2};
    data = importdata(filename,delimiter);
    x = sort(data(:,1)); y = sort(data(:,2));
    pp = spline(x,y);
    pp2 = splinefit(x,y,m(n));
    xx = x(1):0.001:x(end);
    yy = ppval(pp,xx);
    yy2 = ppval(pp2,xx);
    subplot(1,2,n)
    plot(xx,yy,xx,yy2,x,y,'k.')
    leg = legend('Splines','Splinefit');
    set(leg,'Interpreter','latex','FontSize',13)
    xlabel('St\"ochiometrie $x$','Interpreter','latex','FontSize',14)
    ylabel(ylabels{n+2},'Interpreter','latex','FontSize',14)
end


