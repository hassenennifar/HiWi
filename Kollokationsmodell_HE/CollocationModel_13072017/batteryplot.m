function []=batteryplot(result,plotting)

close all;
precomp = result.p.precomp; n = precomp.n; m = precomp.m; r = 1e3*precomp.r;
x = 1e3*precomp.x; 
x_neg = x(1:n(1)); x_sep = x(n(1):n(1)+n(2)-1); x_pos =  x(n(1)+n(2)-1:sum(n)-2);
l_neg = x(n(1));
l_pos = x(n(1)+n(2)-1);
neg = 1:n(1); pos = n(1)+1:n(1)+n(3); sep = n(1):n(1)+n(2)-1;
t = result.t; tjump = result.tjump;

figure('Name','Time Plots');
subplot(3,1,1)
plot(t,result.i(t),tjump,result.i(tjump),'+')
xlabel('t [s]')
ylabel('i [A/m^2]')

subplot(3,1,2)
plot(t,result.phis(:,end))
xlabel('t [s]')
ylabel('U [V]')

subplot(3,1,3)
plot(t,result.T-273.15)
xlabel('t [s]')
ylabel('T [mol/m^3]')

%number of times
N = 10;
time = zeros(N,1);
labels = cell(N,1);
for k = 1:N
    time(k)=find(t >= (k-1)*t(end)/(N-1),1);
    %time(i)=numel(result.t-n+i);
    labels{k}= ['t = ' num2str(t(time(k))) ' s'];
end

if strcmp(plotting,'all')

figure('Name','Internal States');
ColOrd = get(gca,'ColorOrder');

%surface solid state concentration
h1 = subplot(4,2,1);
axes(h1)
hold on
for k = 0:N-1
   plot(x_neg,result.cssurf(time(k+1),neg),'Color',ColOrd(1+mod(k,7),:))
   plot(x_pos,result.cssurf(time(k+1),pos),'Color',ColOrd(1+mod(k,7),:))
end
line([l_neg l_neg],get(h1,'YLim'),'Color',[0.1 0.1 0.1],'LineStyle',':')
line([l_pos l_pos],get(h1,'YLim'),'Color',[0.1 0.1 0.1],'LineStyle',':')
xlabel('x [\mum]')
ylabel('c_s [mol/m^3]')

% electrolyte concentration
h2 = subplot(4,2,2);
axes(h2)
hold on
for k = 0:N-1
   plot(x,result.cl(time(k+1),:),'Color',ColOrd(1+mod(k,7),:))
end
line([l_neg l_neg],get(h2,'YLim'),'Color',[0.1 0.1 0.1],'LineStyle',':')
line([l_pos l_pos],get(h2,'YLim'),'Color',[0.1 0.1 0.1],'LineStyle',':')
xlabel('x [\mum]')
ylabel('c_l [mol/m^3]')

% reaction current
h3 = subplot(4,2,3);
axes(h3)
hold on
for k = 0:N-1
   plot(x_neg,result.ir(time(k+1),neg),'Color',ColOrd(1+mod(k,7),:))
   plot(x_pos,result.ir(time(k+1),pos),'Color',ColOrd(1+mod(k,7),:))
end
line([l_neg l_neg],get(h3,'YLim'),'Color',[0.1 0.1 0.1],'LineStyle',':')
line([l_pos l_pos],get(h3,'YLim'),'Color',[0.1 0.1 0.1],'LineStyle',':')
xlabel('x [\mum]')
ylabel('i_r [A/m^3]')

% electrolyte current
h4 = subplot(4,2,4);
axes(h4)
hold on
for k = 0:N-1
   plot(x_neg,result.il(time(k+1),neg),'Color',ColOrd(1+mod(k,7),:))
   plot(x_sep,repmat(result.i(time(k+1)),n(2)),'Color',ColOrd(1+mod(k,7),:))
   plot(x_pos,result.il(time(k+1),pos),'Color',ColOrd(1+mod(k,7),:))
end
line([l_neg l_neg],get(h4,'YLim'),'Color',[0.1 0.1 0.1],'LineStyle',':')
line([l_pos l_pos],get(h4,'YLim'),'Color',[0.1 0.1 0.1],'LineStyle',':')
xlabel('x [\mum]')
ylabel('i_l [A/m^2]')

% solid-phase potential
h5 = subplot(4,2,5);
axes(h5)
hold on
for k = 0:N-1
   plot(x_neg,result.phis(time(k+1),neg),'Color',ColOrd(1+mod(k,7),:))
   plot(x_pos,result.phis(time(k+1),pos),'Color',ColOrd(1+mod(k,7),:))
end
line([l_neg l_neg],get(h5,'YLim'),'Color',[0.1 0.1 0.1],'LineStyle',':')
line([l_pos l_pos],get(h5,'YLim'),'Color',[0.1 0.1 0.1],'LineStyle',':')
xlabel('x [\mum]')
ylabel('\phi_s [V]')

% electrolyte potential
h6 = subplot(4,2,6);
axes(h6)
hold on
for k = 0:N-1
   plot(x,result.phil(time(k+1),:)) 
end
line([l_neg l_neg],get(h6,'YLim'),'Color',[0.1 0.1 0.1],'LineStyle',':')
line([l_pos l_pos],get(h6,'YLim'),'Color',[0.1 0.1 0.1],'LineStyle',':')
xlabel('x [\mum]')
ylabel('\phi_l [V]')

% anode solid state potential in r-dimension
subplot(4,2,7);
hold on
for k = 0:N-1
   plot(r(:,1),result.cs(time(k+1),1:m)) 
end
xlabel('r [\mum]')
ylabel('c_s [mol/m^3]')

% cathode solid state potential in r-dimension
subplot(4,2,8);
hold on
for k = 0:N-1
   plot(r(:,3),result.cs(time(k+1),end-m+1:end)) 
end
xlabel('r [\mum]')
ylabel('c_s [mol/m^3]')

legend(labels,'Location','Best')

end

figure('Name','Erhaltung der Elektrolytkonzentration');
plot(t,result.intcl)
a = annotation('textbox',[.2 .5 .3 .3],'String',['max. Delta = ' ...
num2str(result.intcldelta)], 'FitBoxToText','on'); a.Interpreter = 'latex'; a.FontSize = 12;

end