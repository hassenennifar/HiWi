function p = current(profile,p)

p.tf = 5000;

if strcmp(profile,'10C on/off')
    
    %input current: 10C dischrage/charge for 1s each
    p.x0 = [0.6 NaN 0.6];
    p.i = @(t)10*p.C*((1<t & t<=2)-(6<t & t<=7));
    p.tf = 20;
    
elseif strcmp(profile,'1C charge')
    
    %input current: constant 1C charge
    p.x0 = [0.6 NaN 0.6];
    p.i = @(t)-p.C;  
    
elseif strcmp(profile,'10C discharge')
    
    %input current: constant 10C discharge
    p.i = @(t)10*p.C; 
    
elseif strcmp(profile,'volatile')
    
    %input current: Switch between between 10C charge/discharge every 100ms
    p.x0 = [0.6 NaN 0.6];
    p.i = @(t)10*p.C*(cos(1*2*pi*t)+cos(10*2*pi*t)+cos(100*2*pi*t));
    p.tf = .2;
    
elseif strcmp(profile,'artemis')
    
    %input current: from artemis drive cycle
    filename = 'artemis.txt';
    delimiter= ' ';
    data = importdata(filename,delimiter);
    x= data(:,1)-data(1,1); y = data(:,2);
    pp = splinefit(x,y,600);
    p.i = @(t)p.C*ppval(pp,t);
    
else
    
    %input current: constant 1C discharge for at most 4000s
    %p.i = @(t)p.C;
    p.i = @(t)p.C*t;          
    p.tend = 5000;

end

end