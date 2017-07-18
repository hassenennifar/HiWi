filename = 'Eeq_pos_NMC_lit.txt';
delimiter = ' ';
data = importdata(filename,delimiter);
x = data(:,1); y = data(:,2);
breaks = [0 0.9  0.975 1];
p.Eeq{3}  = splinefit(x,y,breaks);
p.dEeqdx{3}  = ppdiff(p.Eeq{3},1);