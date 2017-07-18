function [g0, g1] = Eeq_neg(x, P, cur)

% g0 = 0.194+1.5.*exp(-120.0.*sto) ...
%        + 0.0351.*tanh((sto-0.286)./0.083) ...
%        - 0.0045.*tanh((sto-0.849)./0.119) ...
%        - 0.035.*tanh((sto-0.9233)./0.05) ...
%        - 0.0147.*tanh((sto-0.5)./0.034) ...
%        - 0.102.*tanh((sto-0.194)./0.142) ...
%        - 0.022.*tanh((sto-0.9)./0.0164) ...
%        - 0.011.*tanh((sto-0.124)./0.0226) ...
%        + 0.0155.*tanh((sto-0.105)./0.029);
%    
% g1 = -1.5.*(120.0./P.csmax_neg).*exp(-120.0.*sto) ...
%        +(0.0351./(0.083.*P.csmax_neg)).*((cosh((sto-0.286)./0.083)).^(-2)) ...
%        -(0.0045./(P.csmax_neg.*0.119)).*((cosh((sto-0.849)./0.119)).^(-2)) ...
%        -(0.035./(P.csmax_neg.*0.05)).*((cosh((sto-0.9233)./0.05)).^(-2)) ...
%        -(0.0147./(P.csmax_neg.*0.034)).*((cosh((sto-0.5)./0.034)).^(-2)) ...
%        -(0.102./(P.csmax_neg.*0.142)).*((cosh((sto-0.194)./0.142)).^(-2)) ...
%        -(0.022./(P.csmax_neg.*0.0164)).*((cosh((sto-0.9)./0.0164)).^(-2)) ...
%        -(0.011./(P.csmax_neg.*0.0226)).*((cosh((sto-0.124)./0.0226)).^(-2)) ...
%        +(0.0155./(P.csmax_neg.*0.029)).*((cosh((sto-0.105)./0.029)).^(-2));

% if cur < 0
%     g0 = ppval(P.Eeq_neg.ch.pp, sto);
%     g1 = ppval(P.Eeq_neg.ch.qq, sto)/P.csmax_neg;
% else
%     g0 = ppval(P.Eeq_neg.dch.pp, sto);
%     g1 = ppval(P.Eeq_neg.dch.qq, sto)/P.csmax_neg;
% end

g0 = 1.747+0.4682*x-3.449e-07*exp(13.69*x)+0.05748*tanh(-5.228*x+3.585)-...
    0.01585*tanh(30.28*x-16.35)+0.2144*tanh(9.132*x-9.344)-0.9348*...
    tanh(25.18*x-1.854)+0.2918*tanh(4.401*x-0.6359)-0.9978*...
    tanh(43.32*x-0.08442)+0.9335*tanh(24.08*x-1.783)-tanh(0.0706*x)-tanh(3.34*x);

% derivative from Wolfram
g1 = 0.4682-4.721681e-06*exp(13.69*x)-0.479938*sech(16.35-30.28*x).^2 ...
    -23.5383*sech(1.854-25.18*x).^2+1.9579*sech(9.344-9.132*x).^2 ...
    -0.300505*sech(3.585-5.228*x).^2-43.2247*sech(0.08442-43.32*x).^2 ...
    +22.4787*sech(1.783-24.08*x).^2+1.28421*sech(0.6359-4.401*x).^2 ...
    -0.0706*sech(0.0706*x).^2-3.34*sech(3.34*x).^2;
g1 = g1/P.csmax_neg;

% derivative from matlab syms
% g1 = (239969*tanh((757*x)/25 - 327/20)^2)/500000 - (8919003238880972347*exp((1369*x)/100))/1888946593147858085478400 + (1878159*tanh((1307*x)/250 - 717/200)^2)/6250000 + (2942283*tanh((1259*x)/50 - 927/500)^2)/125000 - (561967*tanh((602*x)/25 - 1783/1000)^2)/25000 - (152961*tanh((2283*x)/250 - 1168/125)^2)/78125 - (6421059*tanh((4401*x)/1000 - 6359/10000)^2)/5000000 + (5403087*tanh((1083*x)/25 - 4221/50000)^2)/125000 + (167*tanh((167*x)/50)^2)/50 + (353*tanh((353*x)/5000)^2)/5000 - 1119125271/25000000;
% g1 = g1/P.csmax_neg;