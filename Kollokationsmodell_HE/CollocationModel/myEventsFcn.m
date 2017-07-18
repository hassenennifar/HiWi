function [value,isterminal,direction] = myEventsFcn(t,~,p)
value = double(abs(p.i(t+1e-3)-p.i(t))<1e-7);
isterminal = 1;
direction = 0;
end
