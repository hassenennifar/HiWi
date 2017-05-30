function [x, P, iter] = runModel(k,P,x,cur,ts,sol)

iter = 0;

while 1
    iter = iter+1;
    
    [g, J] = approxBattery(k,P,x,cur,ts,sol);
    
    deltax = J\g;
    

    xknew = x(:,k) + deltax;
    xknew = shoeHorns(xknew,x(:,k),1,P);

    x(:,k) = xknew;

    i = 30;
    if P.n_neg <= i
        i = 1;
    end

    errlim = P.tolabs*ones(P.nx,1);
    errlim(P.idx_jn) = errlim(P.idx_jn)*1e-6; % reduce absolute tolerance for jn
    if sum((abs(deltax(P.nx+1:P.nx*2)) < max(errlim, P.tolrel*abs(x(P.nx+1:P.nx*2,k)))) ...
        | (abs(deltax(P.nx*i+1:P.nx*(i+1)))< max(errlim, P.tolrel*abs(x(P.nx*i+1:P.nx*(i+1),k))))) == P.nx
        break;
    elseif max(abs(deltax)) < P.tolabs
        break;
    elseif iter > P.lim
        error('Solver could not converge: maximum iteration number reached');
    end

end