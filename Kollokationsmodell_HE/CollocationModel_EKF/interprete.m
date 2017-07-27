function result = interprete(p,t,tjump,y)

result.p = p;

n = p.precomp.n;
m = p.precomp.m;
ind = p.precomp.indices;

result.t = t;
result.tjump = tjump;
result.y = y;
result.i = @(t)1e6*p.i(t);
result.T = y(:,ind.T);
result.cs = 1e6*y(:,ind.cs)*kron(eye(n(1)+n(3)),p.precomp.V')*diag([repmat(1./p.precomp.r(:,1),n(1),1); repmat(1./p.precomp.r(:,3),n(3),1)]);
result.cssurf = result.cs(:,m:m:end);
result.phis = y(:,ind.phis);
result.phil = y(:,ind.phil);
result.cl = 1e6*y(:,ind.cl);
result.ir = 1e6*y(:,ind.ir);
result.il = 1e6*y(:,ind.il);

%checking electrolyte conservation
w = [p.epsl(1)*p.precomp.w{1} zeros(1,n(2)-2) p.epsl(3)*p.precomp.w{3}]; w(n(1):n(1)+n(2)-1)=w(n(1):n(1)+n(2)-1)+p.epsl(2)*p.precomp.w{2};
result.intcl = result.cl*w'/ (sum(p.d.*p.epsl)); result.intcldelta = max(result.intcl) - min(result.intcl);
end