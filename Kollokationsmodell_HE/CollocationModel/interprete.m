function result = interprete(p,t,y)

result.p = p;

n = p.precomp.n;
m = p.precomp.m;
ind = p.precomp.indices;

result.t = t;
result.y = y;
result.i = 1e6*p.i(t);
result.T = y(:,ind.T);
result.cs = 1e6*y(:,ind.rcs)*diag([repmat(1./p.precomp.r(:,1),n(1),1);repmat(1./p.precomp.r(:,3),n(3),1)]);
result.cssurf = result.cs(:,m:m:end);
result.phis = y(:,ind.phis);
result.phil = y(:,ind.phil);
result.cl = 1e6*y(:,ind.cl);
result.ir = 1e6*y(:,ind.ir);
result.il = 1e6*y(:,ind.il);

end