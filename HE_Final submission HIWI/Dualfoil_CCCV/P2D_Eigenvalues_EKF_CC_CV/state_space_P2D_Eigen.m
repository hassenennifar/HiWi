function [Fk, Hk] = state_space_P2D_Eigen(x_est, k, P)

%% model matrix Fk
cl = x_est(P.idx_cl:P.nx:P.nx*P.nj, :);

Dl = fun_Dl(cl(:,k-1), P);
tp = fun_tp(cl(:,k-1), P);

M = zeros(P.nj);
Fl = zeros(P.nj);
bl = zeros(P.nj);

j = 1;
epsl = P.epsl_neg;
hx = P.L_neg/P.n_neg;
Fl(j,j:j+1) = [epsl*hx/P.dt*3/8 epsl*hx/P.dt/8];
M(j,j:j+1) = [epsl*hx/P.dt*3/8+epsl^1.5/hx*(Dl(j)+Dl(j+1))/2, epsl*hx/P.dt/8-epsl^1.5/hx*(Dl(j)+Dl(j+1))/2];
bl(j,j) = (1-(tp(j)+tp(j+1))/2)/P.F;
for j = 2:P.nj-1
    
    if j < P.bnd_sep_neg
        epsl1 = P.epsl_neg;
        hx1 = P.L_neg/P.n_neg;
        epsl2 = epsl1;
        hx2 = hx1;
    elseif j == P.bnd_sep_neg
        epsl1 = P.epsl_neg;
        hx1 = P.L_neg/P.n_neg;
        epsl2 = P.epsl_sep;
        hx2 = P.L_sep/P.n_sep;
    elseif j < P.bnd_pos_sep
        epsl1 = P.epsl_sep;
        hx1 = P.L_sep/P.n_sep;
        epsl2 = epsl1;
        hx2 = hx1;
    elseif j == P.bnd_pos_sep
        epsl1 = P.epsl_sep;
        hx1 = P.L_sep/P.n_sep;
        epsl2 = P.epsl_pos;
        hx2 = P.L_pos/P.n_pos;
    else
        epsl1 = P.epsl_pos;
        hx1 = P.L_pos/P.n_pos;
        epsl2 = epsl1;
        hx2 = hx1;
    end
    
    Fl(j, j-1:j+1) = [epsl1*hx1/P.dt/8 (epsl1*hx1/P.dt+epsl2*hx2/P.dt)*3/8 epsl2*hx2/P.dt*1/8];
    M(j, j-1:j+1) = [epsl1*hx1/P.dt/8-epsl1^1.5/hx1*(Dl(j)+Dl(j-1))/2, ...
        (epsl1*hx1/P.dt+epsl2*hx2/P.dt)*3/8+epsl1^1.5/hx1*(Dl(j)+Dl(j-1))/2+epsl2^1.5/hx2*(Dl(j)+Dl(j+1))/2, ...
        epsl2*hx2/P.dt/8-epsl2^1.5/hx2*(Dl(j)+Dl(j+1))/2];
    bl(j,j-1:j) = [-(1-(tp(j)+tp(j-1))/2)/P.F, (1-(tp(j)+tp(j+1))/2)/P.F];
end
j = P.nj;
epsl = P.epsl_pos;
hx = P.L_pos/P.n_pos;
Fl(j, j-1:j) = [epsl2*hx2/P.dt/8, epsl2*hx2/P.dt*3/8];
M(j, j-1:j) = [epsl2*hx2/P.dt/8-epsl^1.5/hx*(Dl(j)+Dl(j-1))/2, ...
        epsl2*hx2/P.dt*3/8+epsl^1.5/hx*(Dl(j)+Dl(j-1))/2];
bl(j,j-1) = -(1-(tp(j)+tp(j-1))/2)/P.F;

Fl = inv(M)*Fl;
bl = inv(M)*bl;


Fk = zeros(P.nj*P.ne);
bk = zeros(P.nj*P.ne, P.nj*2);



j = 1;
Ds = P.Ds_neg;
Rp = P.Rp_neg;
Fs = 1; % 0 1/(1+P.dt*30*Ds/Rp/Rp)];
bs = -3/Rp;% -45/2/Rp/Rp];

Fk((j-1)*P.ne+1,1:P.ne:end) = Fl(j,:);
Fk((j-1)*P.ne+2,(j-1)*P.ne+2) = Fs;
bk((j-1)*P.ne+1,1:P.ne:end) = bl(j,:);
bk((j-1)*P.ne+2, (j-1)*P.ne+2) = bs;

for j = 2:P.nj-1
    
    if j > P.bnd_sep_neg && j < P.bnd_pos_sep
        Fs = 1;% 0 1];
    else
        if j<= P.bnd_sep_neg
            Ds = P.Ds_neg;
            Rp = P.Rp_neg;
        elseif j>= P.bnd_pos_sep
            Ds = P.Ds_pos;
            Rp = P.Rp_pos;
        end
       
        Fs = 1;% 0 1/(1+P.dt*30*Ds/Rp/Rp)];
        bs = -3/Rp;% -45/2/Rp/Rp];
    end
    Fk((j-1)*P.ne+1,1:P.ne:end) = Fl(j,:);
    Fk((j-1)*P.ne+2,(j-1)*P.ne+2) = Fs;
    bk((j-1)*P.ne+1,1:2:end) = bl(j,:);
    bk((j-1)*P.ne+2, (j-1)*2+2) = bs;
end

j = P.nj;

Ds = P.Ds_pos;
Rp = P.Rp_pos;
Fs = 1;% 0 1/(1+P.dt*30*Ds/Rp/Rp)];
bs = -3/Rp;% -45/2/Rp/Rp];

Fk((j-1)*P.ne+1,1:P.ne:end) = Fl(j,:);
Fk((j-1)*P.ne+2,(j-1)*P.ne+2) = Fs;
bk((j-1)*P.ne+1,1:P.ne:end) = bl(j,:);
bk((j-1)*P.ne+2, (j-1)*2+2) = bs;

%% output matrix Hk
cl = x_est(P.idx_cl:P.nx:P.nx*P.nj, :);
phil = x_est(P.idx_phil:P.nx:P.nx*P.nj, :);
jn = x_est(P.idx_jn:P.nx:P.nx*P.nj, :);
css = x_est(P.idx_css:P.nx:P.nx*P.nj, :);

Hk = zeros(1, P.nj*P.ne);
j = 1;
Rp = P.Rp_neg;
Ds = P.Ds_neg;
rka = P.rka_neg;
dcss_dc = 1;
dcss_dq = 8*Rp/35;% - 48/7*P.dt*Ds/Rp;%8*Rp/35/Ds;
[~, dEeq] = Eeq_neg(css(j,k)/P.csmax_neg,P);
i0 = P.F*rka*sqrt(P.csmax_neg-css(j,k))*sqrt(css(j,k))*sqrt(cl(j,k));
dphis_dcl = -P.R*P.T/2/cl(j,k)*jn(j,k)/sqrt(i0^2-(P.F*jn(j,k)/2)^2);
dphis_dcss  = dEeq - P.R*P.T/2*(1/css(j,k) - 1/(P.csmax_neg-css(j,k)))*jn(j,k)/sqrt(i0^2-(P.F*jn(j,k)/2)^2);
dphis_dc  = dcss_dc*dphis_dcss;
% dphis_dq  = dcss_dq*dphis_dcss;
Hk(1:2) = -[dphis_dcl dphis_dc];
Bs = zeros(size(jn(:,k)'));
Bs(1) = -(P.Rf_neg + P.R*P.T/i0/sqrt(1-(P.F*jn(j,k)/2/i0)^2)) - Rp/35/Ds*dphis_dc;

j = P.nj;
Rp = P.Rp_pos;
Ds = P.Ds_pos;
rka = P.rka_pos;
dcss_dc = 1;
dcss_dq = 8*Rp/35;%- 48/7*P.dt*Ds/Rp;%8*Rp/35/Ds;
[~, dEeq] = Eeq_pos(css(j,k)/P.csmax_pos,P);
i0 = P.F*rka*sqrt(P.csmax_pos-css(j,k))*sqrt(css(j,k))*sqrt(cl(j,k));
dphis_dcl = -P.R*P.T/2/cl(j,k)*jn(j,k)/sqrt(i0^2-(P.F*jn(j,k)/2)^2);
dphis_dcss  = dEeq - P.R*P.T/2*(1/css(j,k) - 1/(P.csmax_pos-css(j,k)))*jn(j,k)/sqrt(i0^2-(P.F*jn(j,k)/2)^2);
dphis_dc  = dcss_dc*dphis_dcss;
% dphis_dq  = dcss_dq*dphis_dcss;
Hk(P.nj*P.ne-1:P.nj*P.ne) = [dphis_dcl dphis_dc];
Bs(end) = P.Rf_pos + P.R*P.T/i0/sqrt(1-(P.F*jn(j,k)/2/i0)^2) - Rp/35/Ds*dphis_dc;