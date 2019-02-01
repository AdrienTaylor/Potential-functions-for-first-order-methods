clear all;

verbose = 1;
tolerance = 1e-8;

folder = 'SaveData/';
nname  = 'Full_data1.dat';
pplot  = 0;
ssave  = 0;

% Problem
L = 1;
m = 0;
N = 50;
n = 3;

% R = 1;
% sig = 1;

% P = [ x0 z0 y1 | x11 ... x1n | f'(x0) | G(y1,1) ... G(y1,n) | f'(x11) ... f'(x1n) | w0]
% F = [                           f(x0) |       f(y1)         | f(x1)]

dimG = 4 + 3 * n + 1;
dimF = 3;

nbPts = 3 + n ; % xs x0 y1 x11 ... x1n


x0 = zeros(1, dimG); x0(1) = 1;
z0 = zeros(1, dimG); z0(2) = 1;
w0 = zeros(1, dimG); w0(end) = 1;
y1 = zeros(1, dimG); y1(3) = 1;
x1 = zeros(n, dimG); x1(:,4:3+n) = eye(n);
xs = zeros(1, dimG);

gx0 = zeros(1, dimG); gx0(4+n) = 1;
di  = zeros(n, dimG); di(:,5+n:4+2*n) = eye(n);
gy1 = sum(di,1)/n;
gx1 = zeros(n, dimG); gx1(:,5+2*n:4+3*n) = eye(n);
gxs = zeros(1, dimG);

fx0 = zeros(1, dimF); fx0(1) = 1;
fy1 = zeros(1, dimF); fy1(2) = 1;
fx1 = zeros(1, dimF); fx1(3) = 1;
fxs = zeros(1, dimF);

% Structure of first-step Lyapunov
V0 = [0 0; 0 0];
v0 = 0;
B0 = 1;

% Structure of last-step Lyapunov
tau = sdpvar(1);
VN = [0 0; 0 0];
vN = tau;
BN = 0;

%%

M        = 1/2/(L-m) *[  -L*m, L*m,  m,   -L; L*m,  -L*m,  -m,    L;    m,    -m,  -1,  1 ;    -L,   L,   1,  -1];

V{1} = V0;
v{1} = v0;
B{1} = B0;
V{N+1} = VN;
v{N+1} = vN;
B{N+1} = BN;
acc    = 0;
cons = (tau==N);
for i = 1 : N
    if i < N
        V{i+1} = zeros(2);
        v{i+1} = sdpvar(1);
        B{i+1} = sdpvar(1);
    end
    XX = { xs,  x0,  y1};
    GG = {gxs, gx0, gy1};
    FF = {fxs, fx0, fy1};
    for k = 1:n
        XX{3+k}   = x1(k,:);
        GG{3+k}   = gx1(k,:);
        FF{3+k}   = fx1;
    end
    lam{i}   = sdpvar(nbPts,nbPts,n,'full'); % for pep
    
    Z_def{i} = sdpvar(4);
    W_def{i} = sdpvar(4);
    lam_LSY{i} = sdpvar(3,1);
    lam_SIG{i} = sdpvar(1);
    lam_LSX{i} = sdpvar(5,1);
    
    SP0 = [x0.' gx0.']; SF0 = fx0.';
    
    dist_end = zeros(dimG);
    f_end    = zeros(1, dimF);
    z_dist_e = zeros(dimG);
    w_dist_e = zeros(dimG);
    VARY     = zeros(dimG);
    for k = 1:n
        SP1 = [x1(k,:).' gx1(k,:).'];
        SF1 = fx1.';
        dist_end = dist_end + SP1 * V{i+1} * SP1.';
        f_end    = f_end + v{i+1} * SF1.';
        z_dist_e = z_dist_e + [z0.' (w0-z0).' (y1-z0).' di(k,:).'] * Z_def{i} * [z0.' (w0-z0).' (y1-z0).' di(k,:).'].';
        w_dist_e = w_dist_e + [z0.' (w0-z0).' (y1-z0).' di(k,:).'] * W_def{i} * [z0.' (w0-z0).' (y1-z0).' di(k,:).'].';
        VARY = VARY + (gy1-di(k,:)).'*(gy1-di(k,:));
    end
    dist_end = dist_end/n;
    f_end    = f_end/n;
    z_dist_e = z_dist_e/n;
    VARY     = VARY/n;
    
    ExactLSY1 = gy1.' * (x0 - y1); ExactLSY1 = (ExactLSY1+ExactLSY1.')/2;
    ExactLSY2 = gy1.' * (x0 - z0); ExactLSY2 = (ExactLSY2+ExactLSY2.')/2;
    ExactLSY3 = gy1.' * (x0 - w0); ExactLSY3 = (ExactLSY3+ExactLSY3.')/2;
    
    cons_SDP{i} = - SP0 * V{i} * SP0.' - L/2 * (z0 - xs).' * (z0 - xs) - L/2 * B{i} * (w0 - xs).' * (w0 - xs) - lam_SIG{i} * VARY + dist_end + L/2 * z_dist_e + L/2 * w_dist_e +...
        lam_LSY{i}(1) * ExactLSY1 + lam_LSY{i}(2) * ExactLSY2 + lam_LSY{i}(3) * ExactLSY3;
    
    
    for k = 1:n
        ExactLSX1 = gx1(k,:).' * di(k,:); ExactLSX1 = (ExactLSX1+ExactLSX1.')/2;
        ExactLSX2 = gx1(k,:).' * (y1 - x1(k,:)); ExactLSX2 = (ExactLSX2+ExactLSX2.')/2;
        ExactLSX3 = gx1(k,:).' * (y1 - z0); ExactLSX3 = (ExactLSX3+ExactLSX3.')/2;
        ExactLSX4 = gx1(k,:).' * (y1 - x0); ExactLSX4 = (ExactLSX4+ExactLSX4.')/2;
        ExactLSX5 = gx1(k,:).' * (y1 - w0); ExactLSX5 = (ExactLSX5+ExactLSX5.')/2;
        cons_SDP{i} =  cons_SDP{i}  + lam_LSX{i}(1) * ExactLSX1 + lam_LSX{i}(2) * ExactLSX2 + lam_LSX{i}(3) * ExactLSX3 + lam_LSX{i}(4) * ExactLSX4 + lam_LSX{i}(5) * ExactLSX5;
    end
    
    cons_LIN{i} = - v{i}* SF0.'        + f_end;
    for j = 1:nbPts
        for k = 1:nbPts
            if j ~= k
                xi = XX{j}; xj = XX{k};
                gi = GG{j}; gj = GG{k};
                fi = FF{j}; fj = FF{k};
                TT = [xi; xj; gi; gj];
                
                cons_SDP{i} = cons_SDP{i} + lam{i}(j,k) * TT.' * M * TT;
                cons_LIN{i} = cons_LIN{i} + lam{i}(j,k) * (fi - fj);
            end
        end
    end
    cons = cons + (cons_SDP{i} <= 0);
    cons = cons + (Z_def{i} >= 0);
    cons = cons + (W_def{i} >= 0);
    cons = cons + (Z_def{i}(1,1) == 1);
    cons = cons + (W_def{i}(1,1) == B{i+1});
    cons = cons + (cons_LIN{i} == 0);
    cons = cons + (lam{i} >= 0);
    cons = cons + (lam_SIG{i} >= 0);
    acc  = acc  + lam_SIG{i};
end
obj =  acc;


solver_opt = sdpsettings('solver','mosek','verbose',verbose,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',tolerance);
solverDetails=optimize(cons,obj,solver_opt);
obj_out = double(obj);
%% Try to grasp what happens...
if pplot
    close all;
    dimGL = 2; dimFL = 1;
    VVal = zeros(dimGL, dimGL, N+1);
    vval = zeros(dimFL, 1  , N+1);
    tolp = .5e-4;
    for i = 1:N+1
        VVal(:,:,i) = double(V{i});
        vval(:,:,i) = double(v{i});
    end
    for i = 1:N
        ssig(i)     = double(lam_SIG{i});
        ssig_acc(i) = sum(ssig(1:i));
    end
    for i = 1:dimGL
        for j = 1:dimGL
            tt = sprintf('V%d%d',i,j);
            subplot(dimGL+2,dimGL,(i-1)*dimGL+j); title(tt); hold on;
            if max(abs(VVal(i,j,:))>=tolp)
                plot(squeeze(VVal(i,j,:)),'-b')
                axis([0 N+2 min(squeeze(VVal(i,j,:)))-1 max(squeeze(VVal(i,j,:)))+1]);
            else
                plot(squeeze(VVal(i,j,:)),'--r')
                axis([0 N+2 min(squeeze(VVal(i,j,:)))-1 max(squeeze(VVal(i,j,:)))+1]);
            end
        end
    end
    
    for i = 1:dimFL
        tt = sprintf('v%d', i);
        subplot(dimGL+2,dimGL,(dimGL)^2+i); title(tt); hold on;
        if max(abs(vval(i, 1,:))>=tolp)
            plot(squeeze(vval(i, 1,:)),'-b')
            axis([0 N+2 min(vval(i, 1,:))-1 max(vval(i, 1,:))+1]);
        else
            plot(squeeze(vval(i, 1,:)),'--r')
            axis([0 N+2 min(vval(i, 1,:))-1 max(vval(i, 1,:))+1]);
        end
    end
    if n>1
        tt = sprintf('ek');
        subplot(dimGL+2,dimGL,(dimGL)^2+dimFL+1); title(tt); hold on;
        if max(abs(ssig)>=tolp)
            plot(ssig,'-b');
        else
            plot(ssig,'--r');
        end
        axis([0 N+2 min(ssig)-1 max(ssig)+1]);
        tt = sprintf('ek (accumulated)');
        subplot(dimGL+2,dimGL,(dimGL)^2+dimFL+2); title(tt); hold on;
        if max(abs(ssig_acc)>=tolp)
            plot(ssig_acc,'-b');
        else
            plot(ssig_acc,'--r');
        end
        axis([0 N+2 min(ssig_acc)-1 max(ssig_acc)+1]);
    end
end

% method:
%   y1  = x0 + alpha * (z0 - x0)
%   x1i = y1 + beta  * (y1 - x0) + gamma * (y1 - z0) + delta * G(y1; i)
%   z1i = z1 + eps   * G(y1; i)
if pplot
    figure;
    for i = 1:N
        alpha(i) = double(lam_LSY{i}(2))/double(lam_LSY{i}(1));
        beta(i)  = double(lam_LSX{i}(4))/double(lam_LSX{i}(2));
        gamma(i) = double(lam_LSX{i}(3))/double(lam_LSX{i}(2));
        delta(i) = double(lam_LSX{i}(1))/double(lam_LSX{i}(2));
        eps(i)   = double(Z_def{i}(1,2));
        % check:
        if min(eig(double(Z_def{i})))>1e-5
            fprintf('Care: a matrix has not enough low eigenvalues (index %d, min eig %5.4f)\n',i,min(eig(double(Z_def{i}))))
        end
    end
    subplot(2,3,1); title('alpha_k');hold on;
    if max(abs(alpha)>=tolp)
        plot(alpha,'-b');
    else
        plot(alpha,'--r');
    end
    axis([0 N+2 min(alpha)-1 max(alpha)+1]);
    subplot(2,3,2); title('beta_k');hold on;
    if max(abs(beta)>=tolp)
        plot(beta,'-b');
    else
        plot(beta,'--r');
    end
    axis([0 N+2 min(beta)-1 max(beta)+1]);
    subplot(2,3,3); title('gamma_k');hold on;
    if max(abs(gamma)>=tolp)
        plot(gamma,'-b');
    else
        plot(gamma,'--r');
    end
    axis([0 N+2 min(gamma)-1 max(gamma)+1]);
    subplot(2,3,4); title('delta_k');hold on;
    if max(abs(delta)>=tolp)
        plot(delta,'-b');
    else
        plot(delta,'--r');
    end
    axis([0 N+2 min(delta)-1 max(delta)+1]);
    subplot(2,3,5); title('eps_k');hold on;
    if max(abs(eps)>=tolp)
        plot(eps,'-b');
    else
        plot(eps,'--r');
    end
    axis([0 N+2 min(eps)-1 max(eps)+1]);
    if ssave
        labels{1} = 'k'; labels{2} = 'V11'; labels{3} = 'V12'; labels{4} = 'V22'; labels{5} = 'v1'; labels{6} = 'ek'; labels{7} = 'ekc'; labels{8} = 'alphak'; labels{9} = 'betak'; labels{10} = 'gammak'; labels{11} = 'deltak'; labels{12} = 'epsk';
        data = [(1:size(VVal,3)).' squeeze(VVal(1,1,:))  squeeze(VVal(2,1,:)) squeeze(VVal(2,2,:)) squeeze(vval(1,1,:)) [ 0 ssig].' [ 0 ssig_acc].' [ 0 alpha].' [ 0 beta].' [ 0 gamma].' [ 0 delta].' [ 0 eps].'];
        saveData([folder nname],data,labels);
        fprintf('DATA SAVED\n');
    end
end

%
%
%
%     lam_LSY{i} = sdpvar(2,1);
%     lam_LSX{i} = sdpvar(4,1);
%
%     ExactLSY1 = gy1.' * (x0 - y1); ExactLSY1 = (ExactLSY1+ExactLSY1.')/2;
%     ExactLSY2 = gy1.' * (x0 - z0); ExactLSY2 = (ExactLSY2+ExactLSY2.')/2;
%         ExactLSX1 = gx1(k,:).' * di(k,:); ExactLSX1 = (ExactLSX1+ExactLSX1.')/2;
%         ExactLSX2 = gx1(k,:).' * (y1 - x1(k,:)); ExactLSX2 = (ExactLSX2+ExactLSX2.')/2;
%         ExactLSX3 = gx1(k,:).' * (y1 - z0); ExactLSX3 = (ExactLSX3+ExactLSX3.')/2;
%         ExactLSX4 = gx1(k,:).' * (y1 - x0); ExactLSX4 = (ExactLSX4+ExactLSX4.')/2;





