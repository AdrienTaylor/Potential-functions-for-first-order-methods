clear all;
clc;

verbose = 1; tolerance = 1e-8;
ssave = 0;
pplot = 1;
folder = 'SaveData/';
nname  = 'overparam_design_pavging.dat';
relax  = 0;
% Problem
L = 1;
m = 0;
N = 10;
n = 2;

% P = [ xk zk yk1 | xk1 ... xkn | f1'(xk) ... fn'(xk) | f1'(yk1) ... fn'(xk1) | f1'(xk1) ... fn'(xk1) | ... | f1'(xkn) ... fn'(xkn) | f1'(xs) ... fn'(xs)]
% F = [                    fi(xk)  |  fi(yk1)  | fi(xki) fj(xki) ]

dimG = 3 + 3*n + n^2 + (n-1); dimF = 4;
nbPts = n+3; % xs, xk, yk1, xk1, ..., xkn

xk = zeros(1, dimG); xk(1) = 1;
zk = zeros(1, dimG); zk(2) = 1;
yk1= zeros(1, dimG); yk1(3)= 1;
xs = zeros(1, dimG);
xk1= zeros(n, dimG); xk1(:,4:3+n) = eye(n);

gxk = zeros(n, dimG);gxk(:,4+n:3+2*n) = eye(n);
GXK = sum(gxk,1)/n;
gyk1 = zeros(n, dimG);gyk1(:,4+2*n:3+3*n) = eye(n);
GYK1 = sum(gyk1,1)/n;
gxk1= zeros(n, dimG, n);
GXK1= zeros(n, dimG);
gxs = zeros(n, dimG); gxs(1:n-1,dimG-(n-2):dimG) = eye(n-1);
gxs(end,:) = - sum(gxs,1);
GXS = sum(gxs,1)/n;

fxk = zeros(n, dimF); fxk(:,1) = 1;
FXK = sum(fxk,1)/n;
fyk1 = zeros(n, dimF); fyk1(:,2) = 1;
FYK1 = sum(fyk1,1)/n;
fxk1= zeros(n, dimF, n);
FXK1= zeros(n, dimF);
fxs = zeros(n, dimF);
FXS = sum(fxs,1)/n;

s_index = 4+3*n; e_index = s_index + n-1;
for i = 1:n
    gxk1(:,s_index:e_index,i) = eye(n);
    s_index = s_index + n; e_index = s_index + n-1;
    fxk1(:,4,i) = 1;
    fxk1(i,:,i) = [0 0 1 0];
    FXK1(i,:)   = sum(fxk1(:,:,i),1)/n;
    GXK1(i,:)   = sum(gxk1(:,:,i),1)/n;
end

% potential part in terms of x's
% term 1: || x - xs || ^2
% term 2: || f'(y)  || ^2
% term 3: || f'i(y) -f'i(xs)|| ^2
% term 4: < f'(y); y-xs >

term1_k  = (xk-xs).'*(xk-xs);
term1_k1 = zeros(dimG);
for i = 1:n
    term1_k1 = term1_k1 + (xk1(i,:)-xs).'*(xk1(i,:)-xs)/n;
end

term2_k  = GXK.'*GXK;
term2_k1 = zeros(dimG);
for i = 1:n
    term2_k1 = term2_k1 + GXK1(i,:).'*GXK1(i,:)/n;
end

term3_k  = zeros(dimG);
for i = 1:n
    term3_k = term3_k + (gxk(i,:)-gxs(i,:)).'*(gxk(i,:)-gxs(i,:))/n;
end

term3_k1 = zeros(dimG);
for i = 1:n
    for j = 1:n
        term3_k1 = term3_k1 + (gxk1(i,:,j)-gxs(i,:)).'*(gxk1(i,:,j)-gxs(i,:))/n^2;
    end
end

term4_k  = (xk-xs).'*GXK; term4_k = 1/2 * (term4_k+term4_k.');
term4_k1 = zeros(dimG);
for i = 1:n
    temp     = (xk1(i,:)-xs).'*GXK1(i,:); temp = 1/2 * (temp+temp.');
    term4_k1 = term4_k1 + temp/n;
end

statesKf = (FXK - FXS);
statesK1f = (sum(FXK1,1)/n-FXS);
% line-searches

consLSY{1} = GYK1.' * (xk - yk1); consLSY{1} = (consLSY{1}+consLSY{1}.')/2;
consLSY{2} = GYK1.' * (zk -  xk); consLSY{2} = (consLSY{2}+consLSY{2}.')/2;
for i = 1:n
    consLS1{1,i} = GXK1(i,:).' * (yk1 - xk1(i,:));  consLS1{1,i} = (consLS1{1,i}+consLS1{1,i}.')/2;
    consLS1{2,i} = GXK1(i,:).' * (gyk1(i,:));   consLS1{2,i} = (consLS1{2,i}+consLS1{2,i}.')/2;
end

var = zeros(dimG);
for i = 1:n
    var = var + (gxs(i,:)-GXS).'*(gxs(i,:)-GXS)/n;
end

% potentials' coefficients
Q1_0 = 0;
Q2_0 = 0;
Q3_0 = 0;
Q4_0 = 0; Q0 = [Q1_0; Q2_0; Q3_0; Q4_0];
ap_0 = L/2;
d_0  = 0;

tau = sdpvar(1);

Q1_N = 0;
Q2_N = 0;
Q3_N = 0;
Q4_N = 0; QN = [Q1_N; Q2_N; Q3_N; Q4_N];
ap_N  = L/2;
d_N   = 9/10*N;

%%

M        = 1/2/(L-m) *[  -L*m, L*m,  m,   -L; L*m,  -L*m,  -m,    L;    m,    -m,  -1,  1 ;    -L,   L,   1,  -1];

Vx{1}   = Q0;
vx{1}   = d_0;
Vx{N+1} = QN;
vx{N+1} = d_N;
Vz{1}   = ap_0;
Vz{N+1} = ap_N;

cons = (tau >= 0);
acc_ek = 0;
for kl = 1 : N
    if kl < N
        if ~relax
            Vx{kl+1} = sdpvar(4,1);
            vx{kl+1} = sdpvar(1);
            Vz{kl+1} = sdpvar(1);
        elseif relax==1
            Vx{kl+1} = sdpvar(4,1);
            vx{kl+1} = sdpvar(1);
            Vz{kl+1} = L/2;
        elseif relax==2
            Vx{kl+1} = sdpvar(4,1);
            vx{kl+1} = sdpvar(1);
            Vz{kl+1} = L/2;
        elseif relax == 3
            Vx{kl+1} = [0; sdpvar(1); sdpvar(1); 0];
            vx{kl+1} = sdpvar(1);
            Vz{kl+1} = L/2;
        elseif relax == 4
            Vx{kl+1} = [0; 0; sdpvar(1); 0];
            vx{kl+1} = sdpvar(1);
            Vz{kl+1} = L/2;
        elseif relax == 5
            Vx{kl+1} = [0; 0; 0; 0];
            vx{kl+1} = sdpvar(1);
            Vz{kl+1} = L/2;
            %         Vx{kl+1} = [0 0 ; 0 sdpvar(1)];
            %         vx{kl+1} = sdpvar(1);
            %         Vx{kl+1} = [ 0 0; 0 -vx{kl+1}/2/L];
        end
    end
    
    lambda{kl}    = sdpvar(nbPts,nbPts,n,'full');
    lambdaLSY{kl} = sdpvar(2,1);
    lambdaLSX{kl} = sdpvar(2,1);
    lambdaVAR{kl} = sdpvar(1); acc_ek = acc_ek + lambdaVAR{kl};
    if relax == 2
        cons = cons + (lambdaLSX{kl}(2) == 0);
    end
    
    SStep{kl} = sdpvar(3); cons = cons + (SStep{kl} >= 0); cons = cons + (SStep{kl}(1,1) == Vz{kl+1});
    dsstep_avg = zeros(dimG);
    for i = 1:n
            dsstep = [yk1.' (zk-yk1).' gyk1(i,:).'] * SStep{kl} * [yk1.' (zk-yk1).' gyk1(i,:).'].'; % V2
            dsstep_avg = dsstep_avg + dsstep/n;
    end
    
    
    cons_SDP{kl} = - Vz{kl} * (zk-xs).'*(zk-xs) -lambdaVAR{kl} * var +  dsstep_avg ...
        - Vx{kl}(1)*term1_k-Vx{kl}(2)*term2_k-Vx{kl}(3)*term3_k-Vx{kl}(4)*term4_k...
        + Vx{kl+1}(1)*term1_k1 + Vx{kl+1}(2)*term2_k1 + Vx{kl+1}(3)*term3_k1 + Vx{kl+1}(4)*term4_k1;
    for i = 1:2
            cons_SDP{kl} = cons_SDP{kl} + lambdaLSY{kl}(i) * consLSY{i};
    end
    for i = 1:2
        for j = 1:n
            cons_SDP{kl} = cons_SDP{kl} + lambdaLSX{kl}(i) * consLS1{i,j}/n;
        end
    end
    cons_LIN{kl} = - vx{kl}.'*statesKf + vx{kl+1}.'*statesK1f;
    for k = 1:n
        clear XX FF GG;
        XX = { xs, xk, yk1}; FF = { fxs(k,:), fxk(k,:), fyk1(k,:)}; GG = { gxs(k,:), gxk(k,:), gyk1(k,:)};
        for i = 1:n
            XX{3+i} =  xk1(i, :);
            FF{3+i} = fxk1(k, :, i);
            GG{3+i} = gxk1(k, :, i);
        end
        for i = 1:nbPts
            for j = 1:nbPts
                if j ~= i
                    xi = XX{i}; xj = XX{j};
                    gi = GG{i}; gj = GG{j};
                    fi = FF{i}; fj = FF{j};
                    TT = [xi; xj; gi; gj];
                    
                    cons_SDP{kl} = cons_SDP{kl} + lambda{kl}(i,j,k) * TT.' * M * TT;
                    cons_LIN{kl} = cons_LIN{kl} + lambda{kl}(i,j,k) * (fi - fj);
                end
            end
        end
    end
    cons = cons + (lambdaVAR{kl} >=0);
    cons = cons + (cons_SDP{kl} <= 0);
    cons = cons + (cons_LIN{kl} == 0);
    cons = cons + (lambda{kl} >= 0);
end

obj = -acc_ek;

solver_opt = sdpsettings('solver','mosek','verbose',verbose,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',tolerance);
solverDetails=optimize(cons,-obj,solver_opt);


%% Try to grasp what happens...
if pplot
    close all;
    apk= zeros(1, N+1);
    Q1 = zeros(1, N+1);
    Q2 = zeros(1, N+1);
    Q3 = zeros(1, N+1);
    Q4 = zeros(1, N+1);
    dk = zeros(1, N+1);
    alphak = zeros(1, N+1);
    tauk = zeros(1, N+1);
    deltak = zeros(1, N+1);
    gammak = zeros(1, N+1);
    ek     = zeros(1, N+1);
    
    tolp = 1e-5;
    for i = 1:N+1
        apk(i) = double(Vz{i}(1,1)).';
        Q1(i) = double(Vx{i}(1));
        Q2(i) = double(Vx{i}(2));
        Q3(i) = double(Vx{i}(3));
        Q4(i) = double(Vx{i}(4));
        dk(i) = double(vx{i}(1));
    end
    for i = 1:N
        tauk(i)   = double(lambdaLSY{i}(2)/lambdaLSY{i}(1));
        alphak(i) = -double(lambdaLSX{i}(2)/lambdaLSX{i}(1));
        deltak(i) = double(SStep{i}(1,2)/SStep{i}(1,1));
        gammak(i) = -double(SStep{i}(1,3)/SStep{i}(1,1));
        ek(i)     = double(lambdaVAR{i});
    end
    figure;  hold on;
    subplot(11,1,1); 
    plot(1:N+1,apk,'-b'); title('apk')
    subplot(11,1,2); 
    plot(1:N+1,Q1,'-b'); title('Q1')
    subplot(11,1,3); 
    plot(1:N+1,Q2,'-b');title('Q2')
    subplot(11,1,4);
    plot(1:N+1,Q3,'-b');  title('Q3')
    subplot(11,1,5);
    plot(1:N+1,Q4,'-b'); title('Q4')
    subplot(11,1,6);
    plot(1:N+1,dk,'-b'); title('dk')
    subplot(11,1,7);
    plot(1:N+1,alphak,'-b'); title('alphak')
    subplot(11,1,8);
    plot(1:N+1,tauk,'-b'); title('tauk')
    subplot(11,1,9);
    plot(1:N+1,deltak,'-b'); title('deltak')
    subplot(11,1,10);
    plot(1:N+1,gammak,'-b'); title('gammak')
    subplot(11,1,10);
    plot(1:N+1,ek,'-b'); title('ek')
    
    if ssave
        labels{1} = 'k'; labels{2} = 'Q1'; labels{3} = 'Q2'; labels{4} = 'Q3'; labels{5} = 'Q4'; labels{6} = 'apk'; labels{7} = 'dk'; labels{8} = 'alphak'; labels{9} = 'tauk'; labels{10} = 'deltak'; labels{11} = 'gammak';
        data = [(1:N+1).' Q1.'  Q2.' Q3.' Q4.' apk.' dk.' alphak.' tauk.' deltak.' gammak.'];
        saveData([folder nname],data,labels);
    end
end












