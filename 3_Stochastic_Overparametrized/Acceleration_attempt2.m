clear all;
clc;

verbose = 1; tolerance = 1e-8;
ssave = 0;
pplot = 1;
folder = 'SaveData/';
nname  = 'FGM2_inter.dat';
relax  = 0;
% Problem
L = 1;
m = 0;
N = 100;
n = 2;

% P = [ xk zk | yk1 xk1 | f'(xk)  | f'(yk1) f'(xk1) ]
% F = [                    f(xk)  |  f(yk1)  f(xk1)  ]

dimG = 2 + 2*n + n^2; dimF =  3;
nbPts = n+2; % xs, y0, y1 ... yn

yk = zeros(1, dimG); yk(1) = 1;
zk = zeros(1, dimG); zk(2) = 1;
xs = zeros(1, dimG);
yk1= zeros(n, dimG); yk1(:,3:2+n) = eye(n);

gyk = zeros(n, dimG);gyk(:,3+n:2+2*n) = eye(n);
GYK = sum(gyk,1)/n;
gyk1= zeros(n, dimG, n);
GYK1= zeros(n, dimG);
gxs = zeros(n, dimG);

fyk = zeros(n, dimF); fyk(:,1) = 1;
fyk1= zeros(n, dimF, n);
FYK1= zeros(n, dimF);
fxs = zeros(n, dimF);

s_index = 3+2*n; e_index = s_index + n-1;
for i = 1:n
    gyk1(:,s_index:e_index,i) = eye(n);
    s_index = s_index + n; e_index = s_index + n-1;
    fyk1(:,3,i) = 1;
    fyk1(i,:,i) = [0 1 0];
    FYK1(i,:)   = sum(fyk1(:,:,i),1)/n;
    GYK1(i,:)   = sum(gyk1(:,:,i),1)/n;
end

% potentials
% term 1: || y - xs || ^2
% term 2: || f'(y)  || ^2
% term 3: || f'i(y) || ^2
% term 4: < f'(y); y-xs >

term1_k  = (yk-xs).'*(yk-xs);
term1_k1 = zeros(dimG);
for i = 1:n
    term1_k1 = term1_k1 + (yk1-xs).'*(yk1-xs)/n;
end

term2_k  = GYK.'*GYK;
term2_k1 = zeros(dimG);
for i = 1:n
    term2_k1 = term2_k1 + GYK1(i,:).'*GYK1(i,:)/n;
end

term3_k  = zeros(dimG);
for i = 1:n
    term3_k = term3_k + gyk(i,:).'*gyk(i,:)/n;
end

term3_k1 = zeros(dimG);
for i = 1:n
    for j = 1:n
        term3_k1 = term3_k1 + gyk1(i,:,j).'*gyk1(i,:,j)/n^2;
    end
end

term4_k  = (yk-xs).'*GYK; term4_k = 1/2 * (term4_k+term4_k.');
term4_k1 = zeros(dimG);
for i = 1:n
    temp     = (yk1(i,:)-xs).'*GYK1(i,:); temp = 1/2 * (temp+temp.');
    term4_k1 = term4_k1 + temp/n;
end

statesKf = (sum(fyk,1)/n - sum(fxs,1)/n);
statesK1f = (sum(FYK1,1)/n-sum(fxs,1)/n);
% line-searches

for i = 1:n
    consLS1{1,i} = GYK1(i,:).' * (yk1(i,:) - yk); consLS1{1,i} = (consLS1{1,i}+consLS1{1,i}.')/2;
    consLS1{2,i} = GYK1(i,:).' * (zk       - yk); consLS1{2,i} = (consLS1{2,i}+consLS1{2,i}.')/2;
    consLS1{3,i} = GYK1(i,:).' * (gyk(i,:));      consLS1{3,i} = (consLS1{3,i}+consLS1{3,i}.')/2;
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
d_N   = tau;

%%

M        = 1/2/(L-m) *[  -L*m, L*m,  m,   -L; L*m,  -L*m,  -m,    L;    m,    -m,  -1,  1 ;    -L,   L,   1,  -1];

Vx{1}   = Q0;
vx{1}   = d_0;
Vx{N+1} = QN;
vx{N+1} = d_N;
Vz{1}   = ap_0;
Vz{N+1} = ap_N;

cons = (tau >= 0);
for kl = 1 : N
    if kl < N
        if ~relax
            Vx{kl+1} = sdpvar(4,1);
            vx{kl+1} = sdpvar(1);
            Vz{kl+1} = sdpvar(1);
        elseif relax==1
            Vx{kl+1} = zeros(4,1);
            vx{kl+1} = sdpvar(1);
            Vz{kl+1} = L/2;
        elseif relax == 2
            Vx{kl+1} = [0; sdpvar(1); sdpvar(1); 0];
            vx{kl+1} = sdpvar(1);
            Vz{kl+1} = L/2;
        elseif relax == 3
            Vx{kl+1} = [0; 0; sdpvar(1); 0];
            vx{kl+1} = sdpvar(1);
            Vz{kl+1} = L/2;
            %         Vx{kl+1} = [0 0 ; 0 sdpvar(1)];
            %         vx{kl+1} = sdpvar(1);
            %         Vx{kl+1} = [ 0 0; 0 -vx{kl+1}/2/L];
        end
    end
    
    lambda{kl}    = sdpvar(nbPts,nbPts,n,'full');
    lambdaLSY{kl} = sdpvar(3,1);
    
    SStep{kl} = sdpvar(4); cons = cons + (SStep{kl} >= 0); cons = cons + (SStep{kl}(1,1) == Vz{kl+1});
    dsstep_avg = zeros(dimG);
    for i = 1:n
        for j = 1:n
            dsstep = [zk.' (yk1(i,:)-zk).' (yk-zk).' gyk1(i,:,j).'] * SStep{kl} * [zk.' (yk1(i,:)-zk).' (yk-zk).' gyk1(i,:,j).'].'; % V2
            dsstep_avg = dsstep_avg + dsstep/n^2;
        end
    end
    
    
    cons_SDP{kl} = - Vz{kl} * (zk-xs).'*(zk-xs) +  dsstep_avg ...
        - Vx{kl}(1)*term1_k-Vx{kl}(2)*term2_k-Vx{kl}(3)*term3_k-Vx{kl}(4)*term4_k...
        + Vx{kl+1}(1)*term1_k1 + Vx{kl+1}(2)*term2_k1 + Vx{kl+1}(3)*term3_k1 + Vx{kl+1}(4)*term4_k1;
    for i = 1:3
        for j = 1:n
            cons_SDP{kl} = cons_SDP{kl} + lambdaLSY{kl}(i)/n * consLS1{i,j};
        end
    end
    cons_LIN{kl} = - vx{kl}.'*statesKf + vx{kl+1}.'*statesK1f;
    for k = 1:n
        clear XX FF GG;
        XX = { xs, yk}; FF = { fxs(k,:), fyk(k,:)}; GG = { gxs(k,:), gyk(k,:)};
        for i = 1:n
            XX{2+i} =  yk1(i, :);
            FF{2+i} = fyk1(k, :, i);
            GG{2+i} = gyk1(k, :, i);
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
    cons = cons + (cons_SDP{kl} <= 0);
    cons = cons + (cons_LIN{kl} == 0);
    cons = cons + (lambda{kl} >= 0);
end
obj = tau;

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
    gammakp= zeros(1, N+1);
    
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
        alphak(i) = -double(lambdaLSY{i}(3)/lambdaLSY{i}(1));
        deltak(i) = double(SStep{i}(1,2)/SStep{i}(1,1));
        gammak(i) = -double(SStep{i}(1,3)/SStep{i}(1,1));
        gammakp(i)= -double(SStep{i}(1,4)/SStep{i}(1,1));
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
    subplot(11,1,11);
    plot(1:N+1,gammakp,'-b'); title('gammakp')
    
%     if ssave
%         labels{1} = 'k'; labels{2} = 'ak'; labels{3} = 'apk'; labels{4} = 'bk'; labels{5} = 'ck'; labels{6} = 'dk'; labels{7} = 'alphak'; labels{8} = 'tauk'; labels{9} = 'deltak'; labels{10} = 'gammak'; labels{11} = 'gammakp';
%         data = [(1:N+1).' ak.'  apk.' bk.' ck.' dk.' alphak.' tauk.' deltak.' gammak.' gammakp.'];
%         saveData([folder nname],data,labels);
%     end
end












