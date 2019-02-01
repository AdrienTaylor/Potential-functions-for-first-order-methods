clear all;
clc;

verbose = 1; tolerance = 1e-8;
ssave = 1;
pplot = 1;
folder = 'SaveData/';
nname  = 'FGM1_inter.dat';

% Problem
L = 1;
m = 0;
N = 100;

% P = [ xk zk | yk1 xk1 | f'(xk)  | f'(yk1) f'(xk1) ]
% F = [                    f(xk)  |  f(yk1)  f(xk1)  ]

dimG  = 7;
dimF  = 3;
nbPts = 4; % xs, xk,  yk1, xk1

% Notations for the LMI
xk  = zeros(1, dimG); xk(1)  = 1;
zk  = zeros(1, dimG); zk(2)  = 1;
yk1 = zeros(1, dimG); yk1(3) = 1;
xk1 = zeros(1, dimG); xk1(4) = 1;
xs  = zeros(1, dimG);

gxk = zeros(1, dimG); gxk(5) = 1;
gyk1= zeros(1, dimG);gyk1(6) = 1;
gxk1= zeros(1, dimG);gxk1(7) = 1;
gxs = zeros(1, dimG);

fxk = zeros(1, dimF); fxk(1) = 1;
fyk1= zeros(1, dimF);fyk1(2) = 1;
fxk1= zeros(1, dimF);fxk1(3) = 1;
fxs = zeros(1, dimF);

XX = {  xs,  xk,    yk1,    xk1};
GG = { gxs, gxk,   gyk1,   gxk1};
FF = { fxs, fxk,   fyk1,   fxk1};

% Potential function
statesKx  = [xk-xs; gxk];     statesKf  = fxk-fxs;
statesK1x = [xk1-xs; gxk1];   statesK1f = fxk1-fxs;


% line search for yk1
consLSY{1} = gyk1.'*(xk-yk1);   consLSY{1} = 1/2*(consLSY{1}+consLSY{1}.');
consLSY{2} = gyk1.'*(zk-xk);    consLSY{2} = 1/2*(consLSY{2}+consLSY{2}.');
% line search for xk1
consLSX{1} = gxk1.'*(yk1-xk1);  consLSX{1} = 1/2*(consLSX{1}+consLSX{1}.');
consLSX{2} = gxk1.'*(gyk1);     consLSX{2} = 1/2*(consLSX{2}+consLSX{2}.');


dimQ = 2; % in the same ordering: xk-x* | f'(xk) | zk-x* | f'(zk)
dimq = 1; % in the same ordering f(xk)-f*  | f(zk)-f*
Q0x = zeros(dimQ);  %Q0x(1,1) = L/2;
Q0z = zeros(1); Q0z(1,1) = L/2;
q0 = zeros(dimq,1); q0(1,1) = 0;

tau = sdpvar(1);
% QNx = zeros(dimQ);
QNx = [0 0; 0 tau/2/L];
% QNx(1,1) = 0;
QNz = zeros(1); QNz = L/2;
qN = tau;

%%

M        = 1/2/(L-m) *[  -L*m, L*m,  m,   -L; L*m,  -L*m,  -m,    L;    m,    -m,  -1,  1 ;    -L,   L,   1,  -1];

Vx{1}   = Q0x;
vx{1}   = q0;
Vx{N+1} = QNx;
vx{N+1} = qN;
Vz{1}   = Q0z;
Vz{N+1} = QNz;

cons = (tau >= 0);
for kl = 1 : N
    if kl < N
        %         Vx{kl+1} = sdpvar(2);
        %         vx{kl+1} = sdpvar(1);
        %         Vz{kl+1} = sdpvar(1);
        %         Vx{kl+1} = zeros(2);
        %         Vx{kl+1}(2,2) = sdpvar(1);
        %         Vx{kl+1} = [0 0 ; 0 sdpvar(1)];
        vx{kl+1} = sdpvar(1);
        Vx{kl+1} = [ 0 0; 0 vx{kl+1}/2/L];
        Vz{kl+1} = sdpvar(1);
    end
    
    lambda{kl}    = sdpvar(nbPts,nbPts,'full');
    lambdaLSY{kl} = sdpvar(2,1);
    lambdaLSX{kl} = sdpvar(2,1);
    
    
    SStep{kl} = sdpvar(3); cons = cons + (SStep{kl} >= 0); cons = cons + (SStep{kl}(1,1) == Vz{kl+1});
    cons_SDP{kl} = - Vz{kl} * (zk-xs).'*(zk-xs) - statesKx.'*Vx{kl}*statesKx +  [yk1; zk-yk1; gyk1].'*SStep{kl}*[yk1; zk-yk1; gyk1] + statesK1x.'*Vx{kl+1}*statesK1x;
    for i = 1:2
        cons_SDP{kl} = cons_SDP{kl} + lambdaLSY{kl}(i) * consLSY{i};
    end
    for i = 1:2
        cons_SDP{kl} = cons_SDP{kl} + lambdaLSX{kl}(i) * consLSX{i};
    end
    cons_LIN{kl} = - vx{kl}.'*statesKf + vx{kl+1}.'*statesK1f;
    
    for i = 1:nbPts
        for j = 1:nbPts
            if j ~= i
                xi = XX{i}; xj = XX{j};
                gi = GG{i}; gj = GG{j};
                fi = FF{i}; fj = FF{j};
                TT = [xi; xj; gi; gj];
                
                cons_SDP{kl} = cons_SDP{kl} + lambda{kl}(i,j) * TT.' * M * TT;
                cons_LIN{kl} = cons_LIN{kl} + lambda{kl}(i,j) * (fi - fj);
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
    ak = zeros(1, N+1);
    apk= zeros(1, N+1);
    bk = zeros(1, N+1);
    ck = zeros(1, N+1);
    dk = zeros(1, N+1);
    alphak = zeros(1, N+1);
    tauk = zeros(1, N+1);
    deltak = zeros(1, N+1);
    gammak = zeros(1, N+1);
    
    tolp = 1e-5;
    for i = 1:N+1
        ak(i)  = double(Vx{i}(1,1)).';
        apk(i) = double(Vz{i}(1,1)).';
        bk(i)  = double(Vx{i}(2,2)).';
        ck(i)  = double(Vx{i}(2,1)).';
        dk(i)  = double(vx{i}(1,1)).';
    end
    for i = 1:N
        alphak(i) = -double(lambdaLSX{i}(2)/lambdaLSX{i}(1));
        tauk(i)   = double(lambdaLSY{i}(2)/lambdaLSY{i}(1));
        deltak(i) = double(SStep{i}(1,2)/SStep{i}(1,1));
        gammak(i) = -double(SStep{i}(1,3)/SStep{i}(1,1));
    end
    
    subplot(9,1,1);
    plot(1:N+1,ak,'-b'); hold on
    subplot(9,1,2);
    plot(1:N+1,apk,'-b');
    subplot(9,1,3);
    plot(1:N+1,bk,'-b');
    subplot(9,1,4);
    plot(1:N+1,ck,'-b');
    subplot(9,1,5);
    plot(1:N+1,dk,'-b');
    subplot(9,1,6);
    plot(1:N+1,alphak,'-b');
    subplot(9,1,7);
    plot(1:N+1,tauk,'-b');
    subplot(9,1,8);
    plot(1:N+1,deltak,'-b');
    subplot(9,1,9);
    plot(1:N+1,gammak,'-b');
    
    if ssave
        labels{1} = 'k'; labels{2} = 'ak'; labels{3} = 'apk'; labels{4} = 'bk'; labels{5} = 'ck'; labels{6} = 'dk'; labels{7} = 'alphak'; labels{8} = 'tauk'; labels{9} = 'deltak'; labels{10} = 'gammak';
        data = [(1:N+1).' ak.'  apk.' bk.' ck.' dk.' alphak.' tauk.' deltak.' gammak.'];
        saveData([folder nname],data,labels);
    end
end












