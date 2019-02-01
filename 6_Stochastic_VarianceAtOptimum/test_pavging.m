clear all;
clc;

verbose = 1; tolerance = 1e-8;


% Problem
L = 1;
m = 0;
n = 3;
dk= 1;
dk1=1+1/2;

% P = [ yk zk yk1 | f1'(yk) ... fn'(yk) | f1'(yk1) ... fn'(ykn) | f1'(xs) ... fn'(xs)]
% F = [                    fi(yk)  |  fi(yk1)   ]

dimG = 3 + 2*n + (n-1); dimF = 2;
nbPts = 3; % xs, yk, yk1

yk = zeros(1, dimG); yk(1) = 1;
zk = zeros(1, dimG); zk(2)= 1;
yk1= zeros(1, dimG); yk1(3)= 1;
xs = zeros(1, dimG);

gyk = zeros(n, dimG);gyk(:,4:3+n) = eye(n);
GYK = sum(gyk,1)/n;
gyk1 = zeros(n, dimG);gyk1(:,4+n:3+2*n) = eye(n);
GYK1 = sum(gyk1,1)/n;
gxs = zeros(n, dimG); gxs(1:n-1,dimG-(n-2):dimG) = eye(n-1);
gxs(end,:) = - sum(gxs,1);
GXS = sum(gxs,1)/n;

fyk = zeros(n, dimF); fyk(:,1) = 1;
FYK = sum(fyk,1)/n;
fyk1 = zeros(n, dimF); fyk1(:,2) = 1;
FYK1 = sum(fyk1,1)/n;
fxs = zeros(n, dimF);
FXS = sum(fxs,1)/n;

statesKf = (FYK - FXS);
statesK1f = (sum(FYK1,1)-FXS);
% line-searches

consLSY{1} = GYK1.' * (yk - yk1); consLSY{1} = (consLSY{1}+consLSY{1}.')/2;
consLSY{2} = GYK1.' * (zk -  yk); consLSY{2} = (consLSY{2}+consLSY{2}.')/2;

var = zeros(dimG);
for i = 1:n
    var = var + (gxs(i,:)-GXS).'*(gxs(i,:)-GXS)/n;
end

%%

M        = 1/2/(L-m) *[  -L*m, L*m,  m,   -L; L*m,  -L*m,  -m,    L;    m,    -m,  -1,  1 ;    -L,   L,   1,  -1];
lambda    = sdpvar(nbPts,nbPts,'full');
lambdaLSY = sdpvar(2,1);
lambdaVAR = sdpvar(1);
tau       = sdpvar(1);



SStep = sdpvar(2); cons =  (SStep >= 0); cons = cons + (SStep(1,1) == 1);
cons = cons + (tau==dk1);
dsstep_avg = zeros(dimG);
for i = 1:n
    dsstep = [zk.' gyk1(i,:).'] * SStep * [zk.' gyk1(i,:).'].'; % V2
    dsstep_avg = dsstep_avg + dsstep/n;
end


cons_SDP = - L/2 * (zk-xs).'*(zk-xs) -lambdaVAR * var + L/2 * dsstep_avg;
for i = 1:2
    cons_SDP = cons_SDP + lambdaLSY(i) * consLSY{i};
end

cons_LIN = - dk*statesKf + tau*statesK1f;
for k = 1:n
    clear XX FF GG;
    XX = { xs, yk, yk1}; FF = { fxs(k,:), fyk(k,:), fyk1(k,:)}; GG = { gxs(k,:), gyk(k,:), gyk1(k,:)};
    listi = [1 2];
    listj = [3 3];
    for ll = 1:length(listi)
        i = listi(ll); j = listj(ll);
        if j ~= i
            xi = XX{i}; xj = XX{j};
            gi = GG{i}; gj = GG{j};
            fi = FF{i}; fj = FF{j};
            TT = [xi; xj; gi; gj];
            
            if ll == 1
%                 
%                 matt = gj.'*(xj-xi);
%                 matt = 1/2 * (matt + matt.');
%                 %             matt = gi.'*(xi-xj)+L/2*(xi-xj).'*(xi-xj);
%                 %             matt = 1/2 * (matt + matt.');
cons_SDP = cons_SDP + lambda(i,j) * TT.' * M * TT;
%                 cons_SDP = cons_SDP + lambda(i,j) * matt;
                cons_LIN = cons_LIN + lambda(i,j) * (fi - fj);
            elseif ll == 2
%                 
                matt = gj.'*(xj-xi);
                matt = 1/2 * (matt + matt.');
%                 %             matt = gi.'*(xi-xj)+L/2*(xi-xj).'*(xi-xj);
%                 %             matt = 1/2 * (matt + matt.');
%                             cons_SDP = cons_SDP + lambda(i,j) * TT.' * M * TT;
                cons_SDP = cons_SDP + lambda(i,j) * matt;
                cons_LIN = cons_LIN + lambda(i,j) * (fi - fj);
            end
%                             cons_SDP = cons_SDP + lambda(i,j) * TT.' * M * TT;
%                             cons_LIN = cons_LIN + lambda(i,j) * (fi - fj);
        end
    end
end
cons = cons + (lambdaVAR >=0);
cons = cons + (cons_SDP <= 0);
cons = cons + (cons_LIN == 0);
cons = cons + (lambda >= 0);
obj = -lambdaVAR;

solver_opt = sdpsettings('solver','mosek','verbose',verbose,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',tolerance);
solverDetails=optimize(cons,-obj,solver_opt);


