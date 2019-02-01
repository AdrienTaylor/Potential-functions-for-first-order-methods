clear all;
clc;

% SOLVER OPTIONS
verbose     = 1;
tolerance   = 1e-8;

% OUTPUT OPTIONS                   
ssave   = 1;            % Save the results ?
pplot   = 1;            % Plot the results ?
folder  = 'SaveData/';  % If results saved, name of the saving folder
nname   = 'ParameterSelection.dat';    % Name of the file

% MINIMIZATION PROBLEM SETUP
L = 1;      % Smoothness constant
m = 0;      % Strong convexity constant
n = 2;      % Cardinality of the support for the stochastic gradient
            % (i.e., number of functions in the finite sum)
N = 100;    % Number of iterations

% ALGORITHM SETUP: stochstic gradient descent with step-size policy delta_k.
alpha = 0;
delta = @(k)(1/L/(k+1)^(alpha));   % Step-size (possibly varying function of k)

% POTENTIAL SETUP
% Potential has the form:
%   q1k ||x_k-x*||^2 + q2k ||f'(x_k)||^2 + 2 q3k E_i||f_i'(xk)||^2 
%       + q4k <f'(x_k); x_k-x*>  + dk (f(x_k)-f(x*)) + ak ||z_k-x*||^2

% INTERMEDIARY POTENTIAL SETUPS (aims at reproducing results from
% Section 3.2.2) via option "relaxation".

% Options:
%   Set relax = 0: use the base potential stated in Section 3.2.2.
%   Set relax = 1: force ak = L/2.
%   Set relax = 2: force ak = L/2 and alphak=0 (in the method)

relax = 0;

% INITIAL AND FINAL POTENTIALS SETUP:

a0  = L/2;
q10 = 0;
q20 = 0;
q30 = 0;
q40 = 0;
d0  = 0;

tau = sdpvar(1); % variable to be maximized (see problem (3) in the paper)

aN  = 0;
q1N = 0;
q2N = 0;
q3N = 0;
q4N = 0;
dN  = tau;

if relax == 1
    a0 = L/2; aN = L/2;
elseif relax == 2
    a0 = L/2; aN = L/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              %
% SETTING UP THE NOTATIONS FOR THE LINEAR MATRIX INEQUALITIES  %
%                   (end of editable zone)                     %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recall that x* (optimum) is set to x* = 0 without loss of generality,
% and so does f(x*) = 0. We also have fi'(x*) = 0 by assumption
% (over-parametrized models)


% P = [ xk zk yk1 | xk1 ... xkn | f1'(x_k) ... fn'(x_k) | f1'(y_{k+1}) ... fn'(y_{k+1}) | f1'(x_{k+1}^(1)) ... fn'(x_{k+1}^(1)) | ... | f1'(x_{k+1}^(n)) ... fn'(x_{k+1}^(n))]
% F = [                           f1(x_k)  ... fn(x_k)  | f1(y_{k+1})  ... fn(y_{k+1})  | f1(x_{k+1}^(1))  ... fn(x_{k+1}^(1))  | ... | f1(x_{k+1}^(n))  ... fn(x_{k+1}^(n))]
% NOTE: Symmetry arguments allow simplifying both F and G=P.'*P.

dimG  = 3 + 3*n + n^2;  % dimensions of the Gram matrix P.'*P
dimF  = 2*n + n^2;      % dimensions of F
nbPts = n+3;    % number of points to be incorporated in the discrete
% version of each function f_i:
% x*,y_{k},x_{k},x_{k+1}^{(1)},...,x_{k+1}^{(n)}

xk = zeros(1, dimG); xk(1) = 1; % this is x_k
zk = zeros(1, dimG); zk(2) = 1; % this is z_k
yk1= zeros(1, dimG); yk1(3)= 1; % this is y_{k+1}
xs = zeros(1, dimG);            % this is x*
xk1= zeros(n, dimG); xk1(:,4:3+n) = eye(n); % xk1(i,:) is x_{k+1}^{(i)}

gxk = zeros(n, dimG);gxk(:,4+n:3+2*n) = eye(n); % gxk(i,:) is fi'(x_k)
GXK = sum(gxk,1)/n;                             % GXK is f'(x_k)
gyk1 = zeros(n, dimG);gyk1(:,4+2*n:3+3*n) = eye(n); % gyk(i,:) is fi'(y_{k+1})
GYK1 = sum(gyk1,1)/n;                               % GYK1 is f'(y_{k+1})
gxk1= zeros(n, dimG, n);    % gxk1(i,:,j) is fi'(x_{k+1}^(j))
GXK1= zeros(n, dimG);       % GXK1(j,:)   is f'(x_{k+1}^(j))
gxs = zeros(n, dimG);       % gxs is fi'(x*)
GXS = sum(gxs,1)/n;         % GXS is f(x*)

fxk = zeros(n, dimF); fxk(:,1:n) = eye(n);      % fxk(i,:) is fi(x_k)
FXK = sum(fxk,1)/n;                             % FXK is f(x_k)
fyk1 = zeros(n, dimF); fyk1(:,1+n:2*n) = eye(n);% fyk1(i,:) is fi(y_{k+1})
FYK1 = sum(fyk1,1)/n;                           % FYK1 is f(y_{k+1})
fxk1= zeros(n, dimF, n);  % fxk1(i,:,j) is fi(x_{k+1}^(j))
FXK1= zeros(n, dimF);     % FXK1(j,:)  is f(x_{k+1}^(j))
fxs = zeros(n, dimF);     % fxs(i,:)  is fi(x*)
FXS = sum(fxs,1)/n;       % FXS is f(x*)

Gs_index = 4+3*n; Ge_index = Gs_index + n-1;
Fs_index = 1+2*n; Fe_index = Fs_index + n-1;
for i = 1:n
    gxk1(:,Gs_index:Ge_index,i) = eye(n);
    Gs_index = Gs_index + n; Ge_index = Gs_index + n-1;
    fxk1(:,Fs_index:Fe_index,i) = eye(n);
    Fs_index = Fs_index + n; Fe_index = Fs_index + n-1;
    FXK1(i,:)   = sum(fxk1(:,:,i),1)/n;
    GXK1(i,:)   = sum(gxk1(:,:,i),1)/n;
end

% TERMS IN THE POTENTIAL
% Potential has the form:
%   q1k ||x_k-x*||^2 + q2k ||f'(x_k)||^2 + 2 q3k E_i||f_i'(xk)||^2 
%       + q4k <f'(x_k); x_k-x*>  + dk (f(x_k)-f(x*)) + ak ||z_k-x*||^2

% potential part in terms of x's
% term 1:    ||x-xs|| ^2
% term 2:    ||f'(x)|| ^2
% term 3: E_i||f'i(x)|| ^2
% term 4: < f'(x); x-xs >

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
    term3_k = term3_k + gxk(i,:).'*gxk(i,:)/n;
end

term3_k1 = zeros(dimG);
for i = 1:n
    for j = 1:n
        term3_k1 = term3_k1 + gxk1(i,:,j).'*gxk1(i,:,j)/n^2;
    end
end

term4_k  = (xk-xs).'*GXK; term4_k = 1/2 * (term4_k+term4_k.');
term4_k1 = zeros(dimG);
for i = 1:n
    temp     = (xk1(i,:)-xs).'*GXK1(i,:); temp = 1/2 * (temp+temp.');
    term4_k1 = term4_k1 + temp/n;
end

statesKf = (FXK - FXS); statesK1f = (sum(FXK1,1)/n-FXS);

% line-searches for finding y_{k+1}
A{1} = GYK1.' * (xk - yk1); A{1} = (A{1}+A{1}.')/2;
A{2} = GYK1.' * (zk -  xk); A{2} = (A{2}+A{2}.')/2;
% line-searches for finding x_{k+1}^{(i)} (we average them, using
% appropriate symmetry arguments (no reason for one of them to be
% different).
for i = 1:n
    B{1,i} = GXK1(i,:).' * (yk1 - xk1(i,:)); B{1,i} = (B{1,i}+B{1,i}.')/2;
    B{2,i} = GXK1(i,:).' * (gyk1(i,:));      B{2,i} = (B{2,i}+B{2,i}.')/2;
end

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
            %         Vx{kl+1} = [0 0 ; 0 sdpvar(1)];
            %         vx{kl+1} = sdpvar(1);
            %         Vx{kl+1} = [ 0 0; 0 -vx{kl+1}/2/L];
        end
    end
    
    lambda{kl}    = sdpvar(nbPts,nbPts,n,'full');
    lambdaLSY{kl} = sdpvar(2,1);
    lambdaLSX{kl} = sdpvar(2,1);
    if relax == 2
        cons = cons + (lambdaLSX{kl}(2) == 0);
    end
    
    SStep{kl} = sdpvar(3); cons = cons + (SStep{kl} >= 0); cons = cons + (SStep{kl}(1,1) == Vz{kl+1});
    dsstep_avg = zeros(dimG);
    for i = 1:n
            dsstep = [yk1.' (zk-yk1).' gyk1(i,:).'] * SStep{kl} * [yk1.' (zk-yk1).' gyk1(i,:).'].'; % V2
            dsstep_avg = dsstep_avg + dsstep/n;
    end
    
    
    cons_SDP{kl} = - Vz{kl} * (zk-xs).'*(zk-xs) +  dsstep_avg ...
        - Vx{kl}(1)*term1_k-Vx{kl}(2)*term2_k-Vx{kl}(3)*term3_k-Vx{kl}(4)*term4_k...
        + Vx{kl+1}(1)*term1_k1 + Vx{kl+1}(2)*term2_k1 + Vx{kl+1}(3)*term3_k1 + Vx{kl+1}(4)*term4_k1;
    for i = 1:2
            cons_SDP{kl} = cons_SDP{kl} + lambdaLSY{kl}(i) * A{i};
    end
    for i = 1:2
        for j = 1:n
            cons_SDP{kl} = cons_SDP{kl} + lambdaLSX{kl}(i) * B{i,j}/n;
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
    end
    figure;  hold on;
    subplot(10,1,1); 
    plot(1:N+1,apk,'-b'); title('apk')
    subplot(10,1,2); 
    plot(1:N+1,Q1,'-b'); title('Q1')
    subplot(10,1,3); 
    plot(1:N+1,Q2,'-b');title('Q2')
    subplot(10,1,4);
    plot(1:N+1,Q3,'-b');  title('Q3')
    subplot(10,1,5);
    plot(1:N+1,Q4,'-b'); title('Q4')
    subplot(10,1,6);
    plot(1:N+1,dk,'-b'); title('dk')
    subplot(10,1,7);
    plot(1:N+1,alphak,'-b'); title('alphak')
    subplot(10,1,8);
    plot(1:N+1,tauk,'-b'); title('tauk')
    subplot(10,1,9);
    plot(1:N+1,deltak,'-b'); title('deltak')
    subplot(10,1,10);
    plot(1:N+1,gammak,'-b'); title('gammak')
    
    if ssave
        labels{1} = 'k'; labels{2} = 'Q1'; labels{3} = 'Q2'; labels{4} = 'Q3'; labels{5} = 'Q4'; labels{6} = 'apk'; labels{7} = 'dk'; labels{8} = 'alphak'; labels{9} = 'tauk'; labels{10} = 'deltak'; labels{11} = 'gammak';
        data = [(1:N+1).' Q1.'  Q2.' Q3.' Q4.' apk.' dk.' alphak.' tauk.' deltak.' gammak.'];
        saveData([folder nname],data,labels);
    end
end












