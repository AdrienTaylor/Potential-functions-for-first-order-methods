clear all;
clc;
tic;
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

% POTENTIAL SETUP
% Potential has the form:
%   q1k ||x_k-x*||^2 + q2k ||f'(x_k)||^2 + 2 q3k E_i||f_i'(xk)||^2 
%       + q4k <f'(x_k); x_k-x*>  + dk (f(x_k)-f(x*)) + ak ||z_k-x*||^2

% INTERMEDIARY POTENTIAL SETUPS (aims at reproducing results from
% Section 4 and Appendix E) via option "relaxation".

% Options:
%   Set relax = 0: use the base potential stated in Appendix E.2 (Fig 5).
%   Set relax = 1: force ak = L/2.
%   Set relax = 2: force ak = L/2 and alphak=0 (in the method)
%   Set relax = 3: force ak = L/2, q1k=q2k=q3k=q4k=0 and alphak=0 (in the method)

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
elseif relax == 3
    a0 = L/2; aN = L/2;
    q1N = 0;  q10 = 0;
    q2N = 0;  q20 = 0;
    q3N = 0;  q30 = 0;
    q4N = 0;  q40 = 0;
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


% Matrix encoding interpolation condition for smooth strongly convex
% functions
M = 1/2/(L-m) *[   -L*m,  L*m,   m, -L;
                    L*m, -L*m,  -m,  L;
                      m,   -m,  -1,  1;
                     -L,   L,   1,  -1];
                 
% Notations for potentials
statesKF = (FXK - FXS); statesK1F = (sum(FXK1,1)/n-FXS);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              %
%          SETTING UP THE LINEAR MATRIX INEQUALITIES           %
%     (directly embedded within the larger problem (2))        %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Q{1}   = [q10 q20 q30 q40].';
d{1}   = d0;
a{1}   = a0;
Q{N+1} = [q1N q2N q3N q4N].';
d{N+1} = dN;
a{N+1} = aN;

cons = (tau >= 0);
for k = 1 : N % iteration counter (one LMI per iteration)
    if k < N
        if ~relax % no additional constraints on potentials
            Q{k+1} = sdpvar(4,1);
            d{k+1} = sdpvar(1);
            a{k+1} = sdpvar(1);
        elseif relax == 1 %  ak = L/2
            Q{k+1} = sdpvar(4,1);
            d{k+1} = sdpvar(1);
            a{k+1} = L/2;
        elseif relax == 2 % force ak = L/2 and alphak=0 (in the method)
            Q{k+1} = sdpvar(4,1);
            d{k+1} = sdpvar(1);
            a{k+1} = L/2;
        elseif relax == 3 % force ak = L/2, q1k=q2k=q3k=q4k=0 and alphak=0 (in the method)
            Q{k+1} = zeros(4,1);
            d{k+1} = sdpvar(1);
            a{k+1} = L/2;
        end
    end
    lambda{k}   = sdpvar(nbPts,nbPts,n,'full');
    mu_Y{k}     = sdpvar(2,1); % multipliers for the line-search for y_{k+1}
    mu_X{k}     = sdpvar(2,1); % multipliers for the line-search for x_{k+1}

    if relax == 2
        cons = cons + (mu_X{k}(2) == 0);
    end
    
    S{k} = sdpvar(3); % this is Sk (see Section C.5)
    cons = cons + (S{k} >= 0);
    cons = cons + (S{k}(1,1) == a{k+1});
    
    dsstep_avg = zeros(dimG); % Compute E_i|| z_{k+1}^{(i)}-x_* ||
    for i = 1:n
            % this is the term || z_{k+1}^{(i)}-x_* ||
            dsstep = [yk1.' (zk-yk1).' gyk1(i,:).'] * S{k} * [yk1.' (zk-yk1).' gyk1(i,:).'].'; 
            dsstep_avg = dsstep_avg + dsstep/n; % this is the average
    end
    
    
    cons_SDP{k} = - a{k} * (zk-xs).'*(zk-xs) +  dsstep_avg ...
        - Q{k}(1)*term1_k-Q{k}(2)*term2_k-Q{k}(3)*term3_k-Q{k}(4)*term4_k...
        + Q{k+1}(1)*term1_k1 + Q{k+1}(2)*term2_k1 + Q{k+1}(3)*term3_k1 + Q{k+1}(4)*term4_k1;
    for i = 1:2
            cons_SDP{k} = cons_SDP{k} + mu_Y{k}(i) * A{i};
    end
    for i = 1:2
        for j = 1:n
            cons_SDP{k} = cons_SDP{k} + mu_X{k}(i) * B{i,j}/n;
        end
    end
    cons_LIN{k} = - d{k}.'*statesKF + d{k+1}.'*statesK1F;
    for l = 1:n
    % POINTS IN THE DISCRETE REPRESENTATION OF FUNCTION f_l(.)
    clear X G F;
        X = { xs, xk, yk1}; 
        G = { gxs(l,:), gxk(l,:), gyk1(l,:)};
        F = { fxs(l,:), fxk(l,:), fyk1(l,:)};
        for i = 1:n
            X{3+i} =  xk1(i, :);
            F{3+i} = fxk1(l, :, i);
            G{3+i} = gxk1(l, :, i);
        end
        for i = 1:nbPts
            for j = 1:nbPts
                if j ~= i
                    xi = X{i}; xj = X{j};
                    gi = G{i}; gj = G{j};
                    fi = F{i}; fj = F{j};
                    TT = [xi; xj; gi; gj];
                    
                    cons_SDP{k} = cons_SDP{k} + lambda{k}(i,j,l) * TT.' * M * TT;
                    cons_LIN{k} = cons_LIN{k} + lambda{k}(i,j,l) * (fi - fj);
                end
            end
        end
    end
    cons = cons + (cons_SDP{k} <= 0);
    cons = cons + (cons_LIN{k} == 0);
    cons = cons + (lambda{k} >= 0);
end
obj = tau;

solver_opt = sdpsettings('solver','mosek','verbose',verbose,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',tolerance);
solverDetails=optimize(cons,-obj,solver_opt);
toc

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
        apk(i)= double(a{i}(1,1)).';
        Q1(i) = double(Q{i}(1));
        Q2(i) = double(Q{i}(2));
        Q3(i) = double(Q{i}(3));
        Q4(i) = double(Q{i}(4));
        dk(i) = double(d{i}(1));
    end
    for i = 1:N
        tauk(i)   = double(mu_Y{i}(2)/mu_Y{i}(1));
        alphak(i) = -double(mu_X{i}(2)/mu_X{i}(1));
        deltak(i) = double(S{i}(1,2)/S{i}(1,1));
        gammak(i) = -double(S{i}(1,3)/S{i}(1,1));
    end
    figure;  hold on;
    subplot(4,3,1); 
    plot(1:N+1,Q1,'-b'); title('q1k (coefficient of ||xk-x*||^2)')
    subplot(4,3,4); 
    plot(1:N+1,Q2,'-b');title('q2k (coefficient of ||f''(xk)||^2)')
    subplot(4,3,5);
    plot(1:N+1,Q3,'-b');  title('q3k (coefficient of Ei||fi''(xk)||^2)')
    subplot(4,3,6);
    plot(1:N+1,Q4,'-b'); title('q4k (coefficient of <f''(xk); xk-x*>)')
    subplot(4,3,7);
    plot(1:N+1,dk,'-b'); title('dk (coefficient of f(xk)-f(x*))')
    subplot(4,3,8); 
    plot(1:N+1,apk,'-b'); title('ak (coefficient of ||zk-x*||^2)')
    subplot(4,3,9);
    plot(1:N,tauk(1:end-1),'-b'); title('tauk') 
    subplot(4,3,10);
    plot(1:N,alphak(1:end-1),'-b'); title('alphak')
    subplot(4,3,11);
    plot(1:N,deltak(1:end-1),'-b'); title('deltak')
    subplot(4,3,12);
    plot(1:N,gammak(1:end-1),'-b'); title('gammak')
    
    if ssave
        labels{1} = 'k'; labels{2} = 'Q1'; labels{3} = 'Q2'; labels{4} = 'Q3'; labels{5} = 'Q4'; labels{6} = 'apk'; labels{7} = 'dk'; labels{8} = 'alphak'; labels{9} = 'tauk'; labels{10} = 'deltak'; labels{11} = 'gammak';
        data = [(1:N+1).' Q1.'  Q2.' Q3.' Q4.' apk.' dk.' alphak.' tauk.' deltak.' gammak.'];
        saveData([folder nname],data,labels);
    end
end












