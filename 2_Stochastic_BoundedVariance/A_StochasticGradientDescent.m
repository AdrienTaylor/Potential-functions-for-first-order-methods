clear all;
clc;

% SOLVER OPTIONS
verbose     = 1;
tolerance   = 1e-8;

% OUTPUT OPTIONS
reoptimize  = 1;        % this option should be set to 1 for the two-step
                        % optimization procedure (see (3) in the paper)
tol_reopt   = 1e-5;     % tolerance in the re-optimization process                        
ssave   = 1;            % Save the results ?
pplot   = 1;            % Plot the results ?
folder  = 'SaveData/';  % If results saved, name of the saving folder
nname   = 'SGD.dat';    % Name of the file

% MINIMIZATION PROBLEM SETUP
L = 1;      % Smoothness constant
m = 0;      % Strong convexity constant
n = 2;      % Cardinality of the support for the stochastic gradient
N = 100;    % Number of iterations

% ALGORITHM SETUP: stochstic gradient descent with step-size policy delta_k.
alpha = 0;
delta = @(k)(1/L/(k+1)^(alpha));   % Step-size (possibly varying function of k)

% POTENTIAL SETUP
% Potential has the form:
%   ak ||x_k-x*||^2 + bk ||f'(x_k)||^2 + 2 ck <f'(x_k); x_k-x*>
%       + dk (f(x_k)-f(x*)).

% INTERMEDIARY POTENTIAL SETUPS (aims at reproducing results from
% Section 3.2.2) via option "relaxation".

% Options:
%   Set relax = 0: use the base potential stated in Section 3.2.2.
%   Set relax = 1: force ck = 0.
%   Set relax = 2: force bk = ck = 0.
%   Set relax = 3: force bk = ck = 0, and ak = L/2.

relax = 0;

% INITIAL AND FINAL POTENTIALS SETUP:

a0 = L/2;
b0 = 0;
c0 = 0;
d0 = 0;

tau = sdpvar(1); % variable to be maximized (see problem (3) in the paper)

aN  = 0;
bN  = 0; % bN is to be maximized
cN  = 0;
dN  = tau;

if relax == 1
    c0 = 0; cN = 0;
elseif relax == 2
    b0 = 0; bN = 0;
    c0 = 0; cN = 0;
elseif relax == 3
    a0 = L/2 ; aN = L/2;
    b0 = 0; bN = 0;
    c0 = 0; cN = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              %
% SETTING UP THE NOTATIONS FOR THE LINEAR MATRIX INEQUALITIES  %
%                   (end of editable zone)                     %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recall that x* (optimum) is set to x* = 0 without loss of generality,
% and so does f(x*) = 0.

% P = [ xk  G(x_k;1) ... G(x_k;n) f'(x_{k+1}^{(1)}) ... f'(x_{k+1}^{(n)})]
% F = [           f(x_k)          f(x_{k+1})  ...       f(x_{k+1}^{(n)})   ]
% NOTE: Symmetry arguments allow simplifying both F and G=P.'*P.

dimG  = 1+2*n;  % dimensions of the Gram matrix P.'*P
dimF  = n+1;    % dimensions of F
nbPts = 2+n;    % number of points to be incorporated in the discrete
% version of the function f:
% x*,x_{k},x_{k+1}^{(1)},...,x_{k+1}^{(n)}

% Notations for the LMI
xk  = zeros(1, dimG); xk(1)  = 1; % this is x_{k}
xs  = zeros(1, dimG);             % this is x* = 0

Gxk = zeros(n, dimG); Gxk(:,2:n+1) = eye(n);     % Gxk(i,:) is G(xk;i)
gxk = sum(Gxk,1)/n;                              % this is f'(x_{k})
gxk1= zeros(n, dimG);gxk1(:,2+n:1+2*n) = eye(n); % gxk1(i,:) is f'(x_{k+1}^{(i)})
gxs = zeros(1, dimG);                            % this is f'(x*) = 0

fxk = zeros(1, dimF); fxk(1) = 1;               % this is f(x_{k})
fxk1= zeros(n, dimF);fxk1(:,2:n+1) = eye(n);    % fxk1(i,:) is f(x_{k+1}^{(i)})
fxs = zeros(1, dimF);                           % this is f(x*) = 0

xk1 = @(k)(ones(n,1)*xk - delta(k) * Gxk);      % xk1(i,:) is x_{k+1}^{(i)}

% Matrix encoding interpolation condition for smooth strongly convex
% functions
M = 1/2/(L-m) *[   -L*m,  L*m,   m, -L;
                    L*m, -L*m,  -m,  L;
                      m,   -m,  -1,  1;
                     -L,   L,   1,  -1];
% Matrix encoding the variance condition
Avar = zeros(dimG);
for i = 1:n
    Avar = Avar + (Gxk(i,:)-gxk).'*(Gxk(i,:)-gxk)/n;
end
% Notations for potentials
Q0 = [ a0 c0; c0 b0]; q0 = d0;
QN = [ aN cN; cN bN]; qN = dN;

statesKG  = [(xk-xs).' gxk.'];  statesKF = fxk;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              %
%          SETTING UP THE LINEAR MATRIX INEQUALITIES           %
%     (directly embedded within the larger problem (2))        %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q{1}   = Q0;
d{1}   = q0;
Q{N+1} = QN;
d{N+1} = qN;

cons   = (tau>=0);
for k = 1 : N % iteration counter (one LMI per iteration)
    if k < N
        if ~relax % no constraint on Qk and dk
            d{k+1} = sdpvar(1);
            Q{k+1} = sdpvar(2);
        elseif relax == 1 % force ck = 0
            d{k+1} = sdpvar(1);
            Q{k+1} = [sdpvar(1) 0; 0 sdpvar(1)];
        elseif relax == 2 % force bk = ck = 0
            d{k+1} = sdpvar(1);
            Q{k+1} = [sdpvar(1) 0; 0 0];
        elseif relax == 3 % force bk = ck = 0, and ak = L/2
            d{k+1} = sdpvar(1);
            Q{k+1} = [L/2 0; 0 0];
        end
    end
    
    % POINTS IN THE DISCRETE REPRESENTATION OF FUNCTION f(.)
    clear X G F;
    XK1 = xk1(k);      % xk1(i,:) is x_{k+1}^{(i)}
    X = {  xs,    xk}; % coordinates
    G = { gxs,   gxk}; % gradients
    F = { fxs,   fxk}; % function values
    for i = 1:n
        X{2+i}  = XK1(i,:);
        G{2+i}  = gxk1(i,:);
        F{2+i}  = fxk1(i,:);
    end
    
    
    lambda{k}    = sdpvar(nbPts,nbPts,'full');
    e{k}         = sdpvar(1);
    cons_SDP{k}  = - statesKG * Q{k} * statesKG.' - e{k} * Avar;
    cons_LIN{k}  = - d{k} * statesKF;
    
    for i = 1:n % averaging the states over i of x_{k+1}^{(i)}
        statesK1G = [(XK1(i,:)-xs).' gxk1(i,:).']; statesK1F = fxk1(i,:);
        cons_SDP{k}= cons_SDP{k} + statesK1G * Q{k+1} * statesK1G.'/n;
        cons_LIN{k}= cons_LIN{k} + statesK1F * d{k+1}/n;
    end
    
    for i = 1:nbPts
        for j = 1:nbPts
            if j ~= i
                xi = X{i}; xj = X{j};
                gi = G{i}; gj = G{j};
                fi = F{i}; fj = F{j};
                TT = [xi; xj; gi; gj];
                
                cons_SDP{k} = cons_SDP{k} + lambda{k}(i,j) * TT.' * M * TT;
                cons_LIN{k} = cons_LIN{k} + lambda{k}(i,j) * (fi - fj);
            end
        end
    end
    
    cons = cons + (cons_SDP{k} <= 0);
    cons = cons + (cons_LIN{k} == 0);
    cons = cons + (lambda{k} >= 0);
    cons = cons + (e{k} >= 0);
end
accumulated_ek  = 0;
for k = 1:N
    accumulated_ek = accumulated_ek + e{k};
end
obj = -tau;

solver_opt = sdpsettings('solver','mosek','verbose',verbose,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',tolerance);
solverDetails=optimize(cons,obj,solver_opt);
if reoptimize
    cons = cons+(tau>=-double(obj)-tol_reopt);
    solver_opt = sdpsettings('solver','mosek','verbose',verbose);
    solverDetails=optimize(cons,accumulated_ek,solver_opt);
end

%% Try to grasp what happens by plotting !
if pplot
    close all;
    ak = zeros(1, N+1);
    bk = zeros(1, N+1);
    ck = zeros(1, N+1);
    dk = zeros(1, N+1);
    ek = zeros(1, N+1);
    ekc= zeros(1, N+1);
    for i = 1:N+1
        ak(i)  = double(Q{i}(1,1));
        bk(i)  = double(Q{i}(2,2));
        ck(i)  = double(Q{i}(2,1));
        dk(i)  = double(d{i}(1,1));
    end
    for i = 2:N+1
        ek(i)  = double(e{i-1}(1));
        ekc(i) = sum(ek(1:i));
    end
    
    subplot(3,2,1);
    plot(1:N+1,ak,'-b'); hold on; title('ak (coefficient of ||xk-x*||^2)');
    subplot(3,2,2);
    plot(1:N+1,bk,'-b'); hold on; title('bk (coefficient of ||f''(xk)||^2)');
    subplot(3,2,3);
    plot(1:N+1,ck,'-b'); hold on; title('ck (coefficient of <f''(xk);xk-x*>)');
    subplot(3,2,4);
    plot(1:N+1,dk,'-b'); hold on; title('dk (coefficient of f(xk)-f(x*))');
    subplot(3,2,5);
    plot(1:N+1,ek,'-b'); hold on; title('ek (coefficient of sigma^2)');
    subplot(3,2,6);
    plot(1:N+1,ekc,'-b'); hold on; title('cumulated ek ');
    
    if ssave
        labels{1} = 'k'; labels{2} = 'ak'; labels{3} = 'bk'; labels{4} = 'ck'; labels{5} = 'dk'; labels{6} = 'ek'; labels{7} = 'ekc';
        data = [(1:N+1).' ak.'  bk.' ck.' dk.' ek.' ekc.'];
        saveData([folder nname],data,labels);
    end
end












