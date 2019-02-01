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
nname   = 'SGD_Evaluation_AtAveragedPoint.dat';    % Name of the file

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
%   ak ||y_k-x*||^2 + bk ||f'(y_k)||^2 + 2 ck <f'(y_k); y_k-x*>
%       + dk (f(y_k)-f(x*)) + ak' || x_k-x* ||^2

% INTERMEDIARY POTENTIAL SETUPS (aims at reproducing results from
% Section 3.2.2) via option "relaxation".

% Options:
%   Set relax = 0: use the base potential stated below.
%   Set relax = 1: force ak = bk = ck = 0 and ak' = L/2

relax = 0;

% NOTE: this is for the sake of the example; much more states could be
% incorporated in the potentials (e.g., more dependence on x_k's)

% INITIAL AND FINAL POTENTIALS SETUP:

a0 = 0;
b0 = 0;
c0 = 0;
d0 = 0;
ap0= L/2; % a0'
tau = sdpvar(1); % variable to be maximized (see problem (2) in the paper)

aN  = 0;
bN  = 0; % bN is to be maximized
cN  = 0;
dN  = tau;
apN = 0; % aN'

if relax == 1
    a0 = 0; aN = 0;
    b0 = 0; bN = 0;
    c0 = 0; cN = 0;
    ap0 = L/2; % a0'
    apN = L/2; % aN'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              %
% SETTING UP THE NOTATIONS FOR THE LINEAR MATRIX INEQUALITIES  %
%                   (end of editable zone)                     %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recall that x* (optimum) is set to x* = 0 without loss of generality,
% and so does f(x*) = 0.

% P = [ xk  y_k  f'(y_{k}) G(y_{k+1};1) ... G(y_{k+1};n)]
% F = [           f(y_k)         f(y_{k+1})  ]
% NOTE: Symmetry arguments allow simplifying both F and G=P.'*P.

dimG  = 3+n;  % dimensions of the Gram matrix P.'*P
dimF  = 2;    % dimensions of F
nbPts = 3;    % number of points to be incorporated in the discrete
% version of the function f:
% x*, y_{k}, y_{k+1}

% Notations for the LMI
xk  = zeros(1, dimG); xk(1)  = 1; % this is x_{k}
yk  = zeros(1, dimG); yk(2)  = 1; % this is y_{k}
xs  = zeros(1, dimG);             % this is x* = 0

gyk = zeros(1, dimG); gyk(:,3) = 1;          % this is f'(y_{k})
Gyk1= zeros(n, dimG);Gyk1(:,4:n+3) = eye(n); % Gyk1(i,:) is G(y_{k+1};i)
gyk1= sum(Gyk1,1)/n;                         % gyk1 is f'(y_{k+1})
gxs = zeros(1, dimG);                        % this is f'(x*) = 0

fyk = zeros(1, dimF); fyk(1) = 1;           % this is f(y_{k})
fyk1= zeros(1, dimF);fyk1(2) = 1;           % fyk1 is f(y_{k+1})
fxs = zeros(1, dimF);                       % this is f(x*) = 0

yk1 = @(k)(1/(k)*xk + (k-1)/(k) * yk);
xk1 = @(k)(ones(n,1)*xk - delta(k) * Gyk1);  % xk1(i,:) is x_{k+1}^{(i)}

% Matrix encoding interpolation condition for smooth strongly convex
% functions
M = 1/2/(L-m) *[   -L*m,  L*m,   m, -L;
                    L*m, -L*m,  -m,  L;
                      m,   -m,  -1,  1;
                     -L,   L,   1,  -1];
% Matrix encoding the variance condition
Avar = zeros(dimG);
for i = 1:n
    Avar = Avar + (Gyk1(i,:)-gyk1).'*(Gyk1(i,:)-gyk1)/n;
end
% Notations for potentials
Q0 = [ a0 c0; c0 b0]; q0 = d0; 
QN = [ aN cN; cN bN]; qN = dN;

statesKG  = [(yk-xs).' gyk.'];  statesKF = fyk;  statesK1F = fyk1;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              %
%          SETTING UP THE LINEAR MATRIX INEQUALITIES           %
%     (directly embedded within the larger problem (2))        %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q{1}   = Q0;
d{1}   = q0;
ap{1}  = ap0;
Q{N+1} = QN;
d{N+1} = qN;
ap{N+1}= apN;

cons   = (tau>=0);
for k = 1 : N % iteration counter (one LMI per iteration)
    if k < N
        if ~relax % no constraint on Qk and dk
            d{k+1} = sdpvar(1);
            ap{k+1}= sdpvar(1);
            Q{k+1} = sdpvar(2);
        elseif relax == 1 % force ak = bk = ck = 0; akp = L/2;
            d{k+1} = sdpvar(1);
            ap{k+1}= L/2;
            Q{k+1} = zeros(2);
        end
    end
    
    % POINTS IN THE DISCRETE REPRESENTATION OF FUNCTION f(.)
    clear X G F;
    XK1 = xk1(k);      % xk1(i,:) is x_{k+1}^{(i)}
    X = {  xs,   yk,  yk1(k)}; % coordinates
    G = { gxs,  gyk,  gyk1}; % gradients
    F = { fxs,  fyk,  fyk1}; % function values
    
    lambda{k}    = sdpvar(nbPts,nbPts,'full');
    e{k}         = sdpvar(1);
    cons_SDP{k}  = - statesKG * Q{k} * statesKG.'...
        - ap{k} * (xk-xs).'*(xk-xs) - e{k} * Avar;
    cons_LIN{k}  = - d{k} * statesKF;
    
    
    statesK1G = [(yk1(k)-xs).' gyk1.'];
    cons_SDP{k} = cons_SDP{k} + statesK1G * Q{k+1} * statesK1G.';
    for i = 1:n % averaging the states over i of x_{k+1}^{(i)}
        cons_SDP{k}= cons_SDP{k} + ap{k+1} * (XK1(i,:)-xs).'*(XK1(i,:)-xs)/n;
    end
    cons_LIN{k}= cons_LIN{k} + statesK1F * d{k+1};
    
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
    apk = zeros(1, N+1);
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
        apk(i) = double(ap{i});
    end
    for i = 2:N+1
        ek(i)  = double(e{i-1}(1));
        ekc(i) = sum(ek(1:i));
    end
    
    subplot(4,2,1);
    plot(1:N+1,ak,'-b'); hold on; title('ak (coefficient of ||yk-x*||^2)');
    subplot(4,2,2);
    plot(1:N+1,bk,'-b'); hold on; title('bk (coefficient of ||f''(yk)||^2)');
    subplot(4,2,3);
    plot(1:N+1,ck,'-b'); hold on; title('ck (coefficient of <f''(yk);yk-x*>)');
    subplot(4,2,4);
    plot(1:N+1,dk,'-b'); hold on; title('dk (coefficient of f(yk)-f(x*))');
    subplot(4,2,5);
    plot(1:N+1,ek,'-b'); hold on; title('ek (coefficient of sigma^2)');
    subplot(4,2,6);
    plot(1:N+1,ekc,'-b'); hold on; title('cumulated ek ');
    subplot(4,2,7);
    plot(1:N+1,apk,'-b'); hold on; title('ak'' (coefficient of ||yk-x*||^2)');
    
    if ssave
        labels{1} = 'k'; labels{2} = 'ak'; labels{3} = 'bk'; labels{4} = 'ck'; labels{5} = 'dk'; labels{6} = 'ek'; labels{7} = 'ekc'; labels{8} = 'apk';
        data = [(1:N+1).' ak.'  bk.' ck.' dk.' ek.' ekc.' apk.'];
        saveData([folder nname],data,labels);
    end
end












