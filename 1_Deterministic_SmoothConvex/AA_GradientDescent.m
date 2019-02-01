clear all;
clc;

% SOLVER OPTIONS
verbose     = 1;    
tolerance   = 1e-8;

% OUTPUT OPTIONS
ssave   = 1; % Save the results ? 
pplot   = 1; % Plot the results ?
folder  = 'SaveData/';      % If results saved, name of the saving folder
nname   = 'GradientD.dat';  % Name of the file

% MINIMIZATION PROBLEM SETUP
L = 1;      % Smoothness constant
m = 0;      % Strong convexity constant
N = 100;    % Number of iterations

% ALGORITHM SETUP: gradient descent with step-size policy delta_k.
delta = @(k)(1/L);   % Step-size (possibly varying function of k)

% INTERMEDIARY POTENTIAL SETUPS (aims at reproducing results from
% Section 3.2.1) via option "relaxation".

% Options: 
%   Set relax = 0: use the base potential stated in Section 3.2.1.
%   Set relax = 1: force dk = (2*k+1)*L
%   Set relax = 2: force dk = (2*k+1)*L, ck = 0 and ak = L^2.
%   Set relax = 3: force dk = 0.

relax = 0;

% INITIAL AND FINAL POTENTIALS SETUP:
% Potential has the form: 
%   ak ||x_k-x*||^2 + bk ||f'(x_k)||^2 + 2 ck <f'(x_k); x_k-x*> 
%       + dk (f(x_k)-f(x*)).

a0  = L^2;
b0 = 0;
c0 = 0;
d0 = 0;

tau = sdpvar(1); % variable to be maximized (see problem (2) in the paper)

aN  = 0; 
bN  = tau; % bN is to be maximized
cN  = 0;
dN  = 0;

if relax == 1
    d0 = L;
    dN = (2*N+1)*L;
elseif relax == 2
    d0 = L;
    dN = (2*N+1)*L;
    aN = L^2;
    cN = 0;
elseif relax == 3
    d0 = 0;
    dN = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              %
% SETTING UP THE NOTATIONS FOR THE LINEAR MATRIX INEQUALITIES  %
%                   (end of editable zone)                     %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recall that x* (optimum) is set to x* = 0 without loss of generality,
% and so does f(x*) = 0.

% P = [ xk  f'(x_k)  f'(x_{k+1})    ]
% F = [     f(x_k)   f(x_{k+1})     ]

dimG  = 3; % dimensions of the Gram matrix P.'*P
dimF  = 2; % dimensions of F
nbPts = 3; % number of points to be incorporated in the discrete 
           % version of the function f: x*, x_{k} and x_{k+1}

% Notations for the LMI
xk  = zeros(1, dimG); xk(1)  = 1; % this is x_{k}
xs  = zeros(1, dimG);             % this is x* = 0

gxk = zeros(1, dimG); gxk(2) = 1; % this is f'(x_{k})
gxk1= zeros(1, dimG);gxk1(3) = 1; % this is f'(x_{k+1})
gxs = zeros(1, dimG);             % this is f'(x*) = 0

fxk = zeros(1, dimF); fxk(1) = 1; % this is f(x_{k})
fxk1= zeros(1, dimF);fxk1(2) = 1; % this is f(x_{k+1}) 
fxs = zeros(1, dimF);             % this is f(x*) = 0

xk1 = @(k)(xk - delta(k) * gxk);  % this is x_{k+1} 

% Matrix encoding interpolation condition for smooth strongly convex
% functions
M = 1/2/(L-m) *[  -L*m,  L*m,   m, -L; 
                   L*m, -L*m,  -m,  L;
                     m,   -m,  -1,  1;
                    -L,   L,   1,  -1];

% Notations for potentials
Q0 = [ a0 c0; c0 b0]; q0 = d0;
QN = [ aN cN; cN bN]; qN = dN;

statesKG  = [(xk-xs).' gxk.'];          statesKF = fxk;
statesK1G = @(k)[(xk1(k)-xs).' gxk1.']; statesK1F = fxk1;

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

cons = (tau >= 0);
for k = 1 : N % iteration counter (one LMI per iteration)
    if k < N
        if ~relax % no constraint on Qk and dk
            d{k+1} = sdpvar(1);
            Q{k+1} = sdpvar(2);
        elseif relax == 1 % force dk = (2*k+1)*L and no constraint on Qk
            d{k+1} = (2*k+1)*L;
            Q{k+1} = sdpvar(2);
        elseif relax == 2 % force dk = (2*k+1)*L, ck = 0 and ak = L^2
            d{k+1} = (2*k+1)*L;
            Q{k+1} = [L^2 0; 0 sdpvar(1)];
        elseif relax == 3 % force dk = 0 and no constraint on Qk
            d{k+1} = 0;
            Q{k+1} = sdpvar(2);
        end
    end
    
    lambda{k}    = sdpvar(nbPts,nbPts,'full'); 
    cons_SDP{k}  = - statesKG * Q{k} * statesKG.' ...
        + statesK1G(k) * Q{k+1} * statesK1G(k).';
    
    cons_LIN{k} = - d{k} * statesKF + d{k+1} * statesK1F;
    
    % POINTS IN THE DISCRETE REPRESENTATION OF FUNCTION f(.)
    X = {  xs,    xk,   xk1(k)}; % coordinates
    G = { gxs,   gxk,   gxk1}; % gradients
    F = { fxs,   fxk,   fxk1}; % function values

    % ADD ALL INTERPOLATION CONDITIONS TO THE LMI
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
end
obj = tau;

solver_opt = sdpsettings('solver','mosek','verbose',verbose,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',tolerance);
solverDetails=optimize(cons,-obj,solver_opt);


%% Try to grasp what happens by plotting !
if pplot
    close all;
    ak = zeros(1, N+1);
    bk = zeros(1, N+1);
    ck = zeros(1, N+1);
    dk = zeros(1, N+1);
    
    for i = 1:N+1
        ak(i)  = double(Q{i}(1,1)).';
        bk(i)  = double(Q{i}(2,2)).';
        ck(i)  = double(Q{i}(2,1)).';
        dk(i)  = double(d{i}(1,1)).';
    end
    
    subplot(2,2,1);
    plot(1:N+1,ak,'-b'); hold on; title('ak (coefficient of ||xk-x*||^2)');
    subplot(2,2,2);
    plot(1:N+1,bk,'-b'); hold on; title('bk (coefficient of ||f''(xk)||^2)');
    subplot(2,2,3);
    plot(1:N+1,ck,'-b'); hold on; title('ck (coefficient of <f''(xk);xk-x*>)');
    subplot(2,2,4);
    plot(1:N+1,dk,'-b'); hold on; title('dk (coefficient of f(xk)-f(x*))');
    
    if ssave
        labels{1} = 'k'; labels{2} = 'ak'; labels{3} = 'bk'; labels{4} = 'ck'; labels{5} = 'dk';
        data = [(1:N+1).' ak.'  bk.' ck.' dk.' ];
        saveData([folder nname],data,labels);
    end
end












