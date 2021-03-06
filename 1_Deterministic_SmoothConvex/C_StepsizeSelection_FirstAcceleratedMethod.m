clear all;
clc;
tic;
% SOLVER OPTIONS
verbose     = 1;
tolerance   = 1e-8;

% OUTPUT OPTIONS
ssave   = 1; % Save the results ?
pplot   = 1; % Plot the results ?
folder = 'SaveData/';           % If results saved, name of the saving folder
nname  = 'FastGradientD1.dat';  % Name of the file

% MINIMIZATION PROBLEM SETUP
L = 1;      % Smoothness constant
m = 0;      % Strong convexity constant
N = 100;    % Number of iterations

% INTERMEDIARY POTENTIAL SETUPS (aims at reproducing results from Appendix
% C.3; Figure 3) via option "relaxation".
% Options:
%   Set relax = 0: use the base potential stated in Section 3.2.1.
%   Set relax = 1: force ak = L/2; Qk = 0
%   Set relax = 2: force ak = L/2, Qk(1,1)=Qk(1,2)=0, Qk(2,2)=dk/2/L

relax = 0;

% INITIAL AND FINAL POTENTIALS SETUP:
% Potential has the form:
%   Q(1,1) ||x_k-x*||^2 + Q(2,2) ||F'(x_k)||^2 + 2 Q(1,2) <F'(x_k); x_k-x*>
%       + ak ||z_k-x*||^2 + dk (F(x_k)-f(x*)).

Q0 = zeros(2);
a0 = L/2;
d0 = 0;

tau = sdpvar(1); % variable to be maximized (see problem (12) in the paper)
QN  = zeros(2);
aN  = 0;
dN  = tau;

if relax == 1
    a0 = L/2; aN = L/2;
    Q0 = zeros(2); QN = zeros(2);
elseif relax == 2
    a0 = L/2; aN = L/2;
    Q0 = [0 0; 0 d0/2/L]; QN = [0 0; 0 dN/2/L];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              %
% SETTING UP THE NOTATIONS FOR THE LINEAR MATRIX INEQUALITIES  %
%                   (end of editable zone)                     %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimization of f(x)  with f smooth (possibly strongly) convex.
% Recall that x* (optimum) is set to x* = 0 without loss of generality,
% and so does f(x*) = 0.

% P = [ z_k x_k x_{k+1} y_{k+1} f'(x_k)  f'(x_{k+1}) f'(y_{k+1}) ]
% F = [                         f(x_k)   f(x_{k+1})  f(y_{k+1})]

dimG  = 7; % dimensions of the Gram matrix P.'*P
dimF  = 3; % dimensions of F
nbPts = 4; % number of points to be incorporated in the discrete
% version of the function f: x*, x_k, x_{k+1}, y_{k+1}

% Notations for the LMI
zk  = zeros(1, dimG); zk(1)  = 1; % this is z_k
xk  = zeros(1, dimG); xk(2)  = 1; % this is x_k
xk1 = zeros(1, dimG); xk1(3) = 1; % this is x_{k+1}
yk1 = zeros(1, dimG); yk1(4) = 1; % this is y_{k+1}
xs  = zeros(1, dimG);             % this is x*

gxk = zeros(1, dimG); gxk(5) = 1; % this is f'(x_k)
gxk1= zeros(1, dimG);gxk1(6) = 1; % this is f'(x_{k+1})
gyk1= zeros(1, dimG);gyk1(7) = 1; % this is f'(y_{k+1})
gxs = zeros(1, dimG);             % this is f'(x*) = 0

fxk = zeros(1, dimF); fxk(1) = 1; % this is f(x_k)
fxk1= zeros(1, dimF);fxk1(2) = 1; % this is f(x_{k+1})
fyk1= zeros(1, dimF);fyk1(3) = 1; % this is f(y_{k+1})
fxs = zeros(1, dimF);             % this is f(x*)

% POINTS IN THE DISCRETE REPRESENTATION OF FUNCTION f(.)
X = {  xs,  xk,  xk1,  yk1}; % coordinates (notations w in the paper)
G = { gxs, gxk, gxk1, gyk1}; % gradients
F = { fxs, fxk, fxk1, fyk1}; % function values

% line-search for y_{k+1}
A{1} = gyk1.'*(yk1-xk);   A{1} = 1/2*(A{1}+A{1}.'); % <f'(y_{k+1}); y_{k+1}-x_k> == 0
A{2} = gyk1.'*(zk-xk);    A{2} = 1/2*(A{2}+A{2}.'); % <f'(y_{k+1}); z_k    -x_k> == 0
% line-search for x_{k+1}
A{3} = gxk1.'*(yk1-xk1);  A{3} = 1/2*(A{3}+A{3}.'); % <f'(x_{k+1}); y_{k+1}-x_{k+1}> == 0
A{4} = gxk1.'*gyk1;       A{4} = 1/2*(A{4}+A{4}.'); % <f'(x_{k+1}); f'(y_{k+1})> == 0

% Matrix encoding interpolation condition for smooth strongly convex functions
M = 1/2/(L-m) *[  -L*m,  L*m,   m, -L;
    L*m, -L*m,  -m,  L;
    m,   -m,  -1,  1;
    -L,   L,   1,  -1];

% Notations for potentials
statesKG  = [xk-xs; gxk];     statesKF  = fxk-fxs;
statesK1G = [xk1-xs; gxk1];   statesK1F = fxk1-fxs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              %
%          SETTING UP THE LINEAR MATRIX INEQUALITIES           %
%     (directly embedded within the larger problem (12))       %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q{1}   = Q0;
d{1}   = d0;
a{1}   = a0;
Q{N+1} = QN;
d{N+1} = dN;
a{N+1} = aN;

cons = (tau >= 0);
for k = 1 : N % iteration counter (one LMI per iteration)
    if k < N
        if ~relax % no constraint on Qk, dk and ak
            Q{k+1} = sdpvar(2);
            d{k+1} = sdpvar(1);
            a{k+1} = sdpvar(1);
        elseif relax == 1 % force ak = L/2; Qk = 0
            d{k+1} = sdpvar(1);
            a{k+1} = L/2;
            Q{k+1} = zeros(2);
        elseif relax == 2 % force ak = L/2, Qk(1,1)=Qk(1,2)=0, Qk(2,2)=dk/2/L
            d{k+1} = sdpvar(1);
            a{k+1} = L/2;
            Q{k+1} = [0 0; 0 d{k+1}/2/L];
        end
    end
    
    lambda{k}   = sdpvar(nbPts,nbPts,'full');
    mu{k}       = sdpvar(4,1); % multipliers for line-search conditions
    
    
    S{k} = sdpvar(3); % this is Sk (see Section C.5)
    cons = cons + (S{k} >= 0);
    cons = cons + (S{k}(1,1) == a{k+1});
    
    % Form the SDP constraint
    cons_SDP{k} = - a{k} * (zk-xs).'*(zk-xs) ...
        - statesKG.' * Q{k} * statesKG ...
        + [yk1; zk-yk1; gyk1].' * S{k} * [yk1; zk-yk1; gyk1] ...
        + statesK1G.' * Q{k+1} * statesK1G;
    for i = 1:4
        cons_SDP{k} = cons_SDP{k} + mu{k}(i) * A{i};
    end
    cons_LIN{k} = - d{k} * statesKF + d{k+1} * statesK1F;
    
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

toc

%% Try to grasp what happens by plotting !
if pplot
    close all;
    Q11 = zeros(1, N+1);
    Q12 = zeros(1, N+1);
    Q22 = zeros(1, N+1);
    ak  = zeros(1, N+1);
    dk  = zeros(1, N+1);
    alphak = zeros(1, N+1);
    tauk   = zeros(1, N+1);
    deltak = zeros(1, N+1);
    gammak = zeros(1, N+1);
    
    tolp = 1e-5;
    for i = 1:N+1
        Q11(i) = double(Q{i}(1,1));
        Q12(i) = double(Q{i}(1,2));
        Q22(i) = double(Q{i}(2,2));
        ak(i)  = double(a{i});
        dk(i)  = double(d{i}).';
    end
    for i = 1:N
        tauk(i)   = -double(mu{i}(2)/mu{i}(1));
        alphak(i) =  double(mu{i}(4)/mu{i}(3));
        gammak(i) = -double(S{i}(1,3)/S{i}(1,1));
        deltak(i) =  double(S{i}(1,2)/S{i}(1,1));
    end
    
    subplot(3,3,1);
    plot(1:N+1,Q11,'-b'); hold on; title('Q11 (coefficient of ||xk-x*||^2)');
    subplot(3,3,2);
    plot(1:N+1,Q12,'-b'); hold on; title('Q12 (coefficient of <f''(xk); xk-x*>)');
    subplot(3,3,3);
    plot(1:N+1,Q22,'-b');  hold on; title('Q22 (coefficient of ||f''(xk)||^2)');
    subplot(3,3,4);
    plot(1:N+1,ak,'-b');  hold on; title('ak (coefficient of ||zk-x*||^2)');
    subplot(3,3,5);
    plot(1:N+1,dk,'-b');  hold on; title('dk (coefficient of f(xk)-f(x*))');
    subplot(3,3,6);
    plot(1:N,tauk(1:end-1),'-b'); hold on; title('tauk');
    subplot(3,3,7);
    plot(1:N,alphak(1:end-1),'-b'); hold on; title('alpha');
    subplot(3,3,8);
    plot(1:N,gammak(1:end-1),'-b'); hold on; title('gammak');
    subplot(3,3,9);
    plot(1:N,deltak(1:end-1),'-b'); hold on; title('deltak');
    
    if ssave
        labels{1} = 'k'; labels{2} = 'Q11'; labels{3} = 'Q12'; labels{4} = 'Q22'; labels{5} = 'ak'; labels{6} = 'dk'; labels{7} = 'alphak'; labels{8} = 'tauk'; labels{9} = 'deltak'; labels{10} = 'gammak';
        data = [(1:N+1).' Q11.'  Q12.' Q22.' ak.' dk.' alphak.' tauk.' deltak.' gammak.'];
        saveData([folder nname],data,labels);
    end
end











