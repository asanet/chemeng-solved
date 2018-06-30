function Extractors
%% Steady state train of extractors
%
%   You'll learn:
%       +: How to generalize staged processes
%       +: How to solve linear systems
% 
%% The problem
% 
%   Physical model:
% 
%                Stg (i-1)       Stg (i)        Stg (i+1)
%             -----------------------------------------------
%   X(i-2) -->|             |-->X(i-1)      |-->X(i)        |-->X(i+1)  (R)
%             |             |               |               |
%             |             |               |               |
%             |             |               |               |
%   Y(i-1) <--|      Y(i)<--|      Y(i+1)<--|               |<--Y(i+2)  (E)
%             -----------------------------------------------
%
%   Input:
%   N:  Number of Stages
%   R:  Raffinate flow rate 
%   E:  Extract flow rate  
%   Xf: Molar fration feeding stage 1 at raffinate 
%   Yf: Molar fration feeding stage N at extract 
%   Ki: Equilibrium constants
%
%   Mass balance:
%   Stg 1: E*Y(2) + R*Xf = R*X(1) + E*Y(1)
%   Stg i: E*Y(i+1) + R*X(i-1) = R*X(i) + E*Y(i)
%   Stg N: E*Yf + R*X(N-1) = R*X(N) + E*Y(N)
%
%   Equilibrium constraints
%   Y(i) = K(i)*X(i)
% 
%   Problem: Ax=b
% 
%   ============================================================
%   Author: ataide@peq.coppe.ufrj.br
%   Contact me for help/personal classes!

%% Problem setup
addpath('AuxFunctions')

% Numer of stages
N = 10;
K = 0.25*ones(N,1);
R = 1;
E = 5;
Xf = 0.8;
Yf = 0;

% Building the matrix A
% Prototype
%        X1   X2   X3   Y1   Y2   Y3     
%
% A =    -R    0    0 | -E    E    0    EQ1
%         R   -R    0 |  0   -E    E    EQ2
%         0    R   -R |  0    0   -E    EQ3
%        -------------+-------------
%       -K1    0    0 |  1    0    0    EQ4
%         0  -K2    0 |  0    1    0    EQ5
%         0    0  -K3 |  0    0    1    EQ6  
 
% Generalized
MR = -R*diag(ones(N,1)) + R*diag(ones(N-1,1),-1);
ME = -E*diag(ones(N,1)) + E*diag(ones(N-1,1),1);
MK = -diag(K);
MI = eye(N);

A = [MR ME
     MK MI];

% Vector b
b = zeros(2*N,1);
b(1) = -R*Xf;
b(N) = -E*Yf;

% Solve the linear system
sol = linsolve(A,b);

% Plot the data
close all

Stages = 1:N;
Xi = sol(1:N);
Yi = sol(N+1:2*N);

figured;
h = plot(Stages,Xi,Stages,Yi);
set(h(1),'Marker','o','LineStyle',':','LineWidth',1.5,'MarkerFaceColor',get(h(1),'Color'));
set(h(2),'Marker','o','LineStyle',':','LineWidth',1.5,'MarkerFaceColor',get(h(2),'Color'));
xlabel('Stage')
ylabel('Molar fraction')
axis([1 N 0 1])
legend({'Raffinate','Extract'})












