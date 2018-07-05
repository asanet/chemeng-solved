function Reaction_diffusion
%% Steady steate reaction-diffusion in a spherical particle
%              
%   You'll learn:
%       +: How to discretize by finite differences
%       +: How to solve non linear algebraic problems
%       +: How to use the sparsity algebra package
% 
%% The problem
% 
%   Differential Equation:
%   d^2y/dr^2 + 2/r*dy/dr - phi^2*y^m = 0
%              
%   Boundary Conditions:
%   r = 0 ... dy/dr = 0 
%   r = 1 ... y = ys
% 
%   ============================================================
%   Author: ataide@peq.coppe.ufrj.br
%   homepage: github.com/asanet
%   Date: 2018-07-05
%   Matlab version: R2018a
%   Contact me for help/personal classes!

%% Problem setup
addpath('AuxFunctions')

% Grid Points
n = 5000;

% Reaction order
m = 2;

% Thiele's module
phi = 20;

% Particle Radius
R = 1;

% Delta r
dr = R/(n-1);

% Domain of discretization
r = linspace(0,R,n)';

% Surface concentration
ys = 1;

% Sparsity pattern
B = ones(n,3);
Jp = spdiags(B,-1:1,n,n);

% Initial guess
y0 = ones(n,1);

%% Solution -> fsolve

% Solver Settings
opt = optimoptions(@fsolve,'TolFun',1e-10,'TolX',1e-10,'Display','iter-detailed',...
    'Algorithm','trust-region-reflective');

% Passing the Sparty Pattern (Compare solution time with and without it)
opt.JacobPattern = Jp;

% Solver call
tic
y_sol = fsolve(@model,y0,opt);
toc

% Compute the final residue
res = model(y_sol);

% Plot the data
close all;
figured;
h = plot(r,y_sol);
xlabel('Particle Radius')
ylabel('Concentration')
set(h,'LineWidth',1.5)

figured;
h = plot(r,res);
xlabel('Particle Radius')
ylabel('Residue')
set(h,'LineWidth',1.5)

%% Model
    function f = model(y)
        
        % Memory allocations
        f = zeros(n,1);
        
        % BC in r = 0
        f(1) = -3*y(1) + 4*y(2) -y(3) ;
        
        % Diff. Eq. (inner points)
        for i = 2:n-1
            f(i) = (y(i+1) - 2*y(i) + y(i-1))/dr^2 + 2./r(i).*( y(i+1) - y(i-1) )/2/dr - phi^2*y(i).^m;
        end
        
        % BC in r = 1
        f(n) = y(n) - ys; 
        
    end

end