function Furnace_Wall
%% Transient heat conduction on a composite furnace wall
% 
%   You'll learn:
%       +: How to discretize a domain by orthogonal collocation
%       +: How to handle problems with multiple domains
%       +: How to solve transient problems in DAE formulation
% 
%% The problem
% 
%   Physical model
% 
%                      \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                     |----------------------------------|  
%       h1, Tinf1     |      Body A     |      Body B    |  h2, Tinf2 
%                     |       TA0       |       TB0      |  
%                     |_________________|________________|
%                       \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                  -LA|---------------- 0 -------------->|LB
% 
%   Differential Equations:
%   dTA/dt = (k/rho/cp)_A d^2TA/dx^2
%   dTB/dt = (k/rho/cp)_B d^2TB/dx^2
% 
%   Boundary Conditions:
%   x = -LA ... -kA*dTA/dx = h1*(Tinf1 - TA)
%   x =   0 ... TA = TB
%               kA*dTA/dx = kB*dTB/dx 
%   x =  LB ... -kB*dTB/dx = h2*(TB - Tinf2)
% 
%   Initial Condition:
%   t = 0 ... TA = TA0(x)
%             TB = TB0(x)
% 
%   ============================================================
%   Author: ataide@peq.coppe.ufrj.br
%   homepage: github.com/asanet
%   Contact me for help/personal classes!
%   Date: 2018-07-05
%   Matlab version: R2018a

%% Problem setup
addpath('AuxFunctions')

% Parameters
kA = 50*3600;  kB = 398*3600;
cpA = 900;      cpB = 386;
rhoA = 2697;    rhoB = 8920;
LA = 1;         LB = 1.5;
h1 = 1e2*3600;  h2 = 1e5*3600;
Tinf1 = 298.15; Tinf2 = 1000;

% Inner points
N = 5;

% Total points
n = N + 2;

% Grid and matrices
[M1,M2,xgrid] = collocation(N,'2*x-1',[0 1]);
xB = LB*xgrid;
xA = sort(-LA*xgrid);

% Initial condition
y0 = 298.15*ones(2*n,1);
yp0 = model(0,y0,zeros(2*n,1));

tf = 20;
frames = 900 + 100;
tspan = [ linspace(0, tf/10, 900) linspace( tf/10+0.1, tf, 100)];
[t,y] = ode15i(@model,tspan,y0,yp0);

% Solution at the collocation points
TAs = y(:,1:n); TBs = y(:,n+1:2*n);

% Refined grid
ref = 100;
xAc = linspace(-LA,0,ref);
xBc = linspace(0,LB,ref);

% Lagrange polynomial
TAc = zeros(frames,ref);
TBc = zeros(frames,ref);

labels = cell(frames,1);
for i = 1:frames
    TAc(i,:) = lagrange(xA,TAs(i,:),xAc);
    TBc(i,:) = lagrange(xB,TBs(i,:),xBc);
    labels{i} = sprintf('Temperature profile in t = %2.4f h',t(i));
end

% Plot data
close all
figured(1);
h = plot(t,TAs(:,1:n-1),'b',t,TAs(:,n),'g',t,TBs(:,2:end),'r');
set(h,'LineWidth',1.5)
set(h([2:n-1 n+2:2*n-1]),'HandleVisibility','off')
xlabel('Time (h)')
ylabel('Temperature (K)') 
legend({'Body A','Interface','Body B'},'Location','SouthEast')

pause(1)
figured(2);
axis([-LA LB .95*Tinf1 1.05*Tinf2])
xlabel('Length (m)')
ylabel('Temperature (K)')
opts = { {'LineWidth',1.5,'Color','b'}; 
         {'Marker','*','LineStyle','none','Color','b'};
         {'LineWidth',1.5,'Color','r'}; 
         {'Marker','*','LineStyle','none','Color','r'} };
    
legopt = {'A: Interp','A: Col point','B: Interp','B: Col point','location','NorthOutside','Orientation','Horizontal'};
sldplot({xAc xA xBc xB},{TAc TAs TBc TBs},opts,labels,legopt)

    % The model (discretization by orthogonal collocation)
    function res = model(~,y,yp)
        
        % Variable allocation
        TA = y(1:n);    TB = y(n+1:2*n);
        dTA = yp (1:n); dTB = yp(n+1:2*n);
        res = zeros(2*n,1);
        
        % Boundary Conditions in x = -LA
        res(1) = -kA*M1(1,:)*TA/LA - h1*(Tinf1 - TA(1));

        % Differential Equation for TA (inner points)
        res(2:n-1) = kA/rhoA/cpA*(M2(2:n-1,:)*TA)/LA^2 - dTA(2:n-1);
        
        % Boundary Conditions in x = 0
        res(n) = TA(n) - TB(1);
        res(n+1) = kA*M1(n,:)*TA/LA - kB*M1(1,:)*TB/LB;
        
        % Differential Equation for TB (inner points)
        res(n+2:2*n-1) = kB/rhoB/cpB*(M2(2:n-1,:)*TB)/LB^2 - dTB(2:n-1);
        
        % Boundary Conditions in x = LB
        res(2*n) = -kB*M1(n,:)*TB/LB - h2*(TB(n) - Tinf2);
        
    end
end
