function Steady_state_PFR
%% Steady state PFR
% 
%   You'll learn:
%       +: How to solve steady state problems in DAE formulation
% 
%% The problem
% 
%   Chemical reactions:  A -> B   (R1)
%                       2A -> C   (R2)
% 
%   Differential Equations:
%   dF_i/dz = v_i1*r_1 + v_i2*r_2
%   dT/dz = 1/(cp*Ft)*[ K(Tc - T) - H1*r_1 - H2*r_2]
% 
%   Algebraic Constraints:
%   Ft = sum(F_i)
%   C_i = P/(R*T)*F_i/Ft
%   u*sum(C_i) = Ft
% 
%   Auxiliar equations
%   cp*Ft = sum(F_i*cp_i)
%   
%   ============================================================
%   Author: ataide@peq.coppe.ufrj.br 
%   homepage: github.com/asanet
%   Contact me for help/personal classes!
%   Date: 2018-07-05
%   Matlab version: R2018a

%% Problem setup
addpath('AuxFunctions')

% The parameters
H1 = -20000;        H2 = -60000;
cpa = 90;           cpb = 90;        cpc = 180;
K = 4000;           Tc = 373.15;     Tr = 300;
k1 = 1e-1;          E1 = 4000; 
k2 = 1e-6;          E2 = 9000;
R = 8.3145;         P = 101325;

% The integration interval (reactor length)
tspan = [0 1];

% The consistent initial conditions
Fa0 = 100; Fb0 = 0; Fc0 = 0; T0 = 423; Ft0 = 100;
Ca0 = P/R/T0; Cb0 = 0; Cc0 = 0; u0 = Ft0/(Ca0+Cb0+Cc0); 

y0 = [Fa0 Fb0 Fc0 T0 Ft0 u0 Ca0 Cb0 Cc0]';
yp0 = model(0,y0,zeros(9,1));

% The solution
[x,y] = ode15i(@model,tspan,y0,yp0);

% Plot the results
close all

Fas = y(:,1);   Fbs = y(:,2);   Fcs = y(:,3);  
Ts  = y(:,4);   Fts = y(:,5);   us  = y(:,6);
Cas = y(:,7);   Cbs = y(:,8);   Ccs = y(:,9);   


figured;
h = plot(x,[Fas Fbs Fcs]);
xlabel('Length')
ylabel('Molar flux')
set(h,'LineWidth',1.5);
legend({'F_a','F_b','F_c'})

figured;
h = plot(x,[Cas Cbs Ccs]);
xlabel('Length')
ylabel('Concentration')
set(h,'LineWidth',1.5);
legend({'C_a','C_b','C_c'})

figure;
subplot(311)
plot(x,Fts,'LineWidth',1.5)
xlabel('Length')
ylabel('Total flux')
subplot(312)
plot(x,Ts,'LineWidth',1.5)
xlabel('Length')
ylabel('Temperature')
subplot(313)
plot(x,us,'LineWidth',1.5)
xlabel('Length')
ylabel('Velocity')

    function res = model(~,y,yp)
        
        % Variable allocation
        Fa = y(1);  Fb = y(2);  Fc = y(3); T = y(4);
        
        Ft  = y(5);  u  = y(6);
        Ca = y(7);  Cb = y(8);  Cc = y(9);
        
        dFa = yp(1); dFb = yp(2); dFc = yp(3); dT = yp(4);
        
        % The reaction rates
        r1 = k1*exp(-E1*(1/T - 1/Tr))*Ca;
        r2 = k2*exp(-E2*(1/T - 1/Tr))*Ca^2;
        
        % The mean cp (not a constraint)
        cpFt = Fa*cpa + Fb*cpb + Fc*cpc;
        
        % The differential equations
        res(1,1) = -r1 -2*r2 - dFa;
        res(2,1) = r1 - dFb;
        res(3,1) = r2 - dFc;
        res(4,1) = 1/cpFt*( K*(Tc - T) - H1*r1 - H2*r2) - dT;
        
        % The algebraic constraints
        res(5,1) = Ft - Fa - Fb - Fc;
        res(6,1) = Ca - P/R/T*Fa/Ft;
        res(7,1) = Cb - P/R/T*Fb/Ft;
        res(8,1) = Cc - P/R/T*Fc/Ft;
        res(9,1) = u*(Ca+Cb+Cc) - Ft;
        
        
    end
   
   
end

