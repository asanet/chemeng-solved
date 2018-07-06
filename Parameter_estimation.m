function Parameter_estimation
%% Parameter estimation of a dynamic model
%   
%   You'll learn:
%       +: How to solve optimization problems
%       +: How to estimate parameters of a non linear model with numerical
%          solution
%   
%% The problem
%   
%   Given a data set (xe,ye), find the parameters of a dynamical
%   non linear model. By least squares the problem is:
%   
%   minimize the sum of (ye_i - ycalc_i)^2
%               
%   About the process:
%   CSTR with van de vusse reaction system
%    A -> B
%    B -> C
%   2A -> D
%   
%   Estimate the Arrhenius pre-exponential factor and activation energies
% 
%   ============================================================
%   Author: ataide@peq.coppe.ufrj.br
%   homepage: github.com/asanet
%   Date: 2018-07-05
%   Matlab version: R2018a
%   Contact me for help/personal classes!

%% Problem setup
addpath('AuxFunctions')

% Read the data set from the dataset.xls file
dataset = xlsread('dataset.xls');

% Separate the data set into vectors
te = dataset(:,1);  % The time (independent variable)
Cae = dataset(:,2);  % The temperature
Cbe = dataset(:,3);  % The temperature
Te = dataset(:,4);  % The temperature
                    
% The initial guess for the parameters
k10 = 1e8;  k20 = 1e8;  k30 = 1e3;
E10 = 1e5;  E20 = 1e5;  E30 = 1e5; 

par0 = [k10 k20 k30 E10 E20 E30]';

% The known model parameters
H1 = 4.2e3;         H2 = -11e3;     H3 = -41.85e3;
rho = 934.2;        cp = 3.01e3;    V = 1e-2;      
tau = 80;           Tf = 403.15;    Caf = 1000;
UA = 0.215*1120;    R = 8.3145;     Tk = 402.1;

% Configure the optimization solver
op = optimset('Display','iter','MaxIter',700,'MaxFunEvals',1e5,'TolFun',1e-8,'TolX',1e-8);
odeopt = odeset('Abstol',1e-6,'Reltol',1e-4);

% The initial condition is part of the dataset!!!
y0 = [Cae(1) Cbe(1) Te(1)]';

% Call the optimization solver
parEst = fminsearch(@fobj,par0,op);

% Calculate ycalc with the estimated parameters for comparison
tspan = linspace(0,te(end),100)';
[tc,yc] = ode15s(@model,tspan,y0,odeopt,parEst);
[~,yc2] = ode15s(@model,te,y0,odeopt,parEst);

% Plot data
close all

Cac = yc(:,1);      Cac2 = yc2(:,1);         
Cbc = yc(:,2);      Cbc2 = yc2(:,2);
Tc  = yc(:,3);      Tc2  = yc2(:,3);

% measurement errors
errCa = 20*ones(size(Cae));
errCb = 20*ones(size(Cbe));
errT  = 5*ones(size(Te));

colors = get(0, 'DefaultAxesColorOrder');

figured;
xlabel('Time (s)')
ylabel('Concentration (mol \cdot m^{-3})')
plot(tc,Cac,tc,Cbc,'LineWidth',1.5);
hold on
errorbar(te,Cae,errCa,'-s','MarkerSize',8,'LineStyle','none', ...
            'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),'Color',colors(1,:))

errorbar(te,Cbe,errCb,'-d','MarkerSize',8,'LineStyle','none', ...
            'MarkerEdgeColor',colors(2,:),'MarkerFaceColor',colors(2,:),'Color',colors(2,:))

legend({'Ca_{calc}','Cb_{calc}','Ca_{exp}','Cb_{exp}'},'location','southeast')
hold off

figured;
xlabel('Time (s)')
ylabel('Temperature (K)')
plot(tc,Tc,'LineWidth',1.5,'Color',colors(3,:));
hold on
errorbar(te,Te,errT,'-o','MarkerSize',6,'LineStyle','none', ...
            'MarkerEdgeColor',colors(3,:),'MarkerFaceColor',colors(3,:),'Color',colors(3,:))

legend({'Calculated','Experimental'},'location','southeast')
hold off

figured;
subplot(311)
plot(Cac2,Cac2,'lineWidth',1.5)
hold on
errorbar(Cac2,Cae,errCa,'-s','MarkerSize',8,'LineStyle','none', 'handlevisibility','off', ...
            'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),'Color',colors(1,:))
legend('Ca','location','southeast')
set(gca,'Fontsize',16,'ygrid','on','xlim',[0 max(Cae)],'ylim',[0 max(Cae)])
hold off

subplot(312)
plot(Cbc2,Cbc2,'lineWidth',1.5,'Color',colors(2,:))
ylabel('Experimental values')
hold on
errorbar(Cbc2,Cbe,errCb,'-d','MarkerSize',8,'LineStyle','none', 'handlevisibility','off', ...
            'MarkerEdgeColor',colors(2,:),'MarkerFaceColor',colors(2,:),'Color',colors(2,:))
legend('Cb','location','southeast')
set(gca,'Fontsize',16,'ygrid','on','xlim',[0 max(Cbe)],'ylim',[0 max(Cbe)])
hold off

subplot(313)
plot(Tc2,Tc2,'lineWidth',1.5,'Color',colors(3,:))
xlabel('Predicted values')
hold on
errorbar(Tc2,Te,errT,'-o','MarkerSize',6,'LineStyle','none', 'handlevisibility','off', ...
            'MarkerEdgeColor',colors(3,:),'MarkerFaceColor',colors(3,:),'Color',colors(3,:))
legend('Temperature','location','southeast')
set(gca,'Fontsize',16,'ygrid','on','xlim',[min(Te) max(Te)],'ylim',[min(Te) max(Te)])
hold off


    function f = fobj(par)
        % ycalc must be evaluated at the same experimental time points ti 
        [~,ycalc] = ode15s(@model,te,y0,odeopt,par);
        f = (Cae - ycalc(:,1)).^2 + (Cbe - ycalc(:,2)).^2 + (Te - ycalc(:,3)).^2;
        f = sum(f);

    end

    function dy = model(~,y,par)

        % Van de vusse reaction in a CSTR
        Ca = y(1);  Cb = y(2);  T = y(3);
        
        % The parameters to be estimated
        k1 = par(1);    k2 = par(2);    k3 = par(3);    
        E1 = par(4);    E2 = par(5);    E3 = par(6); 
        
        % reaction rates
        r1 = k1*exp(-E1/R/T)*Ca;
        r2 = k2*exp(-E2/R/T)*Cb;
        r3 = k3*exp(-E3/R/T)*Ca^2;
        
        % mass balances
        dCa = (Caf - Ca)/tau -r1 -2*r3;
        dCb = -Cb/tau + r1 - r2;
        dT  = (Tf - T)/tau -1/rho/cp*( UA/V*(T - Tk) + H1*r1 + H2*r2 +H3*r3 );
        
        % vector of derivatives
        dy = [dCa dCb dT]';

    end


end

