function MPC_tank_reactor
%% Model predictive control of a tank reactor
%   
%   You'll learn:
%       +: how to solve optimization problems
%       +: How to apply a model predictice control with non linear models
%   
%% The problem
%   
%   minimize the integral of (y_pred - y_sp)^2
%   subject to: F(t,y,y') = 0
%               lb <   u  < ub
%                0 < |du| < du_max
%               
%   ============================================================
%   Author: ataide@peq.coppe.ufrj.br
%   homepage: github.com/asanet
%   Date: 2018-07-05
%   Matlab version: R2018a
%   Contact me for help/personal classes!

%% Problem setup
addpath('AuxFunctions')

% The model parameters
k1 = 3.575e8;       k2 = 3.575e8;       k3 = 2.5119e3;
E1 = 8.2101e4;      E2 = 8.2101e4;      E3 = 7.1172e4;
H1 = 4.2e3;         H2 = -11e3;         H3 = -41.85e3;
rho = 0.9342e3;     cp = 3.01e3;        V = 10e-3;      
tau = 80;           Tf = 403.15;        Caf = 1000;
UA = 0.215*1120;    R = 8.3145;

% The objective function parameters
Hp = 1000;      % Prediction horizon
nt = 50;        % Time sampling points
P_sp = 1;       % Set point violetion penalty
P_u = 1e6;      % Saturation penalty
P_du = 1e6;     % Delta u penalty
lb = 273.15;    % lower bound to u
ub = 373.15;    % upper bound to u
du_max = 10;    % max Delta u
p = 20;         % The regularization parameter

% The times for control actions (the size defines the control horizon)
td = [0 10 20 40]';

% The control horizon
% Hc = length(td);

% The Sampling time
Ta = 100;

% The set point
spval = 380;
T_sp = spval*ones(nt,1);

% The current plant state (starting condition)
ynow = [0 0 300]';

% Simulation span
tspan = 1:3*Hp;

% The intial u profile (must have Hc values)
u = [373.15 0 0 0]';
umi = u(1);

% The solver configuration
simulationOpt = odeset('AbsTol',1e-4,'RelTol',1e-3);

% use with fsolve
% controlOpt = optimoptions(@fsolve,'Algorithm','levenberg-marquardt',...
%                             'TolFun',1e-8,'TolX',1e-8,'MaxIter',1e2,...
%                             'MaxFunEvals',1e5,'Display','none');

% use with fminsearch
controlOpt = optimset('TolFun',1e-8,'TolX',1e-8,'MaxIter',100,'Display','none');
                        
% The problem call
nsim = length(tspan);
nsamples = fix(nsim/Ta);

tm = zeros(nsamples,1);
ym = zeros(nsamples,3);
um = zeros(nsamples,1);

ym(1,:) = ynow;
um(1) = u(1);

close all
for i = 1:nsamples-1

    fprintf('Sampling time %d of %d\n',i,nsamples-1)

    % call the controller to find the u predicted (either fsolve or fminsearch)
    %u = fsolve(@mpc_cost_function,u,controlOpt,ynow);
    u = fminsearch(@mpc_cost_function,u,controlOpt,ynow);
    
    % Apply only the first control action
    u(2:end) = 0;
    
    % call the plant with updated u and return t and y measured
    [tmi,ymi] = ode15s(@plant, [tm(i) tm(i)+Ta], ynow, simulationOpt, u);
    
    % Updated the measured variables
    ynow = ymi(end,:)';
    
    tm(i+1) = tmi(end);
    ym(i+1,:) = ymi(end,:);
    um(i+1) = u(1);
    umi = u(1);
    
    % Disturbance in tal
    if i == fix(nsamples/2)
        Tf = 420;
    end
end

%% plot the data
figure;
subplot(211)
h = plot(tm,ym(:,3),tm,spval*ones(nsamples,1),'LineWidth',1.5);
set(h(2),'LineStyle',':')
ylabel('Temperature')
title('Controlled variable')
set(gca,'ygrid','on','xgrid','on','fontsize',16)
legend({'measured','setpoint'},'location','southeast')

subplot(212)
h = plot(tm,um,tm,lb*ones(nsamples,1),tm,ub*ones(nsamples,1),'LineWidth',1.5);
set(h(2:3),'LineStyle',':','Color','k')
ylabel('Jacket temperature')
xlabel('Time')
title('Manipulated variable')
legend({'measured','bounds'},'location','southwest')
set(gca,'ygrid','on','xgrid','on','fontsize',16)

    function f = mpc_cost_function(u,ynow)
        
        % Solve the model until the preditction horizon
        ts1 = linspace(0,Hp,nt); 
        [t,y] = ode15s(@model,ts1,ynow,simulationOpt,u);
             
        % The set point tracking term
        f1 = trapz(t,( y(:,3) - T_sp ).^2 );
        
        % The control variables bound
        min_u = min(cumsum(u));
        max_u = max(cumsum(u));
        f2 = -min(0, ub - max_u) -min(0, min_u - lb);
        
        % Delta u penalty
        f3 = -min(0, du_max - max(abs([umi- u(1); u(2:end)])));
        
        % THe objective function
        f = P_sp*f1 + P_u*f2 + P_du*f3;
    end

    function dy = model(t,y,u)
        % Van de vusse reaction in a CSTR
        Ca = y(1);  Cb = y(2);  T = y(3);
        
        % the regularization function
        reg = 1/2 + tanh(p*(t - td))/2;
        
        % control variable
        Tc = sum(u.*reg);
        
        % reaction rates
        r1 = k1*exp(-E1/R/T)*Ca;
        r2 = k2*exp(-E2/R/T)*Cb;
        r3 = k3*exp(-E3/R/T)*Ca^2;
        
        % mass balances
        dCa = (Caf - Ca)/tau -r1 -2*r3;
        dCb = -Cb/tau + r1 - r2;
        dT  = (Tf - T)/tau -1/rho/cp*( UA/V*(T - Tc) + H1*r1 + H2*r2 +H3*r3 );
        
        % vector of derivatives
        dy = [dCa dCb dT]';

    end

    function dy = plant(t,y,u)
        % The virtual plant: in this case, we consider a perfect model,
        % that is, plant = model
        dy = model(t,y,u);
        
    end
end


