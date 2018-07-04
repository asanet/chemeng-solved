function Heat_transfer_fins
%% Steady state heat transfer in extended surfaces
% 
%   You'll learn:
%       +: How to solve a BVP with the shooting method
% 
%% The problem
% 
%   Differential Equations:
%   -kA*d^2T/dx^2 = hP*(Tinf - T)
% 
%   Boundary Conditions:
%   x = 0 ... T = Tc
%   x = L ... -k*dT/dx = 0
% 
%   ============================================================
%   Author: ataide@peq.coppe.ufrj.br
%   Contact me for help/personal classes!

%% Problem setup
addpath('AuxFunctions')

% The parameters
h = 100;
d = 0.5;
L = 2;
P = pi*d;
A = pi*d^2/4;
k = 237;
Tinf = 298.15;
Tc = 800;

% fsolve configuration
op = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','TolFun',1e-10,'TolX',1e-10, ...
                    'Display','iter-detailed','MaxIter',5000,'MaxFunEvals',10000);

% find the initial condition for u
u0_0 = 0;
u0 = fsolve(@shooting,u0_0,op);

% solve the BVP as a PVI
[x,yc] = ode15s(@model,[0 L],[Tc u0]);

% plot the profile
figured;
plot(x,yc(:,1),'LineWidth',1.5)
xlabel('Length')
ylabel('Temperature')

    function dy = model(~,y)
        
        T = y(1);
        u = y(2);
        
        dy(1,1) = u;
        dy(2,1) = -h*P/k/A*(Tinf - T);
       
    end

    function f = shooting(u0)
        
        y0 = [Tc u0]';
        [~,y] = ode15s(@model,[0 L],y0);
        
        f = y(end,2) - 0;
 
    end

end