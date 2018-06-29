function robertsonKinetics
%% Batch reactor with Robertson's reaction system
%
%   Differential equations:
%   dy1/dt = -k1y1 + k2y2y3
%   dy2/dt = k1y1 - k2y2y3 -k3y2^2
%   dy3/dt = k3y2^2;
%
%   Initial conditions
%   t = 0 ... y1 = 1
%             y2 = 0
%             y3 = 0
% 
%   ============================================================
%   Author: ataide@peq.coppe.ufrj.br
%   Contact me for help/personal classes!

%% Problem setting

% The kinetic constants
k1 = 0.04;
k2 = 1e4;
k3 = 3e7;

% Time interval and initial conditions
tspan = [0 1e8];
y0 = [1 0 0]';

% Pass the analytic jacobian matrix to the odesolver
opt = odeset('Jacobian',@jacobian);

% Solve the model (what if you use a explicit integrator such as ode45?)
tic
[t,y] = ode15s(@model,tspan,y0,opt);
toc

% Plot the data
figured;
h = semilogx(t,y(:,1),t,y(:,2)*1e4,t,y(:,3));
set(h,'LineWidth',1.5);
ylabel('Molar fraction')
xlabel('Time')
axis([t(2) t(end) 0 1])
legend({'y_1','y_2\times 10^4','y_3'})

    % The model
    function dy = model(~,y)

        dy(1,1) = -k1*y(1) + k2*y(2)*y(3);
        dy(2,1) = k1*y(1) - k2*y(2)*y(3) - k3*y(2)^2;
        dy(3,1) = k3*y(2)^2;

    end

    function J = jacobian(~,y)

        J(1,1) = -k1;
        J(1,2) = k2*y(3);
        J(1,3) = k2*y(2);

        J(2,1) = k1;
        J(2,2) = -k2*y(3) - 2*k3*y(2);
        J(2,3) =  -k2*y(2);

        J(3,1) = 0;
        J(3,2) = 2*k3*y(2);
        J(3,3) = 0;

    end

end


