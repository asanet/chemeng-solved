function Tubular_reactor_2D
%% Transient tubular reactor in 2-D
%   
%   You'll learn:
%       +: how to discretize two dimensional domains using finite diff. 
%       +: Auto generate the sparsity pattern for your problem.
%       +: How to solve transient problems in DAE formulation
%   
%% The problem
% 
%   Differential Equation:
%   dC/dt = -d/dz( v_z - D_z*dC/dz ) + D_r/r*d/dr(r*dC/dr) - kC
%       
%   With:
%       C = C(t,r,z)
%       v_z = v_z(r)
%   
%   Boundary Conditions:
%   z = 0 ... v_z*C_f(t,r) = v_z*C(t,r,0) - D_z*dC(t,r,0)/dz 
%   z = L ... dC(t,r,L)/dz = 0
%   r = 0 ... dC(t,0,z)/dr = 0
%   r = R ... dC(t,R,z)/dz = 0
%
%   Initial condition:
%   t = 0 ... C(0,r,z) = C0(r,z) = 0
%
%   ============================================================
%   Author: ataide@peq.coppe.ufrj.br
%   homepage: github.com/asanet
%   Date: 2018-07-05
%   Matlab version: R2018a
%   Contact me for help/personal classes!

%% Problem setup
addpath('AuxFunctions')

% Reactor dimensions
R = 0.1; L = 10;

% Grid dimensions
Nr = 100;
Nz = 100;
dr = R/(Nr-1);
dz = L/(Nz-1);
r = linspace(0,R,Nr)';
z = linspace(0,L,Nz)';
n = Nr*Nz;

% Operational parameters (play around)
vmax = 1;
vz = vmax*(1 - (r/R).^2);   % Velocity profile
Dz = 0.5*0;                 % Diffusion coefficients
Dr = 0.01*0;
k = 0.1*0;                  % Kinetic constant
Cf = ones(Nr,1);            % Feed concentration

% The initial condition
y0 = zeros(n,1);
zz = zeros(n,1);
yp0 = model(0,y0,zz);

% Auto generate the sparsity pattern (don't try to run without it if Nr > 30 and/or Nz > 30)
% Try to build manually and compare with the auto-generated ;)
% 
tic
Jp = spalloc(n,n,3*n);
yr = rand(n,1);
delta = 1e-4;
per = zeros(n,1);
for ii = 1:n
    per(ii) = 1;
    f0 = model(0,yr,zz);
    f1 = model(0,yr+per*delta,zz);
    Jp(:,ii) = (f1 - f0)/delta ~= 0;  %#ok<SPRIX>
    per(ii) = 0;
end

Jyp = speye(n);
etime = toc;

fprintf('Elapsed time building the sparsity pattern = %2.4f s\n',etime);


% The solution
tf = 20;
frames = 2000;
tspan = linspace(0,tf,frames);
op = odeset('AbsTol',1e-6,'RelTol',1e-4,'JPattern',{Jp,Jyp});

tic
[t,y] = ode15i(@model,tspan,y0,yp0,op);
etime = toc;

fprintf('Elapsed time integrating the system = %2.4f s\n',etime);


% Parse the solution matrix
tic
Nt = length(t);
Csol = zeros(Nr,Nz,Nt);
for ii = 1:Nt
    Csol(:,:,ii) = reshape(y(ii,:)',Nr,Nz);
end
etime = toc;

fprintf('Elapsed time parsing the solution = %2.4f s\n',etime);

% Plot the data (to stop the animation, press CTRL+C on command window)
close all

figured;
axis([0 L 0 2*R])
caxis([0 1])
xlabel('Length')
ylabel('Diameter')
hcb = colorbar;
hcb.Label.String = 'Concentration';
colormap jet
for ii = 1:Nt
    imagesc(z,[2*r;2*r],[flip(Csol(:,:,ii));Csol(:,:,ii)])
    title(sprintf('Axial and radial profiles in t = %2.3f s',t(ii)))
    drawnow;
    pause(tf/frames/30)
end

%% The model -> DAE formulation (boundary conditions incorporated in the model)
    function res = model(~,y,yp)
        
        % Memory allocation
        res = zeros(Nr,Nz);
        C = reshape(y,Nr,Nz);
        dC = reshape(yp,Nr,Nz);
        
        % The Diff. Equations (inner points)
        for i = 2:Nr-1
            for j = 2:Nz-1
                res(i,j) = -vz(i)*( C(i,j) - C(i,j-1) )/dz + Dz*( C(i,j+1) - 2*C(i,j) + C(i,j-1) )/dz^2 + ...
                           Dr*( C(i,j) - C(i-1,j) )/2/dr/r(i) + Dr*( C(i+1,j) - 2*C(i,j) + C(i-1,j) )/dr^2  - k*C(i,j) - dC(i,j); 
            end
        end

        % In r = 0 and z = 1:Nz
        res(1,:) = C(2,:) - C(1,:);
        
        % In r = R and z = 1:Nz
        res(Nr,:) = C(Nr,:) - C(Nr-1,:);
        
        % In z = 0 and r = 2:Nr-1
        res(2:Nr-1,1) = Cf(2:Nr-1).*vz(2:Nr-1) - C(2:Nr-1,1).*vz(2:Nr-1) + Dz*( C(2:Nr-1,2) - C(2:Nr-1,1) )/dz;
        
        % In z = L and r = 2:Nr-1
        res(2:Nr-1,Nz) = C(2:Nr-1,Nz) - C(2:Nr-1,Nz-1);
        
        res = reshape(res,n,1);

    end

end