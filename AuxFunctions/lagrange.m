function P = lagrange(xg,yg,x)
%LAGRANGE Gives the lagrange interpolator polynomial.
%
%   P = LAGRANGE(XG,YG,X) Calculates the polynomial P in X with the nodal
%   points given in XG and coefficients in YG. XG and YG must be of the
%   same size.
% 
%   Examples:See the Furnace_wall.m file.
% 
%   ============================================================
%   Author: ataide@peq.coppe.ufrj.br
%   homepage: github.com/asanet
%   Date: 2018-07-05
%   Matlab version: R2018a
%   Contact me for help/personal classes!

% Error handling
if length(xg) ~= length(yg)
    error('xg e yg must have the same lengths!')
end

% Put x in a column vector
x = x(:);

% Get the dimensions
np = length(x);
ng = length(xg);

% Allocation
P = zeros(np,1);

% Calculate the lagrange polynomial in x -> P(x)
for j = 1:ng
    l = ones(np,1);
    for i = 1:ng
        if i ~= j
            l = l.*( x - xg(i) )/( xg(j) - xg(i) );
        end
    end
    P = P + l*yg(j);
end
        
