function [A,B,xgrid] = collocation(n,var,extrax)
%COLLOCATION Calculates the first and second derivatives of the Lagrange
%            interpolator at the roots of the Legendre orthogonal
%            polynomial.
%
%   [A,B,XGRID] = COLLOCATION(N,VAR,EXTRAX) N is the number of roots of the
%   orthogonal polynomial (usually the inner points of a discretization).
%   VAR is a string with the independent variable. Ex: with VAR = 'x', we 
%   are calculating the roots xi in [-1, 1] such as P(x) = 0.
%   With VAR = '2*x-1', we are calculating the roots xi i [0, 1] such as 
%   P(2x-1) = 0. EXTRAX is a 2-by-1 optional vector which pre-append and
%   post-append xgrig by EXTRAX(1) and EXTRAX(2), respectively.
% 
%   Examples:See the Furnace_wall.m file.
% 
%   ============================================================
%   Author: ataide@peq.coppe.ufrj.br
%   homepage: github.com/asanet
%   Contact me for help/personal classes!

if nargin < 3
    extrax = [];
elseif nargin < 2
    var = 'x';
    extrax = [];
elseif nargin < 1
    error('Missing the number of roots.')
end

if n < 1
    n = 1;
    warning('Minimum value for n is 1')
end

if ~ischar(var) || ~any(var == 'x')
    error('Second argument must be a string and contain "x".')
end

% Symbolic variable
syms('x')
var = str2sym(var);

% Polynomial in var
LP = legendreP(n,var); 

% Roots of the polynomial
xgrid = double(vpasolve(LP == 0));

if ~isempty(extrax)
    if length(extrax) ~= 2
        error('extrax must have size two');
    else
        xgrid = [extrax(1); xgrid; extrax(2)];
        n = n + 2;
    end
end

% Matrix A
A = zeros(n,n);
v = zeros(n,1);
for i = 1:n
    p = 1;
    for j = 1:n
        fat = xgrid(i) - xgrid(j);
        v(i) = fat*v(i) + p;
        p = fat*p;
    end
end

for i = 1:n
    for j = 1:n
        if j ~= i
            A(i,j) = v(i)/v(j)/(xgrid(i)-xgrid(j));
            A(i,i) = A(i,i) - A(i,j);
        end
    end
end

B = A^2;





