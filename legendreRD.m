function [A,B,xgrid] = Collocation(n,var,extrax)
% Compute Pn(var), the Legendre polynomial of degree n in var(x) variable.
% Output: Roots and derivative matrixes

if nargin < 3
    extrax = [];
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
    xgrid = [extrax(1); xgrid; extrax(2)];
    n = n + 2;
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





