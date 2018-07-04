function Pl = lagrange(xg,yg,x)

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
Pl = zeros(np,1);

% Calculate the lagrange polinomial in x -> P(x)
for j = 1:ng
    l = ones(np,1);
    for i = 1:ng
        if i ~= j
            l = l.*( x - xg(i) )/( xg(j) - xg(i) );
        end
    end
    Pl = Pl + l*yg(j);
end
        
