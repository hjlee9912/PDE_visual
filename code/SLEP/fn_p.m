function [y,dydx]= fn_p(x)

factor = 1; 
y    = sin(factor*x)+2;
dydx = factor*cos(factor*x);

end
