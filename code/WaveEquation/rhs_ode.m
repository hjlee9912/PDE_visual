function dadt = rhs_ode(t,a,nn, qnt)
% right hand side of ODE a''  = -n^2 a + qnt
dadt(1) = a(2);
dadt(2) = -nn^2*a(1) + qnt(t);
dadt    = dadt'; 
end