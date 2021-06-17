function dadt = waveCoeff(t,a,nn)
dadt(1) = a(2);
dadt(2) = 2*(1/nn-sin(pi*nn)/(pi*nn^2))*sin(t)/pi-nn^2*a(2);
dadt=dadt';
end