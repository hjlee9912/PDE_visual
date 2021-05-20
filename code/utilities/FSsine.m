function [SinSeries, Bn_all] = FSsine(K,f, L,dx,x_mesh)
% evaluate the Fourier sine seires with K modes
% the coefficients are approximated by Riemann sum (can do better by integration functions in MATLAB)
SinSeries = 0; 
Bn_all    = zeros(1,K);
x_mesh    = x_mesh*(pi/L); 
for n=1:K %How many "n" vales? Sine seires starts with n=1
    eigenf     = sin(n*x_mesh); %Eigen function for the Fourier Sine Series
    Bn         = dot(f,eigenf)*dx*2/L; %Coefficient for the Fourier Sine Series
    SinSeries  = SinSeries+Bn*eigenf; %Fourier Sine Series
    Bn_all(n)  = Bn; 
end

end