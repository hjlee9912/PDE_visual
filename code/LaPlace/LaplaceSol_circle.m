function [sol,An,Bn] = LaplaceSol_circle(K,f, a,t,dt,r)
% Solve Laplace equation on circular disk
%The solution has both cosine and sine parts
%evaluate the Fourier seires with K modes
% the coefficients are approximated by Riemann sum (can do better by integration functions in MATLAB)
sinSeries = 0; 
cosSeries = 0;

%for cos part with n=0 
    eigenf    = cos(0*t); %Eigen function 
    A0        = dot(f,eigenf,2)*dt/(2*pi); %Coefficient 

%for cos part
for n=1:K %How many "n" vales? 
    eigenf1    = cos(n*t); %Eigen function 
    An         = dot(f,eigenf1,2)*dt*(a^-n/pi); %Coefficient 
    cosSeries  = cosSeries+An.*r.^n.*eigenf1; %Cosine series
end

cosSeries=cosSeries+A0; %add n=0 term

%for sine part
for nn=1:K %How many "n" vales?
    eigenf2    = sin(nn*t); %Eigen function 
    Bn         = dot(f,eigenf2,2)*dt*(a^-nn/pi); %Coefficient 
    sinSeries  = sinSeries+Bn.*r.^nn.*eigenf2; %Sine Series
end

sol=cosSeries+sinSeries; %add all solutions
end