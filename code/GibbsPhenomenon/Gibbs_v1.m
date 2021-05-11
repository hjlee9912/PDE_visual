% This code demonstrates Gibbs Phenomenon. 
%The Fourier Sine Series of f is used 
%where f=0 @[0,pi/2] & f=1 @[pi/2,pi] & f=-1 @[pi,1.5pi]
%Difference in the frequencies at n=100 and n=300 is shown,
%as well as the difference in the overshoot at the two jump discontinuities. 
%Written by Hojun Lee 

%{
%%%%% TODO (in summer)
- make it interactive;
- animation: n increases 
%}

clear all;  close all; clc

L  = 1.5*pi; 
N  = 1000;  dx = L/(N-1);  % size of the "infinitasimal" x for integration
x  = 0:dx:L;               % x values


%% Define function f
f      = zeros(size(x)); 
bound1 = ceil(length(x)/3);   % First jump discontinuity at x=pi/2
bound2 = ceil(length(x)/3*2); % Second jump discontinuity at x=pi
f(bound1:bound2) = 1; 
f(bound2+1:end)  = -1; %Define f


%% Gibbs phenomenon with n=100
n = 100; 
[SinSeries, ~] = FSsine(n,f,L,dx,x);
xrange = [0.1 4.5]; 
plotFS_Gibbs(x,f,SinSeries,n,xrange); 

%% Gibbs phenomenon with n=300
n = 300; 
[SinSeries, ~] = FSsine(n,f, L,dx,x);
plotFS_Gibbs(x,f,SinSeries,n,xrange); 

% Zoomed in at the first jump discontinuity (n=300)
figure
plot(x,f,'k','LineWidth',2); hold on
plot(x,SinSeries,'g-','LineWidth',1);
xlim([x(bound1-50) x(bound1+50)]); 
title('at L=pi/2 (n=300)')

% Zoomed in at the second jump discontinuity (n=300)
figure
plot(x,f,'k','LineWidth',2);          hold on
plot(x,SinSeries,'c-','LineWidth',1); 
xlim([x(bound2-50) x(bound2+50)]);
title('at L=pi (n=300)')


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

function plotFS_Gibbs(x,f,SinSeries,n,xrange)
% plot the Gibbs phenomenon
figure; 
plot(x,f,'k','LineWidth',2);          hold on; 
plot(x,SinSeries,'g-','LineWidth',1);
yline( (max(f)-min(f))*0.09 + max(f)); 
xlim(xrange); 
legend('Function f','Fourier Series','9% overshot line')
title(sprintf('Gibbs Phenomenon (n=%i)',n) );
end