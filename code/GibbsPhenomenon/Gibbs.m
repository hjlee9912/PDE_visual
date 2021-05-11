%This code demonstrates Gibbs Phenomenon. 
%The Fourier Sine Series of f is used 
%where f=0 @[0,pi/2] & f=1 @[pi/2,pi] & f=-1 @[pi,1.5pi]
%Difference in the frequencies at n=100 and n=300 is shown,
%as well as the difference in the overshoot at the two jump discontinuities. 
%Written by Hojun Lee 

clear all
close all
clc

L=1.5*pi; 
N=1000; 
dx=L/(N-1); %size of the "infinitasimal" x
x=0:dx:L; %x values

f=zeros(size(x)); 
bound1=ceil(length(x)/3); %First jump discontinuity at x=pi/2
bound2=ceil(length(x)/3*2); %Second jump discontinuity at x=pi
f(bound1:bound2)=1; 
f(bound2+1:end)=-1; %Define f
SinSeries=0;

%Now calculate the Fourier sine series with n=100
for n=1:100 %How many "n" vales? Sine seires starts with n=1
    eigenf=sin(n*pi*x/L); %Eigen function for the Fourier Sine Series
    Bn=dot(f,eigenf)*dx*2/L; %Coefficient for the Fourier Sine Series
    SinSeries=SinSeries+Bn*eigenf; %Fourier Sine Series
end

%Gibbs phenomenon with n=100
figure
plot(x,f,'k','LineWidth',4)
hold on
plot(x,SinSeries,'r-','LineWidth',2)
xlim([0.1 4.5])
title('Gibbs Phenomenon (n=100)')

%Now calculate the Fourier sine series with n=300 to show the difference
SinSeries=0;
for n=1:300 %How many "n" vales? Sine seires starts with n=1
    eigenf=sin(n*pi*x/L); %Eigen function for the Fourier Sine Series
    Bn=dot(f,eigenf)*dx*2/L; %Coefficient for the Fourier Sine Series
    SinSeries=SinSeries+Bn*eigenf; %Fourier Sine Series
end

%Gibbs phenomenon with n=300
figure
plot(x,f,'k','LineWidth',4)
hold on
plot(x,SinSeries,'m','LineWidth',2)
xlim([0.1 4.5])
title('Gibbs Phenomenon (n=300)')

%Zoomed in at the first jump discontinuity (n=300)
figure
plot(x,f,'k','LineWidth',4)
hold on
plot(x,SinSeries,'g-','LineWidth',2)
xlim([x(bound1-50) x(bound1+50)])
title('at L=pi/2 (n=300)')

%Zoomed in at the second jump discontinuity (n=300)
figure
plot(x,f,'k','LineWidth',4)
hold on
plot(x,SinSeries,'c-','LineWidth',2)
xlim([x(bound2-50) x(bound2+50)])
title('at L=pi (n=300)')


