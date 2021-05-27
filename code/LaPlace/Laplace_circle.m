%This code demonstrates Laplace equation on circular disk with radius a=1
%BC = f(a,theta)=30

theta = 2*pi;
N  = 300;  dt = theta/(N-1);  % size of the "infinitasimal" t for integration
t  = -pi:dt:pi; 

a=1; %radius
dr = a/(N-1);  % size of the "infinitasimal" r for integration
r  = 0:dr:a; 

[tt,rr]=meshgrid(t,r); %create meshgrid for 3D plot
[xx,yy]=pol2cart(tt,rr); %convert polar to cartesian
f=zeros(size(tt));
f(1:end)=30; %boundary condition f

n=300; %number of summation
[sol,An,Bn] = LaplaceSol_circle(n,f, a,tt,dt,rr); %calculate solution

surf(xx,yy,sol) %plot 3D
