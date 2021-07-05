%This code solves Regular Sturm-Liouville Eigenvalue Problem 
%[p(x)*y']'+q(x)y = lambda*r(x)*y
%where y(x) is eigenfunction and lambda is eigenvalue with domain [a,b]
%This example assumes Dirichlet BC with BC y(a) = 0, y(b) = 0
%uses function 'SLEP' and outputs a plot of: The first three eigenfunctions vs x
%needs functions p,q,and r
%By Hojun Lee and Fei Lu 
%Jul/4/2021

b = pi;
a = 0;
N = 99;  %number of x points excluding both ends
dx = (b-a)/(N+1); %"infitasimal" x
x = a:dx:b; %domain
x_in = x(2:end-1); %x values excluding both ends

[P,dPdx] = p(x_in); %evaluate p(x) and p'(x)
Q = q(x_in); %evaluate Q(x)
R = r(x_in); %evaluate R(x)
[v,e]=SLEP(N,x_in,P,dPdx,Q,R); %solve SLEP

fig=figure;
plot(x',v(:,1:3),'linewidth',2)
darkBackground(fig,[0 0 0],[1 1 1]); set_positionFontsAll;
xlabel('\color{white}x')
legend('\color{white}\phi_1', '\color{white}\phi_2','\color{white}\phi_3')



