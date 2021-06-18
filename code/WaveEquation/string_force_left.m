% This code demostrates the motion of a string that is fixed on the right end, 
% but has a force on the left end. 
%  --  Wave eqaution 
% Future extension: 
%  - kinetic energy and potential energy of the system; exe 4.4.9
%   (since energy is pumped in and the system does not disspate energy, the solution will blow-up in time)
%  - Issue to be address: discontinuity from IC 
% Written by Fei Lu, Hojun Lee; 2021/6/18

%{
================  Algorithm decription ================
 IBVP for 1D wave equation
    u_{tt}  = u_{xx} 
 BC: u(0,t) = F(t);    u(pi,t)  = 0;    % fixed right bdry
 IC: u(x,0) = f(x);    u_t(x,0) = g(x). 

1. Solution via homogenization of BC  u = v + w
   w(x,t) = (1-x/pi)*F(t)
     IBVP v
        v_{tt}  = v_{xx} - w_{xx} + w_t 
     BC: v(0,t) = 0;    v(pi,t)  = 0;    %  homogeneour BC
     IC: v(x,0) = f(x) - w(x,0);    v_t(x,0) = g(x) - w_t(x,t). 
2. solve v by eigenfunction expansion method: sin(nx)
%}

add_my_paths;

%% settings
L     = pi;           % length
tt    = 40;           % time
K     = 50;            % "N" value for fourier sine series

%% BC & IC
kBC  = 1; 
F = @(t) sin(kBC*t); dtF = @(t) kBC*cos(kBC*t);  dttF = @(t) - kBC^2*sin(kBC*t); 
f = @(x) sin(x);  
g = @(x) sin(4*x); 

f_v = @(x) f(x) - (1-x/L).*F(0); 
g_v = @(x) g(x) - (1-x/L).*dtF(0); 
% Q_v = @(x,t) -(1-x/pi).*dttF(t); 

%% space time mesh
NN    = 101;  dx = L/(NN-1);   % infinitasimal x
dt    = tt/400;             % infinitasimal t
x     = 0:dx:L;                % domain
tspan = 0:dt:tt;               % time

%% solution via eigenfunction expansion method; 
v=0; %u=v+w 
[~, an0]     = fs_sine(K,f_v, L,dx/10);          % IC coefficients a_n(0) 
[~, an_dt_0] = fs_sine(K,g_v, L,dx/10);          % IC coefficients a_n'(0)
[~,qn]       = fs_sine(K,@(x)(1-x/L), L,dx/10);  % qnt  = @(t) - qn(nn)*dttF(t); 


for nn = 1:K
    % solve 2nd order diffeq
    ic_ode  = [an0(nn), an_dt_0(nn)]; 
    qnt     = @(t)  - qn(nn)*dttF(t); 
    [t,a]   = ode45( @(t,a) rhs_ode(t,a,nn,qnt), tspan, ic_ode); % solve for a_n(t)
    eigenf = sin(nn*pi*x/L);    % eigenfunction
    v     = v+ a(:,1)*eigenf;
end

w  = sin(tspan)'*( 1- (x/pi)); 
u  = v+w;         %final solution


%% plot the solution 
index = 1;
fig   = figure;
for i=1:NN
    plot(x,u(i,:),'b-','LineWidth',2); 
    xlabel('x'); ylabel('u');   axis([0,pi,-1.5,1.5])
    set_positionFontsAll;
    drawnow;
    frame     = getframe(fig);
    im{index} = frame2im(frame);
    index     = index+1;
    clf;
end
filename = sprintf('string_leftForce.gif');
im_to_gif(filename,im,index-1); %make .gif file for each n value
