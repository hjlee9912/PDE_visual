% This code demonstrates wave equation on domain x[0:L]
% Generates .gif file showing change in the y displacement u for each point
% on x 
%Written by Hojun Lee Jun.04.2021
%PDE Technology Fellow project. 
% 
%{standing wave:
% - a string on [0,L] with boundary condition u(0,t)=u(L,t)=0, 
%    initial condition u(x,0)=f du/dt(x,0)=g
% - Has n-1 nodes
% traveling wave:
% - wave traveling to the right and wave traveling to the left.
% - Superposition of the two traveling waves
%}

L   = 20; 
N   = 301;  dx = L/(N-1);
x   = 0:dx:L;
T0   = 20; %Tension
rho0 = 15; %mass density of the string
K=2; 
tt=40;

c   = (T0/rho0)^(1/2);

f  = @(x) sin(x);
g  = @(x) cos(x);

x_mesh = pi*x/L;

[~, ~, An_all] = computFS_coef(f,L,x,K,'sine'); 
[~, ~, Bn_all] = computFS_coef(g,L,x,K,'sine'); 

for nn = 1:K
    index=1;
    index2=1;
    fig=figure;
    for t = 1:tt %time
        %Standing wave 
        u_t = sin(nn*x_mesh) .* (An_all(nn)*cos(nn*pi*c*t/L) * Bn_all(nn)*(L/(nn*pi*c))*sin(nn*pi*c*t/L));
        plot(x/L,u_t,'b-','LineWidth',2); 
        xlabel('x/L'); ylabel('u'); 
        xlim([0 1]); ylim([-40*10^(-4) 40*10^(-4)]);
        set_positionFontsAll;
        drawnow;
        frame = getframe(fig);
        im{index} = frame2im(frame);
        index=index+1;
        clf;
    end
    filename = sprintf('wave n=%i.gif',nn);
    im_to_gif(filename,im,index-1);
end

for nn = 1:K
    index=1;
    fig=figure;
    for t = 1:tt
        u_tR = 0.5*cos(nn*pi*(x-c*t)/L); %wave traveling to the right
        u_tL = 0.5*cos(nn*pi*(x+c*t)/L); %Wave traveling to the left
        u_tS = u_tR-u_tL; %Superposition of the two waves
        plot(x/L,u_tR,'r-','LineWidth',2); 
        hold on
        plot(x/L,u_tL,'r-','LineWidth',2);
        hold on
        plot(x/L,u_tS,'k-','LineWidth',2);
        xlabel('x/L'); ylabel('u'); 
        xlim([0 1]); ylim([-1 1]);
        set_positionFontsAll;
        drawnow;
        frame = getframe(fig);
        im{index} = frame2im(frame);
        index=index+1;
        clf;
    end
    filename = sprintf('travelingwave n=%i.gif',nn);
    im_to_gif(filename,im,index-1);
end




