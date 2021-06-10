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
% - wave traveling to the right an left.
% - Superposition of the traveling waves
% - sinewave, triangular wave, and rectangular wave
%% Set basic parameters

L   = 10; %length of the string
N   = 301;  dx = L/(N-1);
x   = 0:dx:L;%domain
T0   = 10; %Tension
rho0 = 5; %mass density of the string
K=3; %"n" value at fourier series. 
tt=30; %time
c   = (T0/rho0)^(1/2);

%set functions
f  = @(x) sin(x);
g  = @(x) cos(x);
h  = @(x) sawtooth(x);
p  = @(x) square(x);

%% Evaluate standing wave by solving Wave Equation (PDE)

x_mesh = pi*x/L;

[~, ~, An_all] = computFS_coef(f,L,x,K,'sine'); %compute coefficients
[~, ~, Bn_all] = computFS_coef(g,L,x,K,'sine'); 

for nn = 1:K
    index=1;
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
    im_to_gif(filename,im,index-1); %make .gif file for each n value
end
%% superpose two traveling waves. (show cosine components of f)

% for nn = 3:4
%     index=1;
%     fig=figure;
%     for t = 1:tt
%         u_1 = 0.5*f((x-c*t)/L);
%         u_2 = 0.5*f((x+c*t)/L);
%         u_tR = 0.5*cos(nn*pi*(x-c*t)/L); %wave traveling to the right
%         u_tL = 0.5*cos(nn*pi*(x+c*t)/L); %Wave traveling to the left
%         u_tS = u_tR-u_tL; %Superposition of the two waves
%         plot(x/L,u_tR,'r-','LineWidth',2); 
%         hold on
%         plot(x/L,u_tL,'g-','LineWidth',2);
%         hold on
%         plot(x/L,u_tS,'k-','LineWidth',2);
%         xlabel('x/L'); ylabel('u'); 
%         xlim([0 1]); ylim([-1 1]);
%         set_positionFontsAll;
%         drawnow;
%         frame = getframe(fig);
%         im{index} = frame2im(frame);
%         index=index+1;
%         clf;
%     end
%     filename = sprintf('travelingwave n=%i.gif',nn);
%     im_to_gif(filename,im,index-1);
% end
%% superpose two travleing sine waves. 

    index=1;
    fig=figure;
    for t = 1:0.2:tt
        u_1 = 0.5*f((x-c*t));
        u_2 = 0.5*f((x+c*t));
        u_tS = u_1-u_2; %Superposition of the two waves
        plot(x,u_1,'r-','LineWidth',2); 
        hold on
        plot(x,u_2,'b-','LineWidth',2,'MarkerIndices',[100],'MarkerSize',30);
        hold on
        plot(x,u_tS,'m-p','LineWidth',2,'MarkerIndices',[100],'MarkerSize',30);
        xlabel('x'); ylabel('u'); 
        ylim([-1 1]);
        legend('Traveling wave 1','Traveling wave 2','Standing wave','location','southeast')
        set_positionFontsAll;
        drawnow;
        frame = getframe(fig);
        im{index} = frame2im(frame);
        index=index+1;
        clf;
    end
    filename = sprintf('travelingwave_sine.gif');
    im_to_gif(filename,im,index-1);
    %% superpose two traveling triangular waves. 
    
    index=1;
    fig=figure;    
    for t = 1:0.2:tt
        u_1 = 0.5*h((x-c*t));
        u_2 = 0.5*h((x+c*t));
        u_tS = u_1-u_2; %Superposition of the two waves
        plot(x,u_1,'r-','LineWidth',2); 
        hold on
        plot(x,u_2,'b-','LineWidth',2);
        hold on
        plot(x,u_tS,'k-p','LineWidth',2,'MarkerIndices',[100],'MarkerSize',30);
        xlabel('x'); ylabel('u'); 
        ylim([-1 1])
        legend('Traveling wave 1','Traveling wave 2','Standing wave','location','southeast')
        set_positionFontsAll;
        drawnow;
        frame = getframe(fig);
        im{index} = frame2im(frame);
        index=index+1;
        clf;
    end
    filename = sprintf('travelingwave_sawtooth.gif');
    im_to_gif(filename,im,index-1);
    %% superpose two traveling square waves. 
    
    index=1;
    fig=figure;    
    for t = 1:0.2:tt
        u_1 = 0.5*p((x-c*t));
        u_2 = 0.5*p((x+c*t));
        u_tS = u_1-u_2; %Superposition of the two waves
        plot(x,u_1,'r-','LineWidth',2); 
        hold on
        plot(x,u_2,'b-','LineWidth',2);
        hold on
        plot(x,u_tS,'k-p','LineWidth',2,'MarkerIndices',[100],'MarkerSize',30);
        xlabel('x'); ylabel('u'); 
        ylim([-1 1])
        legend('Traveling wave 1','Traveling wave 2','Standing wave','location','southeast')
        set_positionFontsAll;
        drawnow;
        frame = getframe(fig);
        im{index} = frame2im(frame);
        index=index+1;
        clf;
    end
    filename = sprintf('travelingwave_square.gif');
    im_to_gif(filename,im,index-1);




