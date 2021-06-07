% This code demonstrates wave equation on domain x[0:L]
% Generates .gif file showing change in the y displacement u for each point
% on x 
% Written by Hojun Lee Jun.04.2021
% PDE Technology Fellow project. 
%

%{
More description on the code
 - a string on [0,L] with boundary condition u(0,t)=u(L,t)=0, 
    initial condition u(x,0)=f du/dt(x,0)=g
 - traveling wave and standing wave? 
%}

add_my_paths; 
L   = 20; 
N   = 301;  dx = L/(N-1);
x   = 0:dx:L;
T0   = 20; %Tension
rho0 = 15; %mass density of the string
K=2; 
tt=30;

c   = (T0/rho0)^(1/2);

f  = @(x) sin(x);
g  = @(x) cos(x);

x_mesh = pi*x/L;

[~, ~, An_all] = computFS_coef(f,L,x,K,'sine'); 
[~, ~, Bn_all] = computFS_coef(g,L,x,K,'sine'); 

u=zeros(size(x));
for nn = 1:K
    index=1;
    fig=figure;
    for t = 1:tt
        u_t = sin(nn*x_mesh) .* (An_all(nn)*cos(nn*pi*c*t/L) * Bn_all(nn)*(L/(nn*pi*c))*sin(nn*pi*c*t/L));
        plot(x/L,u_t,'b-','LineWidth',2); 
        xlabel('x/L'); ylabel('u'); 
        xlim([0 1]); ylim([-10*10^(-4) 10*10^(-4)]);
        set_positionFontsAll;
        hold on
        drawnow;
        frame = getframe(fig);
        im{index} = frame2im(frame);
        index=index+1;
        clf;
    end
    filename = sprintf('wave n=%i.gif',nn);
    im_to_gif(filename,im,index-1);
end




