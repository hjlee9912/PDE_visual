%This code demostrates a wave equation with spring-mass system
%using ode45
%with length L
%Solution to the PDE is uploaded in OneDrive. 


%contains error



add_my_paths;
L=pi;  %length
tt=10; %time
K=3; %"N" value for fourier sine series
NN   = 101;  dx = L/(NN-1); %infinitasimal x
dt = tt/(NN-1); %infinitasimal t
x   = 0:dx:L;%domain
tspan   = 0:dt:tt; %time
[X,T]=meshgrid(x,tspan); %generate matrix


index=1;
fig=figure;

v=0; %u=v+w
for nn=1:K
    eigenf=sin(nn*pi*X/L); %eigenfunction
    w = sin(T) + (X/pi)* -1 * sin(T);
    %solve 2nd order diffeq
    if nn==1
        ic=[1,0]; %IC when n=1
    else
        ic=[(-2*sin(pi*nn))/(pi*(nn^2-1)),2*nn*(cos(pi*nn)+1)/(pi*(nn^2-1))];
    end
    [t,a] = ode45(@(t,a) waveCoeff(t,a,nn), tspan, ic); %solve for a, the coefficient for Fourier Sine Series
    v=v+a(:,1)'.*eigenf;
end
u=v+w; %final solution
for i=1:NN
    plot(x,u(i,:),'b-','LineWidth',2); 
    xlabel('x'); ylabel('u'); 
    set_positionFontsAll;
    drawnow;
    frame = getframe(fig);
    im{index} = frame2im(frame);
    index=index+1;
    clf;
end
filename = sprintf('springmass.gif');
im_to_gif(filename,im,index-1); %make .gif file for each n value
