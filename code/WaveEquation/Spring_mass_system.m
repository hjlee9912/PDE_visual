%This code demostrates a wave equation with spring-mass system
%with length L
%Solution to the PDE is uploaded in OneDrive. 
add_my_paths;
L=pi;  %length
tt=10; %time
K=5; %"N" value for fourier sine series
N   = 301;  dx = L/(N-1); %infinitasimal x
x   = 0:dx:L;%domain

index=1;
fig=figure;

for time=1:0.2:tt 
    v=0;
    for n=2:K %error: n=1 does not work in this case. Why? (6/14/2021)
        eigenf=sin(n*pi*x/L); %eigenfunction
        w = sin(time) + (x/pi)*-1*sin(time);
        %solve 2nd order diffeq
        syms a(t)
        Da = diff(a,t);
        ode = diff(a,t,2)+n^2*diff(a,t) == 2*(1/n-sin(pi*n)/(pi*n^2))*sin(t)/pi;
        cond1 = a(0) == (-2*sin(pi*n))/(pi*(n^2-1));
        cond2 = Da(0) == 2*((n^2-1)*sin(pi*n)+(pi*n^3*cos(pi*n)+pi*n))/(pi^2*n^2*(n^2-1));
        conds = [cond1 cond2];
        aSol(t) = dsolve(ode,conds);
        aSolution=matlabFunction(aSol);
        
        
        v=v+aSolution(time)*eigenf;
    end
    u=v+w; %final solution
    plot(x,u,'b-','LineWidth',2); 
    xlabel('x'); ylabel('u'); 
    xlim([-0.5 L]); ylim([-1.3 1.3]);
    set_positionFontsAll;
    drawnow;
    frame = getframe(fig);
    im{index} = frame2im(frame);
    index=index+1;
    clf;
end
filename = sprintf('springmass.gif');
im_to_gif(filename,im,index-1); %make .gif file for each n value
