
%This code demostrates a wave equation with spring-mass system
%This code contains error. Development still in progress. (6/13/2021)
%Solution to the PDE is uploaded in OneDrive. 
add_my_paths;
L=pi;  %length
tt=10; %time
K=10;
N   = 301;  dx = L/(N-1);
x   = 0:dx:L;%domain

index=1;
fig=figure;

for t=1:tt  
    v=0;
    for n=1:K
        eigenf=sin(n*pi*x/L);
        w = sin(t) + (x/pi)*-1*sin(t);
    
        syms a(t)
        Da = diff(a,t);
        ode = diff(a,t,2)+n^2*diff(a,t) == 2*(1/n-sin(pi*n)/(pi*n^2))*sin(t)/pi;
        cond1 = a(0) == (-2*sin(pi*n))/(pi*(n^2-1));
        cond2 = Da(0) == 2*((n^2-1)*sin(pi*n)+(pi*n^3*cos(pi*n)+pi*n))/(pi^2*n^2*(n^2-1));
        conds = [cond1 cond2];
        aSol(t) = dsolve(ode,conds);
        
        v=v+aSol.*eigenf;
    end
    u=v+w;
    plot(x,u,'b-','LineWidth',2); 
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
