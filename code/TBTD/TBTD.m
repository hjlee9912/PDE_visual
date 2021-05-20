%This code demostrates the divergence of Term-By-Term Differentiation(TBTD)
%of a Fourier sine series f=x. 
%The convergence of the FS and the divergence of the TBTD of the FS are shown. 
%Written by Hojun Lee,Fei Lu


clear all;  close all; clc
add_my_paths; 

L  = .5*pi; 
N  = 1000;  dx = L/(N-1);  % size of the "infinitasimal" x for integration
x  = 0:dx:L;               % x values
xrange = [x(1) x(end)];            % x range for plot

%Fourier Sine series of f=x
f=x;
xlength=length(x);
f_prime = ones([1,xlength]);
fig=figure;
index=1;
titlee='Fourier Series of f=x';
leg='Function f';
for n=20*(1:16)
    [SinSeries, ~] = FSsine(n,f,L,dx,x);
    xrange =  [x(200) x(800)]; 
    plotFS_TBTD(x,f,SinSeries,n,xrange,xrange,titlee,leg); set_positionFontsAll;
    drawnow;
    frame = getframe(fig);
    im{index} = frame2im(frame);
    index=index+1;
    clf;
end
close;

fig2=figure;
index=1;
tbtd=0;
titlee2='TBTD of FS for f=x';
leg2="Function f'";
for n=20*(1:16)
    for k=1:n
        tbtd=tbtd+2*(-1)^(k+1)*cos(k*pi*x/L); 
        %formally differentiated sine series of x
    end
    yrange=[-25,25];
    plotFS_TBTD(x,f_prime,tbtd,n,xrange,yrange,titlee2,leg2); set_positionFontsAll;
    drawnow;
    frame2 = getframe(fig2);
    im2{index} = frame2im(frame2);
    index=index+1;
    clf;
end
close;
    
filename1 = 'ConvergenceFSx.gif';
im_to_gif(filename1,im,index-1);
filename2 = 'DivergenceTBTD.gif';
im_to_gif(filename2,im2,index-1);
