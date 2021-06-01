% This code demonstrates Laplace equation. 
clear all;  close all; clc

%Make rectangle with length L and height H
L  = 1; 
N  = 300;  dx = L/(N-1);  % size of the "infinitasimal" x for integration
x  = 0:dx:L;               % x values
xrange = [x(1) x(end)];            % x range for plot

H  = 0.5; 
N2  = 300;  dy = H/(N2-1);  
y  = 0:dy:H;               
yrange = [y(1) y(end)];        


%% Define function f and g
[xx,yy]=meshgrid(x,y);
f1      = zeros(size(xx));
f2      = f1;
f1(1:end) = 10;
f2(1:end) = 30;

g1      = zeros(size(yy));
g2      = g1;
g1(1:end) = 10;
g2(1:end) = 10;
%% Solve the solution for each boundary
n = 200; 
[u1] = LaplaceSol_rec(n,g1,1,L,H,xx,yy,dx,dy);
[u2] = LaplaceSol_rec(n,g2,2,L,H,xx,yy,dx,dy);
[u3] = LaplaceSol_rec(n,f1,3,L,H,xx,yy,dx,dy);
[u4] = LaplaceSol_rec(n,f2,4,L,H,xx,yy,dx,dy);
u=u1+u2+u3+u4;

figure;
surf(xx,yy,u)
colorbar;
xlabel('Length (x)')
ylabel('height (y)')
zlabel('heat u')
%% animation gif file
% index=1;
% fig=figure;
% 
% for n=20*(1:16)
%     [SinSeries, ~] = FSsine(n,f,L,dx,x);
%     xrange =  [x(200) x(800)]; 
%     plotFS_Gibbs(x,f,SinSeries,n,xrange); set_positionFontsAll;
%     drawnow; 
%     frame = getframe(fig);
%     im{index} = frame2im(frame);
%     index=index+1;
%     clf;
% end
% 
% figure;
% for idx=1:16
%     subplot(4,4,idx)
%     imshow(im{idx});
% end
% 
% filename = 'ConvergenceFS.gif';
% im_to_gif(filename,im,idx);