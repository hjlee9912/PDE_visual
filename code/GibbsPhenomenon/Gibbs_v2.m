% This code demonstrates Gibbs Phenomenon. 
% The Fourier Sine Series of f is used 
% where f=0 @[0,pi/2] & f=1 @[pi/2,pi] & f=-1 @[pi,1.5pi]
% Difference in the frequencies at n=100 and n=300 is shown,
% as well as the difference in the overshoot at the two jump discontinuities. 
% Written by Hojun Lee, Fei Lu

%{
%%%%% TODO (in summer)
- make it interactive;
- animation: n increases 
%}

clear all;  close all; clc

L  = 1.5*pi; 
N  = 1000;  dx = L/(N-1);  % size of the "infinitasimal" x for integration
x  = 0:dx:L;               % x values
xrange = [x(1) x(end)];            % x range for plot

%% Define function f
f      = zeros(size(x)); 
bound1 = ceil(length(x)/3);   % First jump discontinuity at x=pi/2
bound2 = ceil(length(x)/3*2); % Second jump discontinuity at x=pi
f(bound1:bound2) = 1; 
f(bound2+1:end)  = -1; %Define f


%% Gibbs phenomenon with n=100
n = 100; 
[SinSeries, Bn_all] = FSsine(n,f,L,dx,x);
figure; plotFS_Gibbs(x,f,SinSeries,n,xrange); 

%% Gibbs phenomenon with n=320
n = 320; 
[SinSeries, Bn_all] = FSsine(n,f, L,dx,x);
figure; plotFS_Gibbs(x,f,SinSeries,n,xrange); 

% Zoomed in at the first jump discontinuity (n=300)
figure
plot(x,f,'k','LineWidth',2); hold on
plot(x,SinSeries,'g-','LineWidth',1);
xlim([x(bound1-50) x(bound1+50)]); 
title('at L=pi/2 (n=300)')

% Zoomed in at the second jump discontinuity (n=300)
figure; 
plot(x,f,'k','LineWidth',2);          hold on
plot(x,SinSeries,'c-','LineWidth',1); 
xlim([x(bound2-50) x(bound2+50)]);
title('at L=pi (n=300)')


FS_l2 =  cumsum(Bn_all.^2); % mean squares
f_l2  = FS_l2(end); 
errL2 = f_l2- FS_l2; 

figure
% plot(FS_l2,'linewidth',2);  hold on; 
plot(errL2,'linewidth',2); 


%% animation gif file
index=1;
fig=figure;

for n=20*(1:16)
    [SinSeries, ~] = FSsine(n,f,L,dx,x);
    xrange =  [x(200) x(800)]; 
    plotFS_Gibbs(x,f,SinSeries,n,xrange); set_positionFontsAll;
    darkBackground(fig,[0 0 0],[1 1 1]);
    drawnow; 
    frame = getframe(fig);
    im{index} = frame2im(frame);
    index=index+1;
    clf;
end

figure;
for idx=1:16
    subplot(4,4,idx)
    imshow(im{idx});
end

filename = 'ConvergenceFS.gif';
im_to_gif(filename,im,idx);
%% pointwise convergence

nn=(1:1:320);
y=zeros(1,320);
y1=y;
ff=ones(1,320);
for n=1:320
    [SinSeries, ~] = FSsine(n,f,L,dx,x);
    y(n)=SinSeries(500);
    y1(n)=SinSeries(bound2-10);
end
fig1 = figure;
plot(nn,ff,'w','LineWidth',2);          hold on; 
plot(nn,y,'c','Linewidth',2)
xlim([1,320]);
ylim([0.8,1.3]);
legend('\color{white}Function f','\color{white}Fourier Series','Color',[0 0 0],'location','best','EdgeColor',[1 1 1])
darkBackground(fig1,[0 0 0],[1 1 1]);
title('\color{white} fuction value at a point');
xlabel('N')
width = 500; height = 200; 
set(gcf, 'Position',  [100, 1000, width, height]);

fig2 = figure;
plot(nn,ff,'w','LineWidth',2);          hold on; 
plot(nn,y1,'c','Linewidth',2)
xlim([1,320]);
ylim([-0.3,1.3]);
legend('\color{white} Function f','\color{white} Fourier Series','Color',[0 0 0],'location','best','EdgeColor',[1 1 1])
darkBackground(fig2,[0 0 0],[1 1 1]);
title('\color{white} x near the jump');
xlabel('N')
width = 500; height = 200;  
set(gcf, 'Position',  [100, 1000, width, height]);



