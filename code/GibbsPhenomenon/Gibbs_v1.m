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

L  = .5*pi; 
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
[SinSeries, ~] = FSsine(n,f,L,dx,x);
figure; plotFS_Gibbs(x,f,SinSeries,n,xrange); 

%% Gibbs phenomenon with n=300
n = 300; 
[SinSeries, ~] = FSsine(n,f, L,dx,x);
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


%% animation
index=1;
fig=figure;
for n=20*(1:16)
    [SinSeries, ~] = FSsine(n,f,L,dx,x);
    xrange =  [x(200) x(800)]; 
    plotFS_Gibbs(x,f,SinSeries,n,xrange); set_positionFontsAll;
    drawnow; 
    frame = getframe(fig);
    im{index} = frame2im(frame); %Create cell array of images. (image for each n)
    index=index+1;
    clf;
end

figure; %Show the image cell array
for idx=1:16
    subplot(4,4,idx)
    imshow(im{idx});
end

filename = 'ConvergenceFS.gif';
im_to_gif(filename,im,idx); %Create gif file




