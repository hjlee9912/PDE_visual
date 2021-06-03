% This code demonstrates Laplace equation on circular disk with radius a=1
% BC = f(a,theta) = \sum_n a_n \cos(n \theta) + b_n \sin(n \theta)
% The solution is 
%       u(r,theta) = \sum_n a_n r^n cos(n*theta) + b_n r^b sin(n*theta)
% We present it cartian coordinates at the end. 

add_my_paths;

N     = 301;  dt = 2*pi/(N-1);  % size of the "infinitasimal" theta for integration
theta = -pi:dt:pi; 

a     = 2; %radius
dr    = a/(N-1);  % size of the "infinitasimal" r for integration
r     = 0:dr:a; 

[tt,rr] = meshgrid(theta,r);   % create meshgrid for 3D plot
[xx,yy] = pol2cart(tt,rr); % convert polar to cartesian

%% Boundary value:  - should not use the whole mesh for f --- only use the boudnary 
% boundary value in polor coordinates: f(theta)
f        = @(x) sin(x).*exp(x);
K        = 10; % number of terms in FS
type     = 'mixed'; % type of FS, 'sine', 'cosine', 'mixed'
[fs_coefs, An_all, Bn_all] = computFS_coef(f,pi,theta,K,type);  % compute Fourier coefficiets of f; 

%% solution u(r,theta): with K FS terms
u = zeros(size(tt)); 
u = u + An_all(1); 
for n = 1:K
    u  = u+ rr.^n .* ( An_all(n+1)* cos(n*theta) + Bn_all(n)*sin(n*theta));
end
figure; surf(xx,yy,u); xlabel('x'); ylabel('y'); title('Solution of Laplace on disk'); 
colormap(parula(100)); shading interp; lighting phong; view([10,32]);
set_positionFontsAll; 
saveas(gcf,'Laplace_Disk.png');


%%  check boundary:  >> plot ring, which can be used in the future 
[tt,rr] = meshgrid(theta,r(end-2:end));   % create meshgrid of boundary
[xx,yy] = pol2cart(tt,rr);                % convert polar to cartesian
fval    = f(tt); 
u_bdry  = u(end-2:end,:); 
figure; 
subplot(121); surf(xx,yy,fval,'linewidth',2);
         xlabel('x'); ylabel('y'); title('Boundary value expand')
         colormap(parula(100)); shading interp; lighting phong; view([10,32]);

subplot(122); surf(xx,yy,u_bdry,'linewidth',2); 
         xlabel('x'); ylabel('y'); title('solution boundary value')
         colormap(parula(100)); shading interp; lighting phong; view([10,32]);
set_positionFontsAll;
saveas(gcf,'Laplace_Disk_bdry.png');  
