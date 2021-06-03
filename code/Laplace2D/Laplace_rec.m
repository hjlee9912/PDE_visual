% This code demonstrates Laplace equation on rectangle 
%with length L and height H

add_my_paths;

%Make rectangle with length L and height H
L  = 1; 
H  = 0.5;
N  = 301;  dx = L/(N-1); dy = H/(N-1); 
NN = N; %for later use
x  = 0:dx:L;               
y  = 0:dy:H;               

[xx,yy] = meshgrid(x,y);   % create meshgrid for 3D plot
%% Boundary Values
f        = @(x) sin(2*pi*x).*exp(x);
K        = 10; % number of terms in FS
type     = 'sine'; % type of FS, 'sine', 'cosine', 'mixed'

%Bdary1: x=0, y=[0:H] with BC=g1
g1       = f;
[fs_coefs1, ~, Bn_all1] = computFS_coef(g1,H,y,K,type);  % compute Fourier coefficiets of g1; 
sinh_part1=(pi*(xx-L)/H);
%Bdary2: x=L, y=[0:H] with BC=g2
g2       = f;
[fs_coefs2, ~, Bn_all2] = computFS_coef(g2,H,y,K,type);  % compute Fourier coefficiets of g2; 
sinh_part2=(pi*xx/H);

%Bdary3: y=0, x=[0:L] with BC=f1
f1       = f;
[fs_coefs3, ~, Bn_all3] = computFS_coef(f1,L,x,K,type);  % compute Fourier coefficiets of f1; 
sinh_part3=(pi*(yy-H)/L);

%Bdary4: y=H, x=[0:L] with BC=f2
f2       = f;
[fs_coefs4, ~, Bn_all4] = computFS_coef(f2,L,x,K,type);  % compute Fourier coefficiets of f2; 
sinh_part4=(pi*yy/L);

%% solution u(x,y): with K FS terms
u0 = zeros(size(xx));
u=u0;
for n = 1:K
    u1  = u0+ sinh(n*sinh_part1) .* sin(n*pi*yy/H) * Bn_all1(n)/sinh(n*pi*(-L)/H);
    u2  = u0+ sinh(n*sinh_part2) .* sin(n*pi*yy/H) * Bn_all2(n)/sinh(n*pi*L/H);
    u3  = u0+ sinh(n*sinh_part3) .* sin(n*pi*xx/L) * Bn_all3(n)/sinh(n*pi*(-H)/L);
    u4  = u0+ sinh(n*sinh_part4) .* sin(n*pi*xx/L) * Bn_all4(n)/sinh(n*pi*H/L);
    u=u+u1+u2+u3+u4;
end
figure; surf(xx,yy,u); xlabel('x'); ylabel('y'); title('Solution of Laplace on rectangle'); 
colormap(parula(100)); shading interp; lighting phong; view([10,32]);
set_positionFontsAll; 
% saveas(gcf,'Laplace_Rectangle.png');
%%  check boundary:  >> plot the boundary condition on the rectangle, which can be used in the future 

xx(3:NN-2,3:NN-3) =   NaN;
yy(3:NN-2,3:NN-2) =   NaN;

fval      = NaN(size(xx));
fval(1,:) = f(x); fval(NN,:) = f(x); fval(2,:) = f(x); fval(NN-1,:) = f(x); 
fval(:,1) = f(y); fval(:,NN) = f(y); fval(:,2) = f(y); fval(:,NN-1) = f(y); 

u_bdry                =   u;
u_bdry(3:NN-2,3:NN-2) =   NaN;

figure; 
subplot(121); surf(xx,yy,fval,'linewidth',2);
         xlabel('x'); ylabel('y'); title('Boundary value expand')
         colormap(parula(100)); shading interp; lighting phong; view([10,32]);

subplot(122); surf(xx,yy,u_bdry,'linewidth',2); 
         xlabel('x'); ylabel('y'); title('solution boundary value')
         colormap(parula(100)); shading interp; lighting phong; view([10,32]);
set_positionFontsAll;
