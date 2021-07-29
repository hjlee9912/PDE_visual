%Vibrating Rectangular Membrane
%textbook chapter 7.3 p280
% d^2u/dt^2=c^2(d^2u/dx^2+d^2u/dy^2)
%BC  u(0,y,t)=u(L,y,t)=u(x,0,t)=u(x,H,t)=0
%IC u(x,y,0)=alpha(x,y)   du/dt(x,y,0)=beta(x,y)

add_my_path
H  = pi;  L = 2*pi;   % domain [a,b]
N  = 99;           % number of x points excluding both ends
dx = L/(N+1); dy = H/(N+1); % "infitasimal" x and y
x  = 0:dx:L; y = 0:dy:H;
c = 2;
alpha = @(x,y)sin(x+y);
beta = @(x,y)sin(x+y);

%Obtained three ODEs after the separation of variables
%u(x,y,t)=h(t)f(x)g(y)
% eq1. d^2h/dt^2=-lambda*c^2*h
% eq2. d^2f/dx^2 = -mu*f
% eq3. d^2g/dy^2 = -(lambda-mu)g
%solve A_nm and B_nm (textbook chapter 7.3)
%evaluate u and create animation
K = 40; %number of iteration
index=1;
fig=figure;
tt = 20; %time
[X,Y]=meshgrid(x,y);
A=zeros(K,K); B=A;
for n=1:K
        for m=1:K
            A_int = @(x,y)alpha(x,y).*sin(m*pi*y/H).*sin(n*pi*x/L);
            B_int = @(x,y)beta(x,y).*sin(m*pi*y/H).*sin(n*pi*x/L);
            A(n,m) = (2/H)*(2/L)*integral2(A_int,0,L,0,H);
            B(n,m) = (2/H)*(2/L)*integral2(B_int,0,L,0,H)/(c*((n*pi/L)^2+(m*pi/H)^2)^0.5);
        end
end


for t = 0:0.2:tt
    u=zeros(size(X));
    for n=1:K
        for m=1:K
            u = u+A(n,m)*sin(m*pi*Y/H).*sin(n*pi*X/L)*cos(c*((n*pi/L)^2+(m*pi/H)^2)^0.5*t)+...
                B(n,m)*sin(m*pi*Y/H).*sin(n*pi*X/L)*sin(c*((n*pi/L)^2+(m*pi/H)^2)^0.5*t);
        end
    end
    surf(x,y,u);
    zlim([-2 2]);
    xlabel('x'); ylabel('y'); zlabel('u'); title('Vibrating Rectangular Membrane')
    colormap(parula(50)); 
    caxis([-2 2])
    s.EdgeColor = 'none';
    colorbar;
    set_positionFontsAll;
    drawnow; 
    frame = getframe(fig);
    im{index} = frame2im(frame);
    index=index+1;
    clf;
end
    filename = 'vibratingRecMembrane_analytic.gif';
im_to_gif(filename,im,index-1);



