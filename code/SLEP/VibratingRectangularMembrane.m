%Vibrating Rectangular Membrane
%textbook chapter 7.3 p280
% d^2u/dt^2=c^2(d^2u/dx^2+d^2u/dy^2)
%BC  u(0,y,t)=u(L,y,t)=u(x,0,t)=u(x,H,t)=0
%IC u(x,y,0)=alpha(x,y)   du/dt(x,y,0)=beta(x,y)

%add_my_path
H  = 1;  L = 2;   % domain [a,b]
N  = 99;           % number of x points excluding both ends
dx = L/(N+1); dy = H/(N+1); % "infitasimal" x and y
x  = 0:dx:L; y = 0:dy:H;
c = 0.05;
[X,Y]=meshgrid(x,y);
alpha = @(x,y)x+sin(y);
beta = @(x,y)x+cos(y);

%Obtained three ODEs after the separation of variables
%u(x,y,t)=h(t)f(x)g(y)
% 1. d^2h/dt^2=-lambda*c^2*h
% 2. d^2f/dx^2 = -mu*f
% 3. d^2g/dy^2 = -(lambda-mu)g
%solve 2 and 3 by SLEP solver. 

%solve 2
p1 = x*0+1;
dp1dx = x*0;
q1 = x*0;
r1 = x*0+1;
[v1,e1]    = SLEP(N,x,p1,dp1dx,q1,r1,'Dirichlet');
%solve 3
p2 = y*0+1;
dp2dy = y*0;
q2 = y*0;
r2 = y*0+1;
[v2,e2]    = SLEP(N,y,p2,dp2dy,q2,r2,'Dirichlet');

%IC
alphaxy = alpha(X,Y);
betaxy = beta(X,Y);

%solve A_nm and B_nm (textbook chapter 7.3)
A = zeros(size(Y));
B = A;
for n = 1:N
    for m=1:N
        A(n,m) = sum((sum(alphaxy.*v2(:,m))*dy)/(dot(v2(:,m)'.^2,r2)*dy).*v1(:,n)')*dx/(dot(v1(:,n)'.^2,r1)*dx);
        B(n,m) = sum(sum(betaxy.*v2(:,m).*v1(:,n))*dy)*dx/((dot(v2(:,m)'.^2,r2)*dy)*(dot(v1(:,n)'.^2,r1)*dx))/(c*(e1(n)+e2(m))^0.5);
    end
end

%evaluate u and create animation
index=1;
fig=figure;
K = 50; %number of iteration
tt = 30; %time
for t = 0:0.5:tt
    u=zeros(size(X));
    for n=1:K
        for m=1:K
            [V1,V2]=meshgrid(v1(:,n),v2(:,m));
            u = u+A(n,m)*V1.*V2*cos(c*(e1(n)+e2(m))^0.5*t)+B(n,m)*V1.*V2*sin(c*(e1(n)+e2(m))^0.5*t);
        end
    end
    surf(x,y,u);
    zlim([-5 5]);
    xlabel('x'); ylabel('y'); zlabel('u'); title('Vibrating Rectangular Membrane')
    colormap(parula(100)); 
    set_positionFontsAll;
    drawnow; 
    frame = getframe(fig);
    im{index} = frame2im(frame);
    index=index+1;
    clf;
end
    filename = 'vibratingRecMembrane.gif';
im_to_gif(filename,im,index-1);



