%vibration on a nonuniform string at x=[a,b]
% rho*du"/d"t = T0*du"/d"x
%assume Dirichlet boundary condition u(a,t) = 0 = u(b,t)
%with IC u(x,0) = f(x); and du/dt(x,0) = g(x);

add_my_path;
a  = 0;  b = pi;   % domain [a,b]
N  = 99;           % number of x points excluding both ends
dx = (b-a)/(N+1);  % "infitasimal" x
x  = a:dx:b;       % domain
f  = @(x) sin(x);
g  = @(x) cos(x);

rho = @(x) sin(x)+3; %density
T0  = 5; %tension (constant)

%set variables for SLEP.m
P = x*0+T0;
dPdx = x*0;
Q = x*0;
R = rho(x);
%solve eigenvalue problem: 
%T0d"+lambda*rho(x)*y=0 where y is eigenfunction
[v,e]    = SLEP(N,x,P,dPdx,Q,R,'Dirichlet'); % solve SLEP, vector and eigenvaluse (ascending)
fig=figure;
plot(x',v(:,1:3),'linewidth',2)%show first three eigenfunction
legend('\phi_1', '\phi_2','\phi_3','location','best')
orthogonal = zeros(length(v(1,:)),length(v(1,:))); %show orthogonality
for i=1:length(v(1,:))
    for j=i:length(v(1,:))
       orthogonal(i,j) = sum(v(:,i).*v(:,j).*R')*dx;
       orthogonal(j,i) = orthogonal(i,j); 
    end
end
figure;
image(orthogonal,'CDataMapping','scaled')
c=colorbar;

%% now combine spatial and time parts
%solution is 
%u = sum(an*sin(lambda^(0.5))*t*y)+sum(bn*cos(lambda^(0.5))*t*y)
fx = f(x); gx = g(x);
an = zeros(1,N); bn = an; %solve coefficients  using orthogonality
for n=1:N
    an(n)=dot(gx.*v(:,n)',R)*dx/(e(n)^(0.5)*(dot(v(:,n)'.^2,R)*dx));
    bn(n)=dot(fx.*v(:,n)',R)*dx/(dot(v(:,n)'.^2,R)*dx);
end

%evaluate u and create animation
index=1;
fig=figure;
K = 40; %number of iteration
tt = 50; %time
for t = 0:0.3:tt
    u = 0; %displacement
    for n = 1:K
        u = u+an(n)*sin(e(n)^(0.5)*t)*v(:,n)+bn(n)*cos(e(n)^(0.5)*t)*v(:,n);
    end
    plot(x,u,'y','linewidth',2)
    title('Vibrations of a Nonuniform String')
    darkBackground(fig,[0 0 0],[1 1 1]); set_positionFontsAll;
    legend('\color{white}u(displacement)','Color',[0 0 0],'EdgeColor',[1 1 1]); 
    ylim([-20 20])
    xlabel('x')
    ylabel('u')
    drawnow; 
    frame = getframe(fig);
    im{index} = frame2im(frame);
    index=index+1;
    clf;
end
filename = 'nonuniformString.gif';
im_to_gif(filename,im,index-1);
