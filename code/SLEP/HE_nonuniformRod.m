%heat flow in nonuniform rod with or without sources at x=[a,b]
% c*rho*du/dt = d/dx(K0*du/dx))+Q
%assume mixed boundary condition u(a,t) = 0 = du/dx(b,t)
%with IC u(x,0) = f(x);
add_my_path;
source = 1; %   0 = no source     1 = put source
a  = 0;  b = 3*pi;   % domain [a,b]
N  = 99;           % number of x points excluding both ends
dx = (b-a)/(N+1);  % "infitasimal" x
x  = a:dx:b;       % domain
f  = @(x) acos(cos(x))+1;
%f  = @(x) x*0+3;  %for f(x) with jump

%set thermal coefficient of the rod
c = @(x) x+6; %specific heat
rho = @(x) sin(x)+2; %mass density
K0 = @(x) 4*x+1; %thermal conductivity
dK0dx = @(x) x*0+4;
alpha = @(x) cos(x); %dependance of Q on u   Q = alpha*u
%% Evaluate the Spatial part
%evaluate SLEP
% -[p(x)*y']'-alpha*y = lambda*r(x)*y      
%y is the eigen function and lambda is eigenvalue
P = K0(x);
dPdx = dK0dx(x);
if source == 1
    Q = alpha(x);
else
    Q = x*0;
end
c_x = c(x);
rho_x = rho(x);
R        = c_x .* rho_x; % evaluate R(x), the weight function
[v,e]    = SLEP(N,x,P,dPdx,Q,R,'mixed'); % solve SLEP, vector and eigenvaluse (ascending)
fig=figure;

plot(x',v(:,1:3),'linewidth',2)
legend('\phi_1', '\phi_2','\phi_3','location','best')
orthogonal = zeros(length(v(1,:)),length(v(1,:)));
for i=1:length(v(1,:))
    for j=i:length(v(1,:))
       orthogonal(i,j) = sum(v(:,i).*v(:,j).*R')*dx;
       orthogonal(j,i) = orthogonal(i,j); 
    end
end
figure;
image(orthogonal,'CDataMapping','scaled')
c=colorbar;
%% Combine the Spatial and time parts
%u = sum(an*v(x)*y_n^(-lambna_n*t)
%first evalulate coefficient 'an'
fx = f(x); %IC
% fx(30:50) = 1; fx(70:90) = 2; %make fx to have jumps
an = zeros(1,N+1); 
for n = 1:N+1
    an(n) = (sum(fx.*R.*v(:,n)')*dx)/(dot(v(:,n)'.^2,R)*dx);
end

%evaluate u and create animation
index=1;
fig=figure;
K = 50; %number of iteration
tt = 30; %time
for t = 0:0.5:tt
    u=0;
    for n = 1:K
        u = u+an(n)*v(:,n)*exp(-e(n)*t);
    end
    plot(x,u,'r','linewidth',2)
    ylim([0 5])
    xlabel('x')
    ylabel('u')
    if source == 0
        title('Heat Flow in a Nonuniform Rod without Sources')
    else
        title('Heat Flow in a Nonuniform Rod with a Source (Q=cos(x)u)')
    end
    darkBackground(fig,[0 0 0],[1 1 1]); set_positionFontsAll;
    legend('\color{white}Heat u','Color',[0 0 0],'EdgeColor',[1 1 1]); 
    xlim([a b])
    set(gca,'XTick',0:pi:3*pi) 
    set(gca,'XTickLabel',{'0','pi','2*pi','3*pi'})
    drawnow; 
    frame = getframe(fig);
    im{index} = frame2im(frame);
    index=index+1;
    clf;
end
if source ==0
    filename = 'nonuniformRod_noSource.gif';
else
    filename = 'nonuniformRod_Source.gif';
end
im_to_gif(filename,im,index-1);




