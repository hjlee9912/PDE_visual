% heat flow in nonuniform rod without sources at x=[a,b]
%          c*rho*du/dt = d/dx(K0*du/dx))
% BC: mixed boundary condition u(a,t) = 0 = du/dx(b,t)
% IC: u(x,0) = f(x);


add_my_path;
a  = 0;  b = pi;   % domain [a,b]
N  = 99;           % number of x points excluding both ends
dx = (b-a)/(N+1);  % "infitasimal" x
x  = a:dx:b;       % domain
f  = @(x) acos(cos(x))+1;
%set thermal coefficient of the rod
c   = @(x) x+6;    % specific heat
rho = @(x) sin(x)+2;  % mass density     % % why it does not work for 1+1.5*cos(x)? --> rho should be positive
K0  = @(x) 4*x+1; % thermal conductivity
dK0dx = @(x) x*0+4;

%% Evaluate the Spatial part
%evaluate SLEP
%%[p(x)*y']'+ = lambda*r(x)*y      
%y is the eigen function and lambda is eigenvalue
P = K0(x);
dPdx = dK0dx(x);
Q = x*0;
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
% u = sum(an*v(x)*y_n^(-lambna_n*t)
%first evalulate coefficient 'an'
fx = f(x); %IC
an = zeros(1,N); 
for n = 1:N
    an(n) = (sum(fx.*R.*v(:,n))*dx)/(dot(v(:,n)'.^2,R)*dx);
end

%evaluate u and create animation
index=1;
fig=figure;
K = 50; %number of iteration
tt = 30; %time
for t = 0:0.5:tt
    u = 0; %heat
    for n = 1:K
        u = u+an(n)*v(:,n)*exp(-e(n)*t);
    end
    plot(x,u,'linewidth',2)
    ylim([0 5])
    xlabel('x')
    ylabel('u')
    title('Heat Flow in a Nonuniform Rod without Sources')
    darkBackground(fig,[0 0 0],[1 1 1]); set_positionFontsAll;
    legend('\color{white}u(heat)','Color',[0 0 0],'EdgeColor',[1 1 1]); 
    drawnow; 
    frame = getframe(fig);
    im{index} = frame2im(frame);
    index=index+1;
    clf;
end
filename = 'nonuniformRod.gif';
im_to_gif(filename,im,index-1);


