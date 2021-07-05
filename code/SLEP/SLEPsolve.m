% Solves the regular Sturm-Liouville Eigenvalue Problem 
%          [p(x)*y']'+q(x)y = lambda*r(x)*y
% where lambda is an eigenvalue and y(x) is an eigenfunction with domain [a,b]
%{
Boundary conditions: 
   - Dirichlet BC with BC y(a) = 0, y(b) = 0
   - Neumann: 
   - Robin: mixed
Regular SLEP
  - p,q,r are continuous on [a,b]; p>0, r>0 --- including boundary points 
     (use p(x) = sin(x), r(x) = abs(sin(x))+1 to show why endpoints are important >> theory refer to book XXX)
  - BC: 
==>> uses function 'SLEP' and outputs 
  - a plot of: The first three eigenfunctions vs x
  - spectral analysis 
%}  
% By Hojun Lee and Fei Lu 
% Jul/4/2021
add_my_path;

a  = 0;  b = pi;   % domain [a,b]
N  = 99;           % number of x points excluding both ends
dx = (b-a)/(N+1);  % "infitasimal" x
x  = a:dx:b;       % domain
x_in = x(2:end-1); % x values excluding both ends

[P,dPdx] = fn_p(x_in); % evaluate p(x) and p'(x)
Q        = fn_q(x_in); % evaluate Q(x)
R        = fn_r(x_in); % evaluate R(x), the weight function
[v,e]    = SLEP(N,x_in,P,dPdx,Q,R); % solve SLEP, vector and eigenvaluse (ascending)
% [e,ind]  = sort(e); v = v(:,ind); (order? ascending/descending? TBD)


fig=figure;
plot(x',v(:,1:3),'linewidth',2)
darkBackground(fig,[0 0 0],[1 1 1]); set_positionFontsAll;
xlabel('\color{white}x')
legend('\color{white}\phi_1', '\color{white}\phi_2','\color{white}\phi_3')


%% Spectral analysis: 
% 1. plot the eigenvalues (order? ascending/descending?
% 2. show eigenfunctions are orthogonal
figure;
subplot(121); plot(e); hold on; plot((1:length(e)).^2,'x'); %figure; semilogy(e);
subplot(122); plot(R); hold on; plot(P); legend('r','p'); xlabel('x'); 







