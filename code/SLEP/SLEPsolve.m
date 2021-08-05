% Solves the regular Sturm-Liouville Eigenvalue Problem 
%          [p(x)*y']'+q(x)y = lambda*r(x)*y
% where lambda is an eigenvalue and y(x) is an eigenfunction with domain [a,b]
%{
Boundary conditions: 
   - Dirichlet BC with BC y(a) = 0, y(b) = 0
   - Neumann BC with BC y'(a) = 0, y'(b) = 0
   - mixed BC with BC y(a) = 0, y'(b) = 0
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

[P,dPdx] = fn_p(x); % evaluate p(x) and p'(x)
Q        = fn_q(x); % evaluate Q(x)
R        = fn_r(x); % evaluate R(x), the weight function
BC_type  = 'Dirichlet';      % BC types are 'Dirichlet', 'Neumann', and 'mixed'
[v,e]    = SLEP(N,x,P,dPdx,Q,R,BC_type); % solve SLEP, vector and eigenvaluse (ascending)


fig=figure; 
% v(:,[3]) = - v(:,[3]);     % flip the sign of the eigenfucntion 
plot(x',v(:,1:3),'linewidth',2); 
darkBackground(fig,[0 0 0],[1 1 1]); set_positionFontsAll;
xlabel('\color{white}x')
legend('\color{white}\phi_1', '\color{white}\phi_2','\color{white}\phi_3','location','best')
% print('eigenFns_fowardDiff.pdf','-dpdf');

%% Spectral analysis: 
% 1. plot the eigenvalues 
% 2. show eigenfunctions are orthogonal 
fig2=figure;
subplot(121); plot(e,'y','LineWidth',2); 
title('\lambda'); xlabel('N'); hold on; % plot((1:length(e)).^2,'x'); 
subplot(122); 
plot(x,R,'b','LineWidth',2); hold on; 
plot(x,P,'r','LineWidth',2); 
legend('\color{white}r','\color{white}p','Color',[0 0 0],'EdgeColor',[1 1 1]); 
xlabel('x'); 
darkBackground(fig2,[0 0 0],[1 1 1]); set_positionFontsAll;

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
ylabel(c,'orthogonality')
title('Orthogonality of SLEP eigenfunctions')
