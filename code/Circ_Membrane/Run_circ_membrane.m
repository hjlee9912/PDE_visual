%% Section 7.7: Vibrating circular membrane and Bessel functions 
%{
The vertical displacement u(r, theta, t) of the circular membrane satisfies:
    \partial^2 u/\partial t^2 = c^2 del^2 u.
We assume BC:
    u(a, theta, t) = 0; (fixed boundary)
and IC:
    u(r,theta, 0) = alpha(r, theta);            ( initial displacement)
    \partial u/\partial t = beta(r, theta) = 0. (zero initial velocity)

Method: analytical series solution with Bessel eigenfunctions 
 - 1. Eigenvalues and eigenfunctions
 - 2. Compute coefficients A_lambda from given initial displacemet alpda(r, theta)
 - 3. Solution with a few modes
%}
clear all; close all;
add_my_paths; 

%% parameter settings
a = 1;
c = 1;
M = 6;
N = 6;

%% 1. Eigenvalues and eigenfunctions

% first N zeros of Jm, m = 0, ..., M
z = besselzero((0:M)', N, 1); % z(m, n) is the n^th zero of  J_{m-1}

% eigenvalues lambda_mn [(7.7.40), HaberMan]
sqrt_lambda = z/a; % sqrt_lambda(m,n) is the n^th eigenvalue of Bessel of order (m-1)

% eig_functions(m,n,i) = u_{m-1, n, i}: n-1 nodal circles, 2m nodal radii. [(7.7.45), HaberMan]
% u_{m,n,1}(r, theta, t) = J_m(sqrt(lambda_{m,n})r)*cos(m*theta)*cos(csqrt(lambda_{m,n})t)
% u_{m,n,1}(r, theta, t) = J_m(sqrt(lambda_{m,n})r)*sin(m*theta)*cos(csqrt(lambda_{m,n})t)
eig_functions = cell(M+1, N, 2); 
for m = 0:M
    for n = 1:N
        eig_functions{m+1, n, 1} = @(r, theta, t) besselj(m, r*sqrt_lambda(m+1,n)).*cos(m*theta).*cos(c*sqrt_lambda(m+1,n)*t);
        eig_functions{m+1, n, 2} = @(r, theta, t) besselj(m, r*sqrt_lambda(m+1,n)).*sin(m*theta).*cos(c*sqrt_lambda(m+1,n)*t);
    end
end

% orthogonality of eigen functions 
% enougn to show orthogonality of J_m(sqrt(lambda_{m,n})r) with weight r)
orth_mat = zeros(M+1, N, N);
for m = 0:M
    for n1 = 1:N
        for n2 = 1:N
            fun = @(r) besselj(m, r*sqrt_lambda(m+1,n1)).*besselj(m, r*sqrt_lambda(m+1,n2)).*r;
            orth_mat(m+1,n1, n2) = integral(fun, 0,a);
        end
    end
end
figure;
for i = 1:6
subplot(2,3,i)
imagesc(1:N, 1:N,squeeze(orth_mat(i,:,:))); colorbar;
title(['m = ', num2str(i-1)]);
xlabel('n'); ylabel('n"');
end
% set_positionFontsAll; 
print(['Orth_mat_Jm.pdf'],'-dpdf', '-fillpage');


%% Movie of a few eig_functions
m_seq = 0:2; n_seq = 1:3; i_seq = 1; %m = 0,...,M, n = 1,...,N, i = 1,2
dt = 0.1; L = 20; steps = numel(0:dt:L);
[r_mesh, theta_mesh, time_mesh] = meshgrid(linspace(0,a,20), linspace(-pi,pi,100), 0:dt:L);
x_mesh = r_mesh(:,:,1).*cos(theta_mesh(:,:,1)); y_mesh = r_mesh(:,:,1).*sin(theta_mesh(:,:,1));
vibration = cell(numel(n_seq), numel(m_seq), numel(i_seq));
gif_Filename = [result_path, sprintf('GIF_eig_functions_c_%i_m_%i_%i_n_%i_%i_i_%i_%i_vibration.gif', c,m_seq(1), m_seq(end),n_seq(1), n_seq(end),i_seq(1),i_seq(end))];
if ~exist(gif_Filename,'file') 
    for mm = m_seq
          for nn = n_seq
               for ii = i_seq
                   vibration{mm+1,nn,ii} = eig_functions{mm+1, nn, ii}(r_mesh, theta_mesh, time_mesh);
                end       
          end
    end
    fig = figure;
    for l = 0:steps-1
        count = 1;
        for mm = m_seq
            for nn = n_seq
                for ii = i_seq
                    time = dt * l;
                    subplot(numel(n_seq), numel(m_seq) * numel(i_seq), count)
                    colormap(autumn(100));  view([-40,40]); 
                    surf(x_mesh, y_mesh, vibration{mm+1,nn,ii}(:,:,l+1), 'EdgeAlpha', 0.3);
                    xlabel('x');  ylabel('y'); zlim([-1,1]); 
                    %darkBackground(fig,[0 0 0],[1 1 1]);
                    title(sprintf('u_{%i,%i, %i}', mm, nn, ii), 'FontSize', 18);
                    count = count + 1;
                end       
            end
        end
        sgtitle(['Time t = ', num2str(time,'%.2f')], 'Color', [1 1 1] , 'FontSize', 20);    
        darkBackground(fig,[0 0 0],[1 1 1]);
        drawnow; 
        frame = getframe(fig);
        im{l+1} = frame2im(frame);
    end
    im_to_gif(gif_Filename,im,1:steps);
end


%% 3 Compute coefficients A_lambda for phi)lambda from given initial displacemet alpda(r, theta)
% initial position of eigenfunctions
phi_lambda = cell(M+1, N, 2); % phi_lambda(m,n) = J_{m-1}(sqrt(lambda(m-1,n)r))*[cos(m*theta + sin(m*theta))]
for m = 0:M
    for n = 1:N
        for i = 1:2
        phi_lambda{m+1, n, i} = @(r, theta) eig_functions{m+1, n, i}(r, theta, 0);
        end
    end
end
% initial displacemet alpda(r, theta) = weighted sum of eig_functions
weights = zeros(M+1, N, 2);
weights(1,4,1) = 1;
weights(6,1,1) = 2;
alpha = @(r, theta) weighted_sum(r, theta, 0, eig_functions, weights);
% plot IC
gif_Filename = [result_path,sprintf('GIF_IC_alpha_14_1_61_2.gif', c, a)];
dt = 0.1; L = 20; steps = numel(0:dt:L);
[r_mesh, theta_mesh] = meshgrid(linspace(0,a,20), linspace(-pi,pi,100));
x_mesh = r_mesh.*cos(theta_mesh); y_mesh = r_mesh.*sin(theta_mesh);
initial = alpha(r_mesh, theta_mesh);
fig = figure;
colormap(autumn(100));  view([-40,60]); 
colorbar;
surf(x_mesh, y_mesh, initial , 'EdgeAlpha', 0.3);
bds = max(weights, [],'all');
xlabel('x');  ylabel('y'); zlim([-bds,bds]);  
darkBackground(fig,[0 0 0],[1 1 1]); set(gcf, 'InvertHardCopy', 'off'); 
print('results/IC.png','-dpng');


% change gif_Filename manully for differnt weights
gif_Filename = [result_path,sprintf('GIF_soln_u_c_%i_a_%i_alpha_14_1_61_2.gif', c, a)];

%% 3 Solution u: movie
dt = 0.1; L = 20; steps = numel(0:dt:L);
[r_mesh, theta_mesh, time_mesh] = meshgrid(linspace(0,a,20), linspace(-pi,pi,100), 0:dt:L);
x_mesh = r_mesh(:,:,1).*cos(theta_mesh(:,:,1)); y_mesh = r_mesh(:,:,1).*sin(theta_mesh(:,:,1));

if ~exist(gif_Filename,'file')   
    % A_lambda : u(r, theta, 0) = alpha(r, theta) = \sum A_lambda * phi_lambda. [(7.7.48/49), HaberMan]
    A_lambda = zeros(M+1, N, 2);
    for m = 0:M
        m
        for n= 1:N
            for i = 1:2
                fun1 = @(r, theta) alpha(r, theta).*phi_lambda{m+1, n,i}(r, theta).*r;
                int_alpha_phi = integral2(fun1, 0,a, -pi,pi);
                fun2 = @(r, theta) (phi_lambda{m+1, n,i}(r, theta).^2).*r;
                int_phi2 = integral2(fun2, 0,a, -pi,pi);
                if int_phi2 == 0
                    A_lambda(m+1,n,i) = 0;
                else
                A_lambda(m+1,n,i) = int_alpha_phi/int_phi2;
                end
            end
        end
    end
    
    % Solution u = sum A_lamda(m,n) *[eig_functions(m,n,1) + eig_functions(m,n,1)]
    u = @(r, theta,t) weighted_sum(r, theta, t, eig_functions, A_lambda);
    
%{
figure;
subplot(121);
surf(x_mesh,y_mesh,alpha(r_mesh(:,:,1),theta_mesh(:,:,1)));
title('alpha')
xlabel('x');  ylabel('y'); 
subplot(122);
surf(x_mesh,y_mesh,u(r_mesh(:,:,1),theta_mesh(:,:,1),0));
title('u(.,.,0)')
xlabel('x');  ylabel('y'); 
%}
    % Compare weights in alpha and computed A_lambda
    %figure;
    %imagesc(0:M, 1:N, weights - A_lambda); colorbar; title('true weights - A_{\lambda}')
    vibration = u(r_mesh, theta_mesh, time_mesh);
    fig = figure(14);
    for l = 0:steps-1
        time = dt * l;
        figure(14);
        colormap(autumn(100));  view([-40,60]); 
        colorbar;
        surf(x_mesh, y_mesh, vibration(:,:,l+1), 'EdgeAlpha', 0.3);
        bds = sum(weights, 'all');
        xlabel('x');  ylabel('y'); zlim([-bds,bds]);  
        title({['Time t = ', num2str(time,'%.2f')]},'FontSize', 20); 
        darkBackground(fig,[0 0 0],[1 1 1]);
        drawnow; 
        frame = getframe(fig);
        im{l+1} = frame2im(frame);
    end
    im_to_gif(gif_Filename,im,1:steps);
end

%% Bessel J functions plot
z = 0:0.1:20;
J = zeros(5,201);
for i = 0:4
    J(i+1,:) = besselj(i,z);
end
figure;
plot(z,J, 'Linewidth',2)
grid on
legend('J_0','J_1','J_2','J_3','J_4','Location','Best')
title('Bessel Functions of the First Kind for $m \in [0, 4]$','interpreter','latex')
xlabel('$z$','interpreter','latex')
ylabel('$J_m(z)$','interpreter','latex');
darkBackground(gcf,[0 0 0],[1 1 1]);set(gcf, 'InvertHardCopy', 'off'); 
set_positionFontsAll;
print(['results\Jm.png'],'-dpng');
