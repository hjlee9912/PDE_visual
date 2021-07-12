%{
The vertical displacement u(r, theta, t) of the circular membrane satisfies:
    \partial^2 u/\partial t^2 = c^2 del^2 u.
We assume BC:
    u(a, theta, t) = 0; (fixed boundary)
and IC:
    u(r,theta, 0) = alpha(r, theta); (given initial displacement)
    \partial u/\partial t = beta(r, theta) = 0. (zero initial velocity)
%}

%% parameter settings
a = 1;
c = 1;
M = 6;
N = 6;

%% Eigenvalues and eigenfunctions

% first N zeros of Jm, m = 0, ..., M
z = besselzero((0:M)', N, 1); % z(m, n) is the n^th zero of  J_{m-1}

% eigenvalues lambda_mn
sqrt_lambda = z/a; % sqrt_lambda(m,n) is the n^th eigenvalue of Bessel of order (m-1)

% eig_functions(m,n,i) = u_{m-1, n, i}: n-1 nodal circles, 2m nodal radii
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
set_positionFontsAll;
print(['Orth_mat_Jm.pdf'],'-dpdf', '-fillpage');


%% Movie of a few eig_functions
m_seq = 0:3; n_seq = 1:4; i_seq = 1; %m = 0,...,M, n = 1,...,N, i = 1,2
dt = 0.1; L = 6; steps = numel(0:dt:L);
[r_mesh, theta_mesh, time_mesh] = meshgrid(0:0.05:a, linspace(-pi,pi,100), 0:dt:L);
x_mesh = r_mesh(:,:,1).*cos(theta_mesh(:,:,1)); y_mesh = r_mesh(:,:,1).*sin(theta_mesh(:,:,1));
vibration = cell(numel(n_seq), numel(m_seq), numel(i_seq));
avi_Filename = sprintf('aMovie_eig_functions_c_%i_m_%i_%i_n_%i_%i_i_%i_%i_vibration.avi', c,m_seq(1), m_seq(end),n_seq(1), n_seq(end),i_seq(1),i_seq(end));
if ~exist(avi_Filename,'file') 
    count = 1;
    for mm = m_seq
          for nn = n_seq
               for ii = i_seq
                   vibration{count} = eig_functions{mm+1, nn, ii}(r_mesh, theta_mesh, time_mesh);
                   count = count + 1;
                end       
          end
    end
    v           = VideoWriter(avi_Filename);   % open file for AVI
    v.FrameRate = 10; open(v);
    for l = 0:steps-1
        count = 1;
        for mm = m_seq
            for nn = n_seq
                for ii = i_seq
                    time = dt * l;
                    figure(13);
                    subplot(numel(n_seq), numel(m_seq) * numel(i_seq), count)
                    colormap(autumn(100));  view([-40,40]); 
                    surf(x_mesh, y_mesh, vibration{count}(:,:,l+1));
                    xlabel('x');  ylabel('y'); zlim([-1,1]);  
                    title(sprintf('u_{%i,%i, %i}', mm, nn, ii))
                    count = count + 1;
                end       
            end
        end
        sgtitle(['Time t = ', num2str(time,'%.2f')]);           
        frame = getframe(gcf);  writeVideo(v,frame);
    end
    close(v);
end


%% Compute coefficients A_lambda for phi)lambda from given initial displacemet alpda(r, theta)
% initial position of eigenfunctions
phi_lambda = cell(M+1, N); % phi_lambda(m,n) = J_{m-1}(sqrt(lambda(m-1,n)r))*[cos(m*theta + sin(m*theta))]
for m = 0:M
    for n = 1:N
        phi_lambda{m+1, n} = @(r, theta) eig_functions{m+1, n, 1}(r, theta, 0) + eig_functions{m+1, n, 2}(r, theta, 0);
    end
end
% initial displacemet alpda(r, theta) = weighted sum of eig_functions
weights = zeros(M+1, N);
weights(1,4) = 1;
weights(6,1) = 2;
alpha = @(r, theta) weighted_sum(r, theta, 0, eig_functions, weights);
% change avi_Filename manully for differnt weights
avi_Filename = sprintf('Movie_soln_u_c_%i_a_%i_alpha_14_1_61_2.avi', c, a);

%% Movie solution u
dt = 0.1; L = 6; steps = numel(0:dt:L);
[r_mesh, theta_mesh, time_mesh] = meshgrid(0:0.05:a, linspace(-pi,pi,100), 0:dt:L);
x_mesh = r_mesh(:,:,1).*cos(theta_mesh(:,:,1)); y_mesh = r_mesh(:,:,1).*sin(theta_mesh(:,:,1));
vibration = u(r_mesh, theta_mesh, time_mesh);
if ~exist(avi_Filename,'file')   
    % A_lambda : u(r, theta, 0) = alpha(r, theta) = \sum A_lambda * phi_lambda
    A_lambda = zeros(M+1, N);
    for m = 0:M
        m
        for n= 1:N
            fun1 = @(r, theta) alpha(r, theta).*phi_lambda{m+1, n}(r, theta).*r;
            int_alpha_phi = integral2(fun1, 0,a, -pi,pi);
            fun2 = @(r, theta) (phi_lambda{m+1, n}(r, theta).^2).*r;
            int_phi2 = integral2(fun2, 0,a, -pi,pi);
            A_lambda(m+1,n) = int_alpha_phi/int_phi2;
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

    v           = VideoWriter(avi_Filename);   % open file for AVI
    v.FrameRate = 10; open(v);
    for l = 0:steps-1
        time = dt * l;
        figure(14);
        colormap(autumn(100));  view([-33,71]); 
        colorbar;
        surf(x_mesh, y_mesh, vibration(:,:,l+1));
        bds = sum(weights, 'all');
        xlabel('x');  ylabel('y'); zlim([-bds,bds]);  
        sgtitle({['Time t = ', num2str(time,'%.2f')]});           
        frame = getframe(gcf);  writeVideo(v,frame);
    end
    close(v);
end




