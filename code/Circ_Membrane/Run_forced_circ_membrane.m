
%% Section 8.5: Add Force 
%{
The vertical displacement u(r, theta, t) of the circular membrane satisfies:
    \partial^2 u/\partial t^2 = c^2 del^2 u + Q(r, theta, t) <-- external force.
We assume BC:
    u(a, theta, t) = 0; (fixed boundary)
and IC:
    u(r,theta, 0) = alpha(r, theta);            (initial displacement)
    \partial u/\partial t = beta(r, theta) = 0. (initial velocity)

Method: analytical series solution with Bessel eigenfunctions 
 - 1. Homogennous eigenvalues and eigenfunctions
 - 2. Compute time-dependent coefficients A_i(t) 
 - 3. Solution u(r, theta, t) = sum A_i(t) eig_fun_i(r, theta)
%}
clear all; close all;
add_my_paths; 

%% parameter settings
a = 1; % a: radius
c = 1;
M = 6;
N = 6;

%% 1. Homogenous eigenvalues and eigenfunctions

% first N zeros of Jm, m = 0, ..., M
z = besselzero((0:M)', N, 1); % z(m, n) is the n^th zero of  J_{m-1}

% eigenvalues lambda_mn [(7.7.40), HaberMan]
sqrt_lambda = z/a; % sqrt_lambda(m,n) is the sqrt of the n^th eigenvalue of Bessel of order (m-1)

% eig_functions(m,n,i) = u_{m-1, n, i}: n-1 nodal circles, 2m nodal radii. [(8.5.5), HaberMan]
% u_{m,n,1}(r, theta, t) = J_m(sqrt(lambda_{m,n})r)*cos(m*theta)
% u_{m,n,2}(r, theta, t) = J_m(sqrt(lambda_{m,n})r)*sin(m*theta)
eig_functions = cell(M+1, N, 2); 
for m = 0:M
    for n = 1:N
        eig_functions{m+1, n, 1} = @(r, theta, t) besselj(m, r*sqrt_lambda(m+1,n)).*cos(m*theta);
        eig_functions{m+1, n, 2} = @(r, theta, t) besselj(m, r*sqrt_lambda(m+1,n)).*sin(m*theta);
    end
end

%% 2. [(8.5.19), HaberMan] A_i(t) = c1*cos(csqrt(lambda_i)t) + c2*sin(csqrt(lambda_i)t) + integral_0^t{qi(s)sin(c*sqrt(lambda_i)(t-s))/c/sqrt(lambda_i) ds}
%% set alpha and Q
% initial displacemet alpda(r, theta) = weighted sum of eig_functions
alpha = @(r, theta) eig_functions{1,4}(r, theta) + 2*eig_functions{6,1}(r, theta);
avi_Filename = [result_path, sprintf('Movie_soln_u_c_%i_a_%i_alpha_14_1_61_2.avi', c, a)];
% force Q(r, theta, t) = weighted sum of eig_functions
force_type = 'resonance';
avi_Filename = strrep(avi_Filename, '_u_', ['_u_Force_',force_type,'_'] );
q_i = cell(M+1, N, 1);
switch force_type
    case 'increase' 
        m = 1; n = 1;
        Q = @(r, theta, t) t .* eig_functions{m,n}(r, theta, t);
        q_i{m,n} = @(t) t;
    case 'periodic'
        m = 1; n = 1;
        Q = @(r, theta, t) cos(t) .* eig_functions{1,1}(r, theta, t);
        q_i{m,n} =@(t) cos(t);
    case 'resonance'
        m = 1; n = 1; omega = c*sqrt_lambda(m,n);
        Q = @(r, theta, t) cos(omega*t) .* eig_functions{m,n}(r, theta, t);
        q_i{m,n} =@(t) cos(omega*t);
end
q_i(cellfun(@isempty,q_i)) = {@(t) 0};

dt = 0.05; L = 6; steps = numel(0:dt:L);
integral_term = zeros(M+1, N, steps);
for m = 0:M
    for n = 1:N
        for h = 1:steps
            t = dt*(h-1);
            mesh_size = dt/20;
            t_mesh = 0:mesh_size:t;
            integral_term(m+1, n, h) = sum(mesh_size*q_i{m+1, n}(t_mesh) .* sin(c*sqrt_lambda(m+1, n)*(t-t_mesh)/c/sqrt_lambda(m+1, n)));
        end
    end
end
c1_mat = zeros(M+1, N);
for m = 0:M
    for n= 1:N
        fun1 = @(r, theta) alpha(r, theta).*eig_functions{m+1, n}(r, theta).*r;
        int_alpha_phi = integral2(fun1, 0,a, -pi,pi);
        fun2 = @(r, theta) (eig_functions{m+1, n}(r, theta).^2).*r;
        int_phi2 = integral2(fun2, 0,a, -pi,pi);
        c1_mat(m+1,n) = int_alpha_phi/int_phi2;
    end
end
% c2 = 0 because initial velocity is set to 0 here

A_t = zeros(M+1, N, steps);
for m = 0:M
    for n = 1:N
        for h = 1:steps
            t = dt*(h-1);
            A_t(m+1, n, h) = integral_term(m+1, n, h) + c1_mat(m+1, n)*cos(c*sqrt_lambda(m+1, n)*t);
        end
    end
end


%% 3. Solution u = sum A_i(t) eig_fun_i(r, theta) 
u = cell(steps); % U{h}(r, theta) is solution at time (h-1)*dt 
for h = 1:steps
    t = dt*(h-1);
    weights =  A_t(:,:,h);
    u{h} =@(r, theta) weighted_sum(r, theta, t, eig_functions, weights);
end

% plot solution u
[r_mesh, theta_mesh] = meshgrid(0:0.05:a, linspace(-pi,pi,100));
x_mesh = r_mesh(:,:,1).*cos(theta_mesh(:,:,1)); y_mesh = r_mesh(:,:,1).*sin(theta_mesh(:,:,1));
if ~exist(avi_Filename,'file')      
    v           = VideoWriter(avi_Filename);   % open file for AVI
    v.FrameRate = 1/dt; open(v);
    vibration = zeros([size(r_mesh),steps]);
    for l = 0:steps-1    
        vibration(:,:,l+1) = u{l+1}(r_mesh, theta_mesh);            
    end
    bds = max(vibration, [],'all');
    for l = 0:steps-1    
        time = dt * l;
        figure(14);
        colormap(autumn(100));  view([-33,71]); 
        colorbar;
        surf(x_mesh, y_mesh, vibration(:,:,l+1));
        xlabel('x');  ylabel('y'); zlim([-bds,bds]);  
        sgtitle({['Time t = ', num2str(time,'%.2f')]});           
        frame = getframe(gcf);  writeVideo(v,frame);
    end
    close(v);
end


