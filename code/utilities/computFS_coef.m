function [fs_coefs, An_all, Bn_all] = computFS_coef(f,L,x_mesh,K,type)
% compute Fourier coefficients of 1D function f with mesh x. 
% Input
%      f    - a funciton handle 
%      x    - space mesh for integration: interval [0,L], [-L, L]; 
%      K    - number of coefficients to be computed; if not provides, the number of Fourier coefficients = length(x) 
%      type - sine cosine or mix
% Output 
%      Fourier coeffcients of f;   NOTE: An(1) -->> a_0; An(2) = a_1;  

fval    = f(x_mesh); 
dx      = x_mesh(2)-x_mesh(1);  
x_mesh  = x_mesh*(pi/L); 

if ~exist('K','var');     K    = length(x_mesh);  end

An_all = zeros(1,K+1);  Bn_all = zeros(1,K);
switch type
    case 'sine'    % evaluate Bn
        for n = 1:K
            eigenfn = sin(n*x_mesh);    % eigen-function
            Bn_all(n)  = dot(fval,eigenfn)*dx*2/L;   
        end
    case 'cosine' % evaluate An
        for n = 1:K
            eigenfn = cos(n*x_mesh);    % eigen-function
            An_all(n)  = dot(fval,eigenfn)*dx*2/L;   
        end
    case 'mixed'  % evaluate An & Bn
        An_all(1) = mean(fval); 
        for n=1:K
            sin_n       = sin(n*x_mesh);    % eigen-function
            Bn_all(n)   = dot(fval,sin_n)*dx*2/L; 
            cos_n       = cos(n*x_mesh);    % eigen-function
            An_all(n+1) = dot(fval,cos_n)*dx*2/L; 
        end
end
fs_coefs.An_cosine = An_all;  
fs_coefs.Bn_sine   = Bn_all; 
    
 