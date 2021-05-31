function [FS] = LaplaceSol_rec(K,f,f_name,L,H,x,y,dx,dy)
% evaluate the Fourier seires with n modes
% the coefficients are approximated by Riemann sum (can do better by integration functions in MATLAB)
FS=0;
% An_all    = zeros(a(1),K);
for n=1:K %How many "n" vales? Sine seires starts with n=1
    
    if f_name==1
        sinh_mesh=(pi*(x-L)/H);
        sin_mesh=(pi*y/H);
        coeff=2/(H*sinh(n*pi*(-L)/H));
        dvar=dy;
        eigenf     = sin(n*sin_mesh); %Eigen function for the sine term
        h          = sinh(n*sinh_mesh); %sinh term
        An         = dot(f,eigenf,1).*dvar*coeff; %Coefficient 
        FS         = FS+An.*eigenf.*h; %Solution
%         An_all(n)  = An; 
    
    elseif f_name==2
        sinh_mesh=(pi*x/H);
        sin_mesh=(pi*y/H);
        coeff=2/(H*sinh(n*pi*L/H));
        dvar=dy;
        eigenf     = sin(n*sin_mesh); 
        h          = sinh(n*sinh_mesh);
        An         = dot(f,eigenf,1).*dvar*coeff; 
        FS         = FS+An.*eigenf.*h; 
%         An_all(n)  = An; 
        
    elseif f_name==3
        sinh_mesh=(pi*(y-H)/L);
        sin_mesh=(pi*x/L);
        coeff=2/(L*sinh(n*pi*(-H)/L));
        dvar=dx;
        eigenf     = sin(n*sin_mesh); 
        h          = sinh(n*sinh_mesh);
        An         = dot(f,eigenf,2).*dvar*coeff;
        FS         = FS+An.*eigenf.*h; 
%         An_all(n)  = An; 
        
    elseif f_name==4
        sinh_mesh=(pi*y/L);
        sin_mesh=(pi*x/L);
        coeff=2/(L*sinh(n*pi*H/L));
        dvar=dx;
        eigenf     = sin(n*sin_mesh); 
        h          = sinh(n*sinh_mesh);
        An         = dot(f,eigenf,2).*dvar*coeff; 
        FS         = FS+An.*eigenf.*h; 
%         An_all(n)  = An; 
    end


end

end