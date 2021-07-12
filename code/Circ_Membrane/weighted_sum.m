function s = weighted_sum(r, theta, t, eig_functions, A_lambda)
u_mni = cellfun(@(F) F(r,theta, t), eig_functions, 'UniformOutput', false);
[M, N] = size(A_lambda);
s = 0;
for m = 1:M
    for n =1:N
        if A_lambda(m,n) == 0 
            continue
        else
            s = s + A_lambda(m,n).*(u_mni{m,n,1} + u_mni{m,n,2});
        end
    end
end




end