function [v,e]=SLEP(N,x,P,dPdx,Q,R,type)
%This function solves Regular Sturm-Liouville Eigenvalue Problem 
%[p(x)*y']'+q(x)y = lambda*r(x)*y
%where y(x) is eigenfunction and lambda is eigenvalue
%uses finite difference method and matrix linear operations

%v is (N+1)x(N+1) matrix where the 1st and last rows satisfy BC 
    %and each column is a eigen vector corresponding to certain eigenvalue
%e is a column vector containing eigenvalues


    
dx = x(2)-x(1); 
switch type
    case 'Dirichlet'
        PP = zeros(N,N); 
        dPPdx = PP;
        %make PP as NxN matrix where each column is p(x)
        %make dPPdx as NxN matrix where each column is dp/dx
        for i=1:N
            PP(:,i) = P(2:end-1)';
            dPPdx(:,i) = dPdx(2:end-1)';
        end
        RR = diag(R((2:end-1))); 
        QQ = diag(Q((2:end-1)));
        %represent [p(x)*y']'+q(x)y = lambda*r(x)*y as
        % Ay+By+Cy = lambda*RR*y 
        %A contains y" part
        A = diag(2*ones(1,N))-diag(ones(1,N-1),-1)-diag(ones(1,N-1),1);
        A = dx^(-2)*A.*PP;
        %B contains y' part
        B = diag(-1*ones(1,N))+diag(ones(1,N-1),1);
        B = dx^(-1)*B.*dPPdx;
        %C is q(x)y part
        C = QQ;
        D = A+B+C; %linear operation
        %solve eigenvalue problem D*y = lambda*y
        [v,e] = eig(D,RR);
        v = [zeros(1,N);v;zeros(1,N)];
    case 'Neumann'
        PP = zeros(N+1,N+1); 
        dPPdx = PP;
        %make PP as (N+1)x(N+1) matrix where each column is p(x)
        %make dPPdx as (N+1)x(N+1) matrix where each column is dp/dx
        for i=1:N+1
            PP(:,i) = P(2:end)';
            dPPdx(:,i) = dPdx(2:end)';
        end
        RR = diag(R(2:end)); 
        QQ = diag(Q(2:end));
        %represent [p(x)*y']'+q(x)y = lambda*r(x)*y as
        % Ay+By+Cy = lambda*RR*y 
        %A contains y" part
        A = diag(2*ones(1,N+1))-diag(ones(1,N),-1)-diag(ones(1,N),1);
        A(N+1,N+1)=1;
        A = dx^(-2)*A.*PP;
        %B contains y' part
        B = diag(-1*ones(1,N+1))+diag(ones(1,N),1);
        B(N+1,N+1) = 0; 
        B = dx^(-1)*B.*dPPdx;
        %C is q(x)y part
        C = QQ;
        D = A+B+C; %linear operation
        %solve eigenvalue problem D*y = lambda*y
        [v,e] = eig(D,RR);
        v = [zeros(1,N+1);v];
        
        case 'Robin'
            %stil working on it;
end

e = diag(e);
[e,ind]  = sort(e,'ascend'); v = v(:,ind); % (ascending)
end
