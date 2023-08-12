function [W, time] = Parareal_v2(U0,Tj,Fine,Coarse,Riemann)
% PARAREAL implementation of the parareal algorithm

N = length(Tj) - 1;

W = cell(N+1,N+1);
Go = cell(1,N+1);
Gn = cell(1,N+1);
F = cell(1,N+1);
W{1,1} = U0;

t1 = 0;
t2 = 0;
t3 = 0;
time = zeros(1,N);

% initial guess with G
for n=1:N
    tic
    Go{n+1} = Coarse(W{1,n},Tj(n),Tj(n+1));
    t1 = t1 + toc;
    W{1,n+1} = Go{n+1}{end};
end
% parareal iteration

for k=1:N
    
    fprintf('k = %d \n',k)
    
    for n=1:N
        tic
        Wkn = W{k,n};
        Ukn = Proj(Wkn(1,:),zeros(3,Nx) + (Wkn(1,:) > dx^2).*(Wkn(2:4,:)./Wkn(1,:)));
        Wkn = Wkn(1,:).*[ones(1,Nx);Ukn];
        F{n+1} = Fine(Wkn,Tj(n),Tj(n+1));     % parallel with F
        t2 = t2 + toc;
    end
    
    for j=1:k
        W{k+1,j} = W{k,j};
    end
    
    for n=1:N
        tic
        Gn{n+1} = Coarse(W{k+1,n},Tj(n),Tj(n+1));       % sequential with G
        W{k+1,n+1} = F{n+1}{end} + Gn{n+1}{end} - Go{n+1}{end};  % parareal update
        t3 = t3 + toc;
    end
    
    Go = Gn;
    
    time(k) = t1 + t2/N + t3;
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function U = Proj(hs,Us)

    % WARNING : spdiags() call sparse() so maybe better to code everything
    % with sparse or, with only 1 call of spdiags()

    global Nx
    global dx
    
    kappa = ([hs(1),hs,hs(end)].^3)./3;
    A = sparse(1:Nx,1:Nx,(kappa(1:end-2) + kappa(3:end)),Nx,Nx) +...
    sparse(1:Nx-2,3:Nx,-kappa(3:end-2),Nx,Nx) +...
    sparse(3:Nx,1:Nx-2,-kappa(3:end-2),Nx,Nx) +...
    sparse([1, 2, Nx-1, Nx],[2, 1, Nx, Nx-1],[kappa(1), kappa(2), kappa(end-1), kappa(end)],Nx,Nx);
    A = A./(4*dx^2);
    A = sparse(1:Nx,1:Nx,hs,Nx,Nx) + A;

    us = Us(1,:); ws = Us(2,:); ss = Us(3,:);
    
%     h_kp1 = [hs(2:end),hs(end)];
%     h_km1 = [hs(1),hs(1:end-1)];
%     w_kp1 = [ws(2:end),ws(end)];
%     w_km1 = [ws(1),ws(1:end-1)];
%     s_kp1 = [ss(2:end),ss(end)];
%     s_km1 = [ss(1),ss(1:end-1)];
    
    beta = hs.*us + ([hs(2:end),hs(end)].^2.*([ws(2:end),ws(end)] + [ss(2:end),ss(end)]./sqrt(3)) - [hs(1),hs(1:end-1)].^2.*([ws(1),ws(1:end-1)] + [ss(1),ss(1:end-1)]./sqrt(3)))./(4*dx);
    
    u = mldivide(A,beta')';
    
    [w, s] = get_ws(hs,u);
    
    U = [u;w;s];
end