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