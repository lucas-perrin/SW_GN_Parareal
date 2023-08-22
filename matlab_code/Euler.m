function [W_np1] = Euler(W_n,F_n,dt)
    global Nx
    global dx
    W_np1=W_n;
    % case f = 1
    W_np1(:,1) = W_np1(:,1) + (dt/dx).*F_n(:,1);
    
    % case f = Nx+1
    W_np1(:,Nx) = W_np1(:,Nx) - (dt/dx).*F_n(:,Nx+1);
    
    %case f in 2:Nx
    for f = 2:1:Nx
        W_np1(:,f-1) = W_np1(:,f-1) - (dt/dx).*F_n(:,f);
        W_np1(:,f) = W_np1(:,f) + (dt/dx).*F_n(:,f);
    end
end