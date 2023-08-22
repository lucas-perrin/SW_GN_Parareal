function [F_n,dt] = Solver(W_n,W_border,t,T_end,Riemann)

    dt = T_end - t;

    global dx
    global Nx
        
    % case f = 1
    [F_n(:,1), lambda] = Riemann(W_border(:,1),W_n(:,1));
    dt               = min(dt, dx/(2*lambda));
    
    % case f = Nx+1
    [F_n(:,Nx+1), lambda] = Riemann(W_n(:,Nx),W_border(:,2));
    dt               = min(dt, dx/(2*lambda));
    
    % case f in 2:Nx
    for f = 2:1:Nx
       [F_n(:,f), lambda] = Riemann(W_n(:,f-1),W_n(:,f));
       dt               = min(dt, dx/(2*lambda));
    end
end