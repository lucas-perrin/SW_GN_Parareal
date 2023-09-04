function [F_n,dt] = Solver(W_n,W_border,t,T_end,Riemann)
    % Solves the Riemann Problem with function Riemann. It computes all the flux
    % corresponding to the solution W_n at time t
    % this a depreciated version, new version: Solver_vec
    warning('you are using function Solver now deprciated, new version: Solver_vec')
    % 
    % Inputs:
    % -> W_n: solution at time t
    % -> W_border: _ ?
    % -> t: _ ?
    % -> T_end: _ ?
    % -> Riemann: the Riemann solver, for now available : HLL, HLLC, Rusanov
    %
    % Outputs:
    % -> F_n: flux at time t
    % -> dt: optimal time step for time integration

    % global variables
    global dx
    global Nx

    % time till end
    dt = T_end - t;
    
    % Computing at the two borders, and then in the middle (improved in Solver_vec)
    % case f = 1 (first border of space domain)
    [F_n(:,1), lambda] = Riemann(W_border(:,1),W_n(:,1));
    dt               = min(dt, dx/(2*lambda));
    
    % case f = Nx+1 (second border of space domain)
    [F_n(:,Nx+1), lambda] = Riemann(W_n(:,Nx),W_border(:,2));
    dt               = min(dt, dx/(2*lambda));
    
    % case f in 2:Nx (middle of space domain)
    for f = 2:1:Nx
       [F_n(:,f), lambda] = Riemann(W_n(:,f-1),W_n(:,f));
       dt               = min(dt, dx/(2*lambda));
    end

end