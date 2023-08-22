function [W_j, dt_list, p] = Shallow_Water(W_0, T_start, Tj)
    % Shallow Water model computation
    %
    % Inputs:
    % -> W_0    : Initial condition, of the form (h, hu, hw, hs)
    % -> T_start: Scalar, starting time
    % -> Tj     : Vector of waypoints for the computation
    %
    % Outputs:
    % -> W_j     : Solution at times T_start + Tj
    % -> dt_list : List of all timesteps taken
    % -> p       : Number of timesteps taken (can be obtained using numel(dt_list))

    p = 0; % Initialize the timestep counter
    J = length(Tj); % Number of waypoints
    W_j = cell(1, J + 1); % Cell array to store solutions at different times
    W_j{1} = W_0; % Store the initial condition in the first cell
    
    W_sol = W_0; % Initialize the current solution to the initial condition
    t = T_start; % Initialize the current time
    
    dt_list = []; % Initialize the list of time steps
    
    % Loop through the waypoints
    for i = 1:J
        % Continue until the current time reaches the waypoint
        while t < Tj(i)
            dt_end = Tj(i) - t; % Remaining time until the waypoint

            % Solve Riemann problem on every cell face
            [F_sol, dt_r] = Solver_vec(W_sol, t, Tj(i));
            dt = min(dt_r, dt_end); % Choose the smaller of two time steps

            dt_list = [dt_list, dt]; % Store the chosen time step

            % Time integration on every cell face using the Euler method
            W_sol = Euler(W_sol, F_sol, dt);

            t = t + dt; % Update the current time
            p = p + 1; % Increment the timestep counter
        end
        
        W_j{i + 1} = W_sol; % Store the computed solution at the current waypoint
    end
    
end