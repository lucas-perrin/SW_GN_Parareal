function [W_j, dt_list, p] = Green_Nagdhi(W_0, T_start, Tj)
    global Nx
    global dx
    
    % W_0 where:
    % W_0(1,:) = h (water height)
    % W_0(2,:) = hu (momentum density in x-direction)

    % U_0 = Proj(W_0(1,:), [W_0(2,:)./W_0(1,:);zeros(2,Nx)]);
    % W_0 = W_0(1,:).*[ones(1,Nx);U_0];
        
    % p=0;
    % [W_0(3,:),W_0(4,:)] = get_ws(W_0(1,:),W_0(2,:)./W_0(1,:));
    % W_0(3,:) = W_0(3,:).*W_0(1,:);
    % W_0(4,:) = W_0(4,:).*W_0(1,:);
    
    p = 0; % Initialize the timestep counter
    J = length(Tj); % Number of waypoints
    W_j = cell(1, J + 1); % Cell array to store solutions at different times
    W_j{1} = W_0; % Store the initial condition in the first cell

    W_sol = W_0; % Initialize the current solution to the initial condition
    t = T_start; % Initialize the current time

    dt_list = []; % Initialize the list of time steps
    
    for i = 1:J
        % Continue until the current time reaches the waypoint
        while t < Tj(i)
            dt_end = Tj(i) - t; % Remaining time until the waypoint
            
            [F_sol, dt_r] = Solver_vec(W_sol, t, Tj(i));
            dt = min(dt_r, dt_end); % Choose the smaller of two time steps

            dt_list = [dt_list, dt]; % Store the chosen time step

            % Time integration on every cell face using the Euler method
            W_sol = Euler(W_sol, F_sol, dt);

            % Projection: Update momentum based on water height
            U_sol = Proj(W_sol(1,:), zeros(3, Nx) + (W_sol(1,:) > dx^2).*(W_sol(2:4,:)./W_sol(1,:)));
            W_sol = W_sol(1,:) .* [ones(1, Nx); U_sol];

            t = t + dt; % Update the current time
            p = p + 1; % Increment the timestep counter
        end
        
        W_j{i + 1} = W_sol; % Store the computed solution at the current waypoint
    end
    
end