function [W_j, dt_list, p] = Green_Nagdhi(W_0,T_start,Tj)
    global Nx
    global dx
    
    %W_0 tel que :
    %W_0(1,:) = h
    %W_0(2,:) = hu
    
%     U_0 = Proj(W_0(1,:), [W_0(2,:)./W_0(1,:);zeros(2,Nx)]);
%     W_0 = W_0(1,:).*[ones(1,Nx);U_0];
%     
%     p=0;
%    [W_0(3,:),W_0(4,:)] = get_ws(W_0(1,:),W_0(2,:)./W_0(1,:));
%     W_0(3,:) = W_0(3,:).*W_0(1,:);
%     W_0(4,:) = W_0(4,:).*W_0(1,:);
    
    p = 0;
    J = length(Tj);
    W_j = cell(1,J+1);
    W_j{1} = W_0;

    W_sol = W_0;
    t = T_start;
    
    for i = 1:J
    
        while t < Tj(i)

            dt_end = Tj(i) - t;
            
            [F_sol,dt_r] = Solver_vec(W_sol,t,Tj(i));
            
            dt    = min(dt_r, dt_end);

            % Time integration on evry faces
            W_sol = Euler(W_sol,F_sol,dt);

            % Projection
             U_sol = Proj(W_sol(1,:),zeros(3,Nx) + (W_sol(1,:) > dx^2).*(W_sol(2:4,:)./W_sol(1,:)));
             W_sol = W_sol(1,:).*[ones(1,Nx);U_sol];

            % update 
            t     = t + dt;
            
            p = p+1;
            
        end
        
    W_j{i+1} = W_sol;
    
    end
    %fprintf('p green = %d \n',p)
end