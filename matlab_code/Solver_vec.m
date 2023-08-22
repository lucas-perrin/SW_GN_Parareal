function [F_faces,dt] = Solver_vec(W,t,T_end)

    global dx
    global Nx
    global g

    dt = T_end - t;
    
    % compute w on border
    W_ = [[1;-1;1;1].*W(:,1), W, [1;-1;1;1].*W(:,Nx)];

    % compute flux on each cell
    F_ = [W_(2,:); (W_(2,:).^2)./W_(1,:) + (g.*(W_(1,:)).^2)./2; W_(2,:).*W_(3,:)./W_(1,:);  W_(2,:).*W_(4,:)./W_(1,:)];

    % compute lambda on faces
    H_ = W_(1,:);
    U_ = W_(2,:)./W_(1,:);
    lambda_faces = max([abs(U_(1:Nx+1)) + sqrt(g.*H_(1:Nx+1)) ; abs(U_(2:Nx+2)) + sqrt(g.*H_(2:Nx+2))], [], 1);

    % compute flux on faces
    F_faces = (F_(:,1:Nx+1) + F_(:,2:Nx+2) - lambda_faces.*(W_(:,2:Nx+2) - W_(:,1:Nx+1)))/2;

    % compute dt
    dt = min([dt, dx./(2*lambda_faces)]);
    
end