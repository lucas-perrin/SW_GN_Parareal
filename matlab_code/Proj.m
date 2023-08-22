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