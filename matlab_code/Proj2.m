function W_out = Proj2(W_in)
    global Nx
    global dx
    
    h = W_in(1,:);
    U = zeros(3,Nx) + (h > dx^2).*(W_in(2:4,:)./h);
    Pi_h = Proj(h,U);
    W_out = h.*[ones(1,Nx);Pi_h];
    
end