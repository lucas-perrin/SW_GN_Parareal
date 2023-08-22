function [Flux, lambda] = Riemann_Rusanov(W_L,W_R)
    global g
    global dx
    
    u_l = 0 + (W_L(1) > (dx)^2)*W_L(2)/W_L(1);
    u_r = 0 + (W_R(1) > (dx)^2)*W_R(2)/W_R(1);
    
    lambda = max(abs(u_l) + sqrt(g*W_L(1)),abs(u_r) + sqrt(g*W_R(1)));
        
    Flux = (Flux_SW(W_L) + Flux_SW(W_R) - lambda.*(W_R - W_L))/2;
end