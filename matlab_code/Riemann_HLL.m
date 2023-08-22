function [Flux, lambda] = Riemann_HLL(W_L,W_R)
    global g
    global dx
    
    u_l = 0 + (W_L(1) > (dx)^2)*W_L(2)/W_L(1);
    u_r = 0 + (W_R(1) > (dx)^2)*W_R(2)/W_R(1);
    
    lambda_p = max(0,max(u_l + sqrt(g*W_L(1)),u_r + sqrt(g*W_R(1))));
    lambda_m = max(0,-min(u_l - sqrt(g*W_L(1)),u_r - sqrt(g*W_R(1))));
        
    Flux = (lambda_p*Flux_SW(W_L) + lambda_m*Flux_SW(W_R) - lambda_m*lambda_p.*(W_R - W_L))/(lambda_p+ lambda_m);
    
    lambda = max(lambda_m,lambda_p);
end