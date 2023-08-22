function [w, s] = get_ws(h,u)
    global dx
    
    w = - h.*([u(2:end),-u(end)] - [-u(1),u(1:end-1)])./(4*dx);
    s = w./sqrt(3);
    
end