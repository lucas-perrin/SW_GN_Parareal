function Flux = Flux_SW(W)
    global dx
    global g
    Flux = [0;0;0;0];
    if W(1) > dx^2
        Flux = [W(2); W(2)^2/W(1) + (g*(W(1))^2)/2; W(2)*W(3)/W(1);  W(2)*W(4)/W(1)];
    end      
end