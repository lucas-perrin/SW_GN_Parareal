function E = entropy(W)

    global g
    global Nx
    
    h = W(1,:);
    u = zeros(1,Nx) + (W(1,:) > 0).*W(2,:)./W(1,:);
    w = zeros(1,Nx) + (W(1,:) > 0).*W(3,:)./W(1,:);
    s = zeros(1,Nx) + (W(1,:) > 0).*W(4,:)./W(1,:);
    E = (h.^2).*g/2 + (h./2).*(u.^2 + w.^2 + s.^2);
end