function plot_water(x,h,u)
    Nx = length(x);
    X = [x(1:Nx-1);x(2:Nx);x(2:Nx);x(1:Nx-1)];
    Y = [h(1:Nx-1);h(2:Nx);zeros(2,Nx-1)];
    C = (u(:,Nx-1)+u(2:Nx))/2;
    %figure(i)
    patch(X,Y,C,'EdgeColor','none')
    colorbar
end