%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   Work with M. Parisot : Parareal, SW, GN   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is version 2 : trying working parareal with good energy

% check this if you want to compute
do_computations = 0;
if do_computations
    clearvars -except do_computations parallel
end
close all

% define global variables

% global T_end
global g
global dx
global Nx
global L
global x
global epsilon

% for the plots

Fontsize = 18;
Fontsize_label = 18;
Fontsize_axes = 18;
Linesize = 2;
Marksize = 9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===== Parameters to tune =====

L  = 1;    % length of intervals
Nx = 5.e3; % Number of cells % // good : 1.e3 / 5.e3
%Nx = 5.e2;
%Nx = 1.e3;

epsilon  = 5.e-3; % // good : 5.e-3 / 10 1.e3

initial_condition = 'gaussian'; % choices : 'gaussian' or 'riemann'
gaussian = 1;     % initial condition (Left Gaussian (1) or Riemann problem (0))


% ===== Other Parameters =====

dx = L/Nx;    % space step
g  = 9.18;    % gravity constant

% space (for plots)
x  = linspace(0+dx/2,1-dx/2,Nx);


if do_computations
% ===== Initial condition =====

riemann = 1 - gaussian;

W_0 = zeros(4,Nx);

if initial_condition == 'gaussian'
    T_end = .6/sqrt(g*epsilon);
    W_0(1,:) = epsilon.*(1 + exp(-(x/1e-1).^2));
    W_0(2,:) = zeros(1,Nx);
else
    T_end = 0.4/sqrt(g*epsilon);
    W_0(1,1:floor(Nx/2)) = epsilon*1;
    W_0(2,1:floor(Nx/2)) = 0;
    W_0(1,floor(Nx/2)+1:end) = epsilon*0.5;
    W_0(2,floor(Nx/2)+1:end) = 0;
end

U_0 = Proj(W_0(1,:),zeros(3,Nx) + (W_0(1,:) > dx^2).*(W_0(2:4,:)./W_0(1,:)));
W0 = W_0(1,:).*[ones(1,Nx);U_0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Run models
fprintf('========================================\n')
fprintf('computing straightforward models... \n \n')

N = 10;
Tj_list = T_end.*[1/N:1/N:1];

% Run Shallow-Water model

tic
[SW_sol, dt_list, p] = Shallow_Water(W0,0,Tj_list);
SW_time = toc;

fprintf(['Elapsed time Shallow-Water : ',num2str(SW_time),'\n'])

fprintf(['-> initial total SW entropy : ',num2str(sum(entropy(SW_sol{1}))),'\n'])
fprintf(['-> final   total SW entropy : ',num2str(sum(entropy(SW_sol{end}))),'\n'])

% Green-Nagdhi

tic
[GN_sol] = Green_Nagdhi(W0,0,Tj_list);
GN_time = toc;
fprintf('\n')
fprintf(['Elapsed time Green-Nagdhi : ',num2str(GN_time),'\n'])

fprintf(['-> initial total GN entropy : ',num2str(sum(entropy(GN_sol{1}))),'\n'])
fprintf(['-> final   total GN entropy : ',num2str(sum(entropy(GN_sol{end}))),'\n'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parareal procedure

fprintf('========================================\n')
fprintf('computing Parareal... \n')

[U_para, t_para] = Parareal_v4(W0,[0, Tj_list],@Green_Nagdhi,@Shallow_Water);
fprintf(['Theoritical elapsed time Parareal : ',num2str(t_para),'\n\n'])

% end do_computations
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===== plot initial condition =====

call_plot(W_0,1)

% =====

figure(2)
subplot(2,1,1)
plot(x,GN_sol{end}(1,:),x,SW_sol{end}(1,:))
ylim([0,2*epsilon])
subplot(2,1,2)
plot(x,GN_sol{end}(2,:),x,SW_sol{end}(2,:))

% =====

figure(3)
subplot(2,1,1)
plot_water(x,SW_sol{end}(1,:),SW_sol{end}(2,:)./SW_sol{end}(1,:))
title('Shallow Water','interpreter','latex')
subplot(2,1,2)
plot_water(x,GN_sol{end}(1,:),GN_sol{end}(2,:)./GN_sol{end}(1,:))
title('Green-Naghdi','interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Entropy polts

lines_plt = round(sqrt(N+1));
colum_plt = round(sqrt(N+1)) + 1 - round(sqrt(N+1) - floor(sqrt(N+1)));

Ent_1 = zeros(N+1,N+1);
Ent_2 = Ent_1;

figure(7)
for j = 1:N+1
    for k = 1:N+1
    Ent_1(k,j) = sum(entropy(U_para{k,j}));
    end
    subplot(lines_plt,colum_plt,j)
    semilogy(1:N+1,Ent_1(:,j),'-o',...
        1:N+1,ones(1,N+1).*sum(entropy(W0)),'-r',...
        'LineWidth',Linesize,'MarkerSize',Marksize)
    set(gca,'FontSize',Fontsize_axes)
    xlabel('Parareal iteration nb','interpreter','latex','Fontsize',Fontsize_label)
    title(['entropy at $t = T_', num2str(j),'$'],'interpreter','latex','Fontsize',Fontsize)
    grid on
end

% ===== entropy at each steap

figure(8)
for j = 1:N+1
    for k = 1:N+1
        Ent_2(j,k) = sum(entropy(U_para{j,k}));
    end
    subplot(lines_plt,colum_plt,j)
    semilogy(0:N,Ent_2(j,:),'-o',...
        1:N+1,ones(1,N+1).*sum(entropy(W0)),'-r',...
        'LineWidth',Linesize,'MarkerSize',Marksize)
    set(gca,'FontSize',Fontsize_axes)
    xlim([0.14,inf])
    xlabel('Parareal iteration nb','interpreter','latex','Fontsize',Fontsize_label)
    grid on
end
sgtitle(['entropy at step ', num2str(j)],'interpreter','latex','Fontsize',Fontsize)

% ===== error at each time

init_zeros = zeros(N+1,N+1);

error_L2_h        = init_zeros;
error_L2_u        = init_zeros;
error_entropy     = init_zeros;
SWerror_L2_h      = init_zeros;
SWerror_L2_u      = init_zeros;
baseline_err_L2_h = init_zeros;
baseline_err_L2_u = init_zeros;

for j = 1:N+1 %(k : iteration, j : timestamp)
    for k = 1:N+1
        error_L2_h(k,j) = norm(GN_sol{j}(1,:) - U_para{k,j}(1,:),2);
        error_L2_u(k,j) = norm(GN_sol{j}(2,:)./GN_sol{j}(1,:) - U_para{k,j}(2,:)./U_para{k,j}(1,:),2);
        SWerror_L2_h(k,j) = norm(SW_sol{j}(1,:) - U_para{k,j}(1,:),2);
        SWerror_L2_u(k,j) = norm(SW_sol{j}(2,:)./SW_sol{j}(1,:) - U_para{k,j}(2,:)./U_para{k,j}(1,:),2);
        error_entropy(k,j) = err_entropy(GN_sol{j},U_para{k,j});
        %error_Linf_h(k,j) = norm(GN_sol{j}(1,:) - U{k,j}(1,:),Inf);
        %error_Linf_u(k,j) = norm(GN_sol{j}(2,:)./GN_sol{j}(1,:) - U{k,j}(2,:)./U{k,j}(1,:),Inf);
    end
    baseline_err_L2_h(j) = norm(GN_sol{j}(1,:) - SW_sol{j}(1,:),2);
    baseline_err_L2_u(j) = norm(GN_sol{j}(2,:)./GN_sol{j}(1,:) - SW_sol{j}(2,:)./SW_sol{j}(1,:),2);
    %baseline_err_Linf_h(j) = norm(GN_sol{j}(1,:) - SW_sol{j}(1,:),inf);
    %baseline_err_Linf_u(j) = norm(GN_sol{j}(2,:)./GN_sol{j}(1,:) - SW_sol{j}(2,:)./SW_sol{j}(1,:),inf);
end

figure(9)
for j = 1:N+1
    subplot(lines_plt,colum_plt,j)
    semilogy(0:N,error_L2_h(:,j),'-o',...
        1:N+1,ones(1,N+1).*baseline_err_L2_h(j),'-r',...
        1:N+1,SWerror_L2_h(:,j),'-x',...
        'LineWidth',Linesize,'MarkerSize',Marksize)
    set(gca,'FontSize',Fontsize_axes)
    xlabel('Parareal iteration nb','interpreter','latex','Fontsize',Fontsize_label)
    legend('$\|h_{para} - h_{GN}\|_2$',...
        '$\|h_{SW} - h_{GN}\|_2$',...
        '$\|h_{para} - h_{SW}\|_2$',...
        'Fontsize',Fontsize,...
        'interpreter','latex','Location','southwest')
    title(['at $t = T_',num2str(j),'$'],'interpreter','latex','Fontsize',Fontsize) 
    grid on
end
sgtitle('error on $h$','interpreter','latex','Fontsize',Fontsize) 

figure(10)
for j = 1:N+1
    subplot(floor(sqrt(N+1))+1,floor(sqrt(N+1))+1,j)
    semilogy(0:N,error_entropy(:,j),'-o',...
        'LineWidth',Linesize,'MarkerSize',Marksize)
    set(gca,'FontSize',Fontsize_axes)
    % title('error entropy : $\mathcal{E}$(GN) - $\mathcal{E}$(Parareal)','interpreter','latex')
    xlabel('Parareal iteration nb','interpreter','latex','Fontsize',Fontsize_label)
    title(['at $t = T_{',num2str(j-1),'}$'],'interpreter','latex','Fontsize',Fontsize) 
    grid on
end
sgtitle('error entropy : $\|\mathcal{E}$(GN) - $\mathcal{E}$(Parareal)$\|_1$','interpreter','latex','Fontsize',Fontsize)

for j = 1:N+1
    figure(10+j)
    for k = 1:N+1
        subplot(floor(sqrt(N+1))+1,floor(sqrt(N+1))+1,k)
        plot(x,GN_sol{j}(1,:),x,U_para{k,j}(1,:),...
            'LineWidth',Linesize,'MarkerSize',Marksize)
        set(gca,'FontSize',Fontsize_axes)
        grid on
        title(['Parareal iteration $',num2str(k-1),'$'],'interpreter','latex','Fontsize',Fontsize)
    end
    sgtitle(['solution on $h$ at $t = T_{',num2str(j-1),'}$'],'interpreter','latex','Fontsize',Fontsize)
end

for j = 1:N+1
    figure(30+j)
    for k = 1:N+1
        subplot(floor(sqrt(N+1))+1,floor(sqrt(N+1))+1,k)
        plot(x,GN_sol{j}(2,:)./GN_sol{j}(1,:),x,U_para{k,j}(2,:)./U_para{k,j}(1,:),...
            'LineWidth',Linesize,'MarkerSize',Marksize)
        set(gca,'FontSize',Fontsize_axes)
        grid on
        title(['Parareal iteration $',num2str(k-1),'$'],'interpreter','latex','Fontsize',Fontsize)
    end
    sgtitle(['solution on $u$ at $t = T_{',num2str(j-1),'}$'],'interpreter','latex','Fontsize',Fontsize)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [W_j, dt_list, p] = Shallow_Water(W_0,T_start,Tj)

% % W_0 is (h,hu,hw,hs)

%     p = 0;
%     J = length(Tj);
%     W_j = cell(1,J+1);
%     W_j{1} = W_0;
    
%     W_sol = W_0;
%     t = T_start;
    
%     dt_list = [];
    
%     for i = 1:J
        
%         while t < Tj(i)

%             dt_end = Tj(i) - t;

%             % Solve Riemann on every faces
%             [F_sol,dt_r] = Solver_vec(W_sol,t,Tj(i));           
%             dt    = min(dt_r, dt_end);
            
%             dt_list = [dt_list, dt];
            
%             % Time integration on every faces
%             W_sol = Euler(W_sol,F_sol,dt);
                        
%             % update 
%             t     = t + dt;
            
%             p = p+1;
            
%         end
        
%         W_j{i+1} = W_sol;
        
%     end
%     %fprintf('p Shallow-Water = %d \n',p)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W_j, dt_list, p] = Shallow_Water_d(W_0,T_start,Tj)

    % W_0 is (h,hu,hw,hs)

    p = 0;
    J = length(Tj);
    W_j = cell(1,J+1);
    W_j{1} = W_0;
    
    W_sol = W_0;
    t = T_start;
    
    dt_list = [];
    
    for i = 1:J
        
        while t < Tj(i)

            dt_end = Tj(i) - t;

            % Solve Riemann on every faces
            [F_sol,dt_r] = Solver_vec(W_sol,t,Tj(i));           
            dt    = min(dt_r, dt_end);
            
            dt_list = [dt_list, dt];
            
            % Time integration on every faces
            W_sol = Euler(W_sol,F_sol,dt);
                        
            % update 
            t     = t + dt;
            
            p = p+1;
            
        end
        
        W_j{i+1} = W_sol;
        
    end
    %fprintf('p Shallow-Water = %d \n',p)
    W_j = W_j{end};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [W_j, dt_list, p] = Green_Nagdhi(W_0,T_start,Tj)

%     global Nx
%     global dx

%     %%%%%%%%%%%%%%%%%%%%%
%     % Computes the __
%     % 
%     % W_0 tel que :
%     % W_0(1,:) = h
%     % W_0(2,:) = hu
%     %%%%%%%%%%%%%%%%%%%%%
%     % U_0 = Proj(W_0(1,:), [W_0(2,:)./W_0(1,:);zeros(2,Nx)]);
%     % W_0 = W_0(1,:).*[ones(1,Nx);U_0];
%     %     
%     % p=0;
%     % [W_0(3,:),W_0(4,:)] = get_ws(W_0(1,:),W_0(2,:)./W_0(1,:));
%     % W_0(3,:) = W_0(3,:).*W_0(1,:);
%     % W_0(4,:) = W_0(4,:).*W_0(1,:);
 
%     p = 0;
%     J = length(Tj);
%     W_j = cell(1,J+1);
%     W_j{1} = W_0;

%     W_sol = W_0;
%     t = T_start;
    
%     for i = 1:J
    
%         while t < Tj(i)

%             dt_end = Tj(i) - t;
            
%             [F_sol,dt_r] = Solver_vec(W_sol,t,Tj(i));
            
%             dt    = min(dt_r, dt_end);

%             % Time integration on evry faces
%             W_sol = Euler(W_sol,F_sol,dt);

%             % Projection
%              U_sol = Proj(W_sol(1,:),zeros(3,Nx) + (W_sol(1,:) > dx^2).*(W_sol(2:4,:)./W_sol(1,:)));
%              W_sol = W_sol(1,:).*[ones(1,Nx);U_sol];

%             % update 
%             t     = t + dt;
            
%             p = p+1;
            
%         end
        
%     W_j{i+1} = W_sol;
    
%     end
%     %fprintf('p green = %d \n',p)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W_j, dt_list, p] = Green_Nagdhi_d(W_0,T_start,Tj)
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
    W_j = W_j{end};  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function U = Proj(hs,Us)

%     % WARNING : spdiags() call sparse() so maybe better to code everything
%     % with sparse or, with only 1 call of spdiags()

%     global Nx
%     global dx
    
%     kappa = ([hs(1),hs,hs(end)].^3)./3;
%     A = sparse(1:Nx,1:Nx,(kappa(1:end-2) + kappa(3:end)),Nx,Nx) +...
%     sparse(1:Nx-2,3:Nx,-kappa(3:end-2),Nx,Nx) +...
%     sparse(3:Nx,1:Nx-2,-kappa(3:end-2),Nx,Nx) +...
%     sparse([1, 2, Nx-1, Nx],[2, 1, Nx, Nx-1],[kappa(1), kappa(2), kappa(end-1), kappa(end)],Nx,Nx);
%     A = A./(4*dx^2);
%     A = sparse(1:Nx,1:Nx,hs,Nx,Nx) + A;

%     us = Us(1,:); ws = Us(2,:); ss = Us(3,:);
    
% %     h_kp1 = [hs(2:end),hs(end)];
% %     h_km1 = [hs(1),hs(1:end-1)];
% %     w_kp1 = [ws(2:end),ws(end)];
% %     w_km1 = [ws(1),ws(1:end-1)];
% %     s_kp1 = [ss(2:end),ss(end)];
% %     s_km1 = [ss(1),ss(1:end-1)];
    
%     beta = hs.*us + ([hs(2:end),hs(end)].^2.*([ws(2:end),ws(end)] + [ss(2:end),ss(end)]./sqrt(3)) - [hs(1),hs(1:end-1)].^2.*([ws(1),ws(1:end-1)] + [ss(1),ss(1:end-1)]./sqrt(3)))./(4*dx);
    
%     u = mldivide(A,beta')';
    
%     [w, s] = get_ws(hs,u);
    
%     U = [u;w;s];
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [Flux, lambda] = Riemann_Rusanov(W_L,W_R)
%     global g
%     global dx
    
%     u_l = 0 + (W_L(1) > (dx)^2)*W_L(2)/W_L(1);
%     u_r = 0 + (W_R(1) > (dx)^2)*W_R(2)/W_R(1);
    
%     lambda = max(abs(u_l) + sqrt(g*W_L(1)),abs(u_r) + sqrt(g*W_R(1)));
        
%     Flux = (Flux_SW(W_L) + Flux_SW(W_R) - lambda.*(W_R - W_L))/2;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [Flux, lambda] = Riemann_HLL(W_L,W_R)
%     global g
%     global dx
    
%     u_l = 0 + (W_L(1) > (dx)^2)*W_L(2)/W_L(1);
%     u_r = 0 + (W_R(1) > (dx)^2)*W_R(2)/W_R(1);
    
%     lambda_p = max(0,max(u_l + sqrt(g*W_L(1)),u_r + sqrt(g*W_R(1))));
%     lambda_m = max(0,-min(u_l - sqrt(g*W_L(1)),u_r - sqrt(g*W_R(1))));
        
%     Flux = (lambda_p*Flux_SW(W_L) + lambda_m*Flux_SW(W_R) - lambda_m*lambda_p.*(W_R - W_L))/(lambda_p+ lambda_m);
    
%     lambda = max(lambda_m,lambda_p);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [Flux, lambda] = Riemann_HLLC(W_L,W_R)
%     global g
%     global dx
    
%     u_l = 0 + (W_L(1) > (dx)^2)*W_L(2)/W_L(1);
%     u_r = 0 + (W_R(1) > (dx)^2)*W_R(2)/W_R(1);
    
%     lambda_p = max(0,max(u_l + sqrt(g*W_L(1)),u_r + sqrt(g*W_R(1))));
%     lambda_m = max(0,-min(u_l - sqrt(g*W_L(1)),u_r - sqrt(g*W_R(1))));
        
%     Flux = (lambda_p*Flux_SW(W_L) + lambda_m*Flux_SW(W_R) - lambda_m*lambda_p.*(W_R - W_L))/(lambda_p + lambda_m);
%     lambda = max(lambda_m,lambda_p);
    
%     Flux(3:4) = [W_L(3);  W_L(4)]/W_L(1)*max(0,Flux(1))+[W_R(3);  W_R(4)]/W_R(1)*min(0,Flux(1));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function Flux = Flux_SW(W)
%     global dx
%     global g
%     Flux = [0;0;0;0];
%     if W(1) > dx^2
%         Flux = [W(2); W(2)^2/W(1) + (g*(W(1))^2)/2; W(2)*W(3)/W(1);  W(2)*W(4)/W(1)];
%     end      
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [F_n,dt] = Solver(W_n,W_border,t,T_end,Riemann)

%     dt = T_end - t;

%     global dx
%     global Nx
        
%     % case f = 1
%     [F_n(:,1), lambda] = Riemann(W_border(:,1),W_n(:,1));
%     dt               = min(dt, dx/(2*lambda));
    
%     % case f = Nx+1
%     [F_n(:,Nx+1), lambda] = Riemann(W_n(:,Nx),W_border(:,2));
%     dt               = min(dt, dx/(2*lambda));
    
%     % case f in 2:Nx
%     for f = 2:1:Nx
%        [F_n(:,f), lambda] = Riemann(W_n(:,f-1),W_n(:,f));
%        dt               = min(dt, dx/(2*lambda));
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [F_faces,dt] = Solver_vec(W,t,T_end)

%     global dx
%     global Nx
%     global g

%     dt = T_end - t;
    
%     % compute w on border
%     W_ = [[1;-1;1;1].*W(:,1), W, [1;-1;1;1].*W(:,Nx)];

%     % compute flux on each cell
%     F_ = [W_(2,:); (W_(2,:).^2)./W_(1,:) + (g.*(W_(1,:)).^2)./2; W_(2,:).*W_(3,:)./W_(1,:);  W_(2,:).*W_(4,:)./W_(1,:)];

%     % compute lambda on faces
%     H_ = W_(1,:);
%     U_ = W_(2,:)./W_(1,:);
%     lambda_faces = max([abs(U_(1:Nx+1)) + sqrt(g.*H_(1:Nx+1)) ; abs(U_(2:Nx+2)) + sqrt(g.*H_(2:Nx+2))], [], 1);

%     % compute flux on faces
%     F_faces = (F_(:,1:Nx+1) + F_(:,2:Nx+2) - lambda_faces.*(W_(:,2:Nx+2) - W_(:,1:Nx+1)))/2;

%     % compute dt
%     dt = min([dt, dx./(2*lambda_faces)]);
    
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [W_np1] = Euler(W_n,F_n,dt)
%     global Nx
%     global dx
%     W_np1=W_n;
%     % case f = 1
%     W_np1(:,1) = W_np1(:,1) + (dt/dx).*F_n(:,1);
    
%     % case f = Nx+1
%     W_np1(:,Nx) = W_np1(:,Nx) - (dt/dx).*F_n(:,Nx+1);
    
%     %case f in 2:Nx
%     for f = 2:1:Nx
%         W_np1(:,f-1) = W_np1(:,f-1) - (dt/dx).*F_n(:,f);
%         W_np1(:,f) = W_np1(:,f) + (dt/dx).*F_n(:,f);
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [w, s] = get_ws(h,u)
%     global dx
    
%     w = - h.*([u(2:end),-u(end)] - [-u(1),u(1:end-1)])./(4*dx);
%     s = w./sqrt(3);
    
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function E = entropy(W)

%     global g
%     global Nx
    
%     h = W(1,:);
%     u = zeros(1,Nx) + (W(1,:) > 0).*W(2,:)./W(1,:);
%     w = zeros(1,Nx) + (W(1,:) > 0).*W(3,:)./W(1,:);
%     s = zeros(1,Nx) + (W(1,:) > 0).*W(4,:)./W(1,:);
%     E = (h.^2).*g/2 + (h./2).*(u.^2 + w.^2 + s.^2);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function err = err_entropy(W1,W2)

    global dx

    err = norm(entropy(W1) - entropy(W2),1)*dx;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_water(x,h,u)
    global epsilon

    Nx = length(x);
    X = [x(1:Nx-1);x(2:Nx);x(2:Nx);x(1:Nx-1)];
    Y = [h(1:Nx-1);h(2:Nx);zeros(2,Nx-1)];
    C = (u(:,Nx-1)+u(2:Nx))/2;
    %figure(i)
    patch(X,Y,C,'EdgeColor','none')
    colorbar
    ylim([0,2*epsilon])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [W, time] = Parareal_v3(U0,Tj,Fine,Coarse)
% % PARAREAL implementation of the parareal algorithm
% global Nx
% global dx

% N = length(Tj) - 1;

% W = cell(N+1,N+1);
% Go = cell(1,N+1);
% Gn = cell(1,N+1);
% F = cell(1,N+1);
% W{1,1} = U0;

% t1 = 0;
% t2 = 0;
% t3 = 0;
% time = zeros(1,N);

% % initial guess with G
% for n=1:N
%     tic
%     Go{n+1} = Coarse(W{1,n},Tj(n),Tj(n+1));
%     t1 = t1 + toc;
%     W{1,n+1} = Go{n+1}{end};
% end
% % parareal iteration

% for k=1:N
    
%     fprintf('k = %d \n',k)
    
%     for n=k:N
%         tic
%         F{n+1} = Fine(Proj2(W{k,n}),Tj(n),Tj(n+1));     % parallel with F
%         t2 = t2 + toc;
%     end
    
%     for j=1:k
%         W{k+1,j} = W{k,j};
%     end

% %     W{k+1,1} = W{1,1};

%     for n=k:N
%         tic
%         Gn{n+1} = Coarse(Proj2(W{k+1,n}),Tj(n),Tj(n+1));       % sequential with G
%         W{k+1,n+1} = F{n+1}{end} + Proj2(Gn{n+1}{end}) - Proj2(Go{n+1}{end});  % parareal update
%         t3 = t3 + toc;
%     end
    
%     Go = Gn;
    
%     time(k) = t1 + t2/N + t3;
    
% end

% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function W_out = Proj2(W_in)
%     global Nx
%     global dx
    
%     h = W_in(1,:);
%     U = zeros(3,Nx) + (h > dx^2).*(W_in(2:4,:)./h);
%     Pi_h = Proj(h,U);
%     W_out = h.*[ones(1,Nx);Pi_h];
    
% end