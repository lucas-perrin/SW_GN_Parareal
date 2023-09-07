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