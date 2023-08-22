function [] = call_plot(W,i)

global x

figure(i)
leg = ['h','u','w','s'];
for j = 1:4
subplot(2,2,j)
plot(x,W(j,:),'-')
title(leg(j),'Interpreter','latex')
%ylim([-0.5,1.5])
grid on
end

end