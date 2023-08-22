function [] = call_2plot(W_1,W_2,i)

global x

figure(i)
leg = ['h','u','w','s'];
for j = 1:4
subplot(2,2,j)
plot(x,W_1(j,:),'b-',x,W_2(j,:),'r-')
title(leg(j),'Interpreter','latex')
%ylim([-0.5,1.5])
grid on
end

end