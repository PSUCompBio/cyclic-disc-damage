% h=semilogx(n_3,d_3); set(h, 'Linewidth', 3); 
% set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'XMinorTick', 'on','YMinorTick','on');
% ylim([0, 1]);
% hold all
% box on

close all;
alpha=zeros(1, 9); mult=zeros(1, 9);
a=4.8898; beta=2.9357; M0=21.9419; UTS=6.25; b=1.5;
for i=1:9
    alpha(i)=1-a*maxS(i)/(UTS-maxS(i));
    mult(i)=ampS(i)/(M0*(1-b*meanS(i)))^beta;
end
figure
plot(alpha);
set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'XMinorTick', 'on','YMinorTick','on');
figure
plot(mult);
set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'XMinorTick', 'on','YMinorTick','on');