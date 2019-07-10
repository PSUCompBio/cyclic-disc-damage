%********************************************************************
%Damage Evolution code for Intervertebral Discs
%Created for CFDRC project
%Author: Shruti Motiwale
%********************************************************************

clear all; 
tic
disc_num=2;

%% Load input data for civilians
n_days1=365*30; %Number of days
load('von_ach');
n_steps=22800; %Number of steps per day  
%j=1;
%legend_string=strings;
%for i=2:1788
%von=mises_ach{2:end, i};
%von=von.*0.000001;
von=von_ach; 
[meanS,  maxS, ampS, cyc] = RainflowCounting2(von, disc_num, time);
cycleData=[maxS ampS meanS cyc];
% Run Damage Model
[Dmech1, Dtotal1, Dage1, d1] = damage9(cycleData, n_days1, n_steps, 1.5, 0.5, 'mean');
[Dmech1_m, Dtotal1_m, Dage1_m, d1] = damage9(cycleData, n_days1, n_steps, 1.5, 0.5, 'lower limit');
[Dmech1_p, Dtotal1_p, Dage1_p, d1] = damage9(cycleData, n_days1, n_steps, 1.5, 0.5, 'upper limit');
%% Plot Data for intial D=0
hold all;
h1=plot((1:100:n_days1)/365,Dtotal1(1:100:end)); set(h1, 'Linewidth', 3); % Total Damage
h2=plot((1:100:n_days1)/365,Dmech1(1:100:end)); set(h2, 'Linewidth', 3); % Mechanical Damage
h3=plot((1:100:n_days1)/365,Dage1(1:100:end)); set(h3, 'Linewidth', 3); % Aging Damage
hold on;
ylim([0,1])
set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'XMinorTick', 'on','YMinorTick','on');
box on
ylabel('Damage')
xlabel('Number of years')
%% Plot Confidence Intervals
% plot_ci((1:100:n_days1)'/365,[Dage1(1:100:end)',Dage1_m(1:100:end)',Dage1_p(1:100:end)'],...
%     'PatchColor', 'k', 'PatchAlpha', 0.1, 'MainLineWidth', 3, 'MainLineStyle', '-', 'MainLineColor', [0.93,0.69,0.13],'LineWidth', 2, 'LineStyle','--', 'LineColor', [0.93,0.69,0.13]);
% plot_ci((1:100:n_days1)'/365, [Dmech1(1:100:end)'+Dage1(1:100:end)',...
%     Dmech1_m(1:100:end)'+Dage1_m(1:100:end)',Dmech1_p(1:100:end)'+Dage1_p(1:100:end)'],...
%     'PatchColor', 'k', 'PatchAlpha', 0.1, 'MainLineWidth', 3, 'MainLineStyle', '-', 'MainLineColor', [0,0.45,0.74],'LineWidth', 2, 'LineStyle','--', 'LineColor', [0,0.45,0.74]);

%plot_ci((1:10:n_days1+n_days2)'/365,[[Dage1(1:10:end),Dage2(1:10:end)]',...
%    [Dage1_m(1:10:end),Dage2_m(1:10:end)]',[Dage1_p(1:10:end),Dage2_p(1:10:end)]'],...
%    'PatchColor', 'k', 'PatchAlpha', 0.1, 'MainLineWidth', 3, 'MainLineStyle', '-', 'MainLineColor', [0.93,0.69,0.13],'LineWidth', 2, 'LineStyle','--', 'LineColor', [0.93,0.69,0.13]);
%plot_ci((1:10:n_days1+n_days2)'/365,...
%    [[Dmech1(1:10:end),Dmech2(1:10:end)]'+[Dage1(1:10:end),Dage2(1:10:end)]',...
%    [Dmech1_m(1:10:end),Dmech2_m(1:10:end)]'+[Dage1_m(1:10:end),Dage2_m(1:10:end)]',...
%    [Dmech1_p(1:10:end),Dmech2_p(1:10:end)]'+[Dage1_p(1:10:end),Dage2_p(1:10:end)]'],...
%    'PatchColor', 'k', 'PatchAlpha', 0.1, 'MainLineWidth', 3, 'MainLineStyle', '-', 'MainLineColor', [0,0.45,0.74],'LineWidth', 2, 'LineStyle','--', 'LineColor', [0,0.45,0.74]);
%% Set Plot Properties
% h3=plot((1:100:n_days1)/365,Dtotal1(1:100:end)-Dmech1(1:100:end)); set(h3, 'Linewidth', 3);
% ylim([0,1]);
%set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'XMinorTick', 'on','YMinorTick','on');
%ylabel('Damage')
% xlabel('Age (No. of years)')
%box on
%xlabel('Number of years')
% legend('Total Damage','Mechanical Damage','Aging')
%% R square correlation
%ED=damage;
%TD=Dtotal1(int16(Age*365));
%RSquared=corrcoef(TD,ED)
%dam_array(:, i)=Dtotal1;
%end