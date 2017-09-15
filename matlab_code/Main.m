%********************************************************************
%Damage Evolution code for Intervertebral Discs
%Created for CFDRC project
%Author: Shruti Motiwale
%********************************************************************

clear all; %#ok<CLSCR>
tic
disc_num=2;

%% Load input data for civilians
n_days1=365*80; %Number of days
% var=0;
% if(var)
load('von_noach.mat'); 
n_steps=9000; %Number of steps in a day for a civilian
[meanS,  maxS, ampS, cyc] = RainflowCounting2(von_noach, disc_num, time);
% else
%n_miles=12;
%n_steps=n_miles*1900; %Number of steps in a day
%[meanS, maxS, ampS, cyc] = RainflowCounting('Walk_ACH.mat', disc_num);
% end
cycleData=[maxS ampS meanS cyc];

% Run Damage Model
[Dmech1, Dtotal1, Dage1, d1] = damage9(cycleData, n_days1, n_steps, 1.5, 0.5, 'mean');
[Dmech1_m, Dtotal1_m, Dage1_m, d1] = damage9(cycleData, n_days1, n_steps, 1.5, 0.5, 'lower limit');
[Dmech1_p, Dtotal1_p, Dage1_p, d1] = damage9(cycleData, n_days1, n_steps, 1.5, 0.5, 'upper limit');

%% Load input data for soldiers
%n_miles=12;
%n_steps=n_miles*1900; %Number of steps in a day
%n_days2=365*22; %Number of days
%[meanS, maxS, ampS, cyc] = RainflowCounting('Walk_ACH.mat', disc_num);
%cycleData=[maxS ampS meanS cyc];

%% Run Damage Model for non-zero initial damage
%[Dmech2, Dtotal2, Dage2, d2] = damage10(d1, cycleData, n_days2, n_steps, 1.5, 0.5, 'mean');
%[Dmech2_m, Dtotal2_m, Dage2_m, d2] = damage10(d1, cycleData, n_days2, n_steps, 1.5, 0.5, 'lower limit');
%[Dmech2_p, Dtotal2_p, Dage2_p, d2] = damage10(d1, cycleData, n_days2, n_steps, 1.5, 0.5, 'upper limit');
%% Plot Data for intial D=0

hold all;
h1=plot((1:100:n_days1)/365,Dtotal1(1:100:end)); set(h1, 'Linewidth', 3); % Total Damage
h2=plot((1:100:n_days1)/365,Dmech1(1:100:end)); set(h2, 'Linewidth', 3); % Mechanical Damage
h3=plot((1:100:n_days1)/365,Dage1(1:100:end)); set(h3, 'Linewidth', 3); % Aging Damage
%% Plot Data for non-0 initial D

% figure
hold on;
% h1=plot((1:10:n_days1+n_days2)/365,[Dtotal1(1:10:end),Dtotal2(1:10:end)],'color',[0,0.45,0.74]); set(h1, 'Linewidth', 3);
%h2=plot((1:10:n_days1+n_days2)/365,[Dmech1(1:10:end),Dmech2(1:10:end)],'color',[0.85,0.33,0.1]); set(h2, 'Linewidth', 3);
% h3=plot((1:n_days1+n_days2)/365,[Dage1(1:end),Dage2(1:end)],'k'); set(h3, 'Linewidth', 3);
% h3=plot((1:10:n_days1+n_days2)/365,[Dage1(1:10:end),Dage2(1:10:end)],'color',[0.93,0.69,0.13]); set(h3, 'Linewidth', 3);
ylim([0,1])
set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'XMinorTick', 'on','YMinorTick','on');
box on
ylabel('Damage')
xlabel('No. of years')
legend('Total Damage','Mechanical Damage','Aging Damage')

% figure
% hold on;
% h1=plot((1:100:n_days2)/365+n_days1/365,Dtotal2(1:100:end),'b'); set(h1, 'Linewidth', 3); 
% h2=plot((1:100:n_days2)/365+n_days1/365,Dmech2(1:100:end),'r'); set(h2, 'Linewidth', 3); 
% h3=plot((1:100:n_days2)/365+n_days1/365,Dage2(1:100:end),'k'); set(h3, 'Linewidth', 3);
% ylim([0,1])
% % h3=plot((1:100:n_days2)/365+n_days1/365,Dtotal2(1:100:end)-Dmech2(1:100:end),'k:');
% set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'XMinorTick', 'on','YMinorTick','on');
% box on
% ylabel('Damage')
% xlabel('No. of years')
% legend('Total Damage','Mechanical Damage','Aging')
% title('Damage due to ACH')
%% Experimental data
% load CervDegen1.mat
% h1=plot(Age, Damage,'o','color',[0.49,0.18,0.56]); set(h1, 'Linewidth' , 3);
% load LumbDegen1.mat
% h1=plot(Age, Damage,'o','color',[0.49,0.18,0.56]); set(h1, 'Linewidth' , 3);
% set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'XMinorTick', 'on','YMinorTick','on');
% load LumbDegen2.mat
% h1=plot(Age, Damage,'o','color',[0.49,0.18,0.56]); set(h1, 'Linewidth' , 3);
% set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'XMinorTick', 'on','YMinorTick','on');
% load Thomson1.mat
% h1=plot(Age, Damage,'o','color',[0.49,0.18,0.56]); set(h1, 'Linewidth' , 3);
% set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'XMinorTick', 'on','YMinorTick','on');
% load Thomson2.mat
% h1=plot(Age, Damage,'o','color',[0.49,0.18,0.56]); set(h1, 'Linewidth' , 3);
% ylim([0,1])
% set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'XMinorTick', 'on','YMinorTick','on');
% legend('Total Damage','Mechanical Damage','Aging','Disc Degeneration Data')
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
set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'XMinorTick', 'on','YMinorTick','on');
ylabel('Damage')
% xlabel('Age (No. of years)')
box on
xlabel('Number of years')
% legend('Total Damage','Mechanical Damage','Aging')
%%

%% R square correlation
%ED=damage;
%TD=Dtotal1(int16(Age*365));
%RSquared=corrcoef(TD,ED)
% hold off
% toc