clear all;%close all; %#ok<CLSCR>

% n_days1=365*1; %Number of days
UTS=6.25;
a=4.8898; beta=2.9357; M0=21.9419; 
% b=1.5;
% n_steps=9000; %Number of steps in a day

% disc_num=2;
% [Smean, Smax, Samp, cyc] = RainflowCounting('Walk_noHelmet.mat', disc_num);
% [Smean, Smax, Samp, cyc] = RainflowCounting('Walk_ACH.mat', disc_num);
Smax=1; Samp=1; Smean=0; %cyc=6;
% cycleData=[Smax Samp Smean cyc];
alpha=1-a*Smax./(UTS-Smax);
% alpha=1-a*Samp./(UTS-Smax);
% m1=(Samp./(M0*(1-b*Smean))).^beta;
% m2=sum(cyc)*n_steps*365;%#cycles in a step*#steps in a day*#days in a year
% m=m1*m2;
% Compare=[Smax Smean Samp alpha m];
dam=0;
n=linspace(0, 2500000, 2500001);
dam=zeros(1, 2500001);

for i=1:2500001
    Nf=1/(1-alpha)*(Smax/M0)^(-beta);
    dam(i)=((n(i)-1)/Nf).^(1/(1-alpha));
    
end


hold all
h2=semilogx(n,dam, '--'); set(h2, 'Linewidth', 3); 
set(gca, 'FontSize', 22, 'LineWidth', 1.5, 'XMinorTick', 'on','YMinorTick','on');
ylim([0, 1]);
hold all
box on
