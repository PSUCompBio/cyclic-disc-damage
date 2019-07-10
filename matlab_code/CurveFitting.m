%********************************************************************
%Curve fitting with Stress life data for damage Evolution function 
%based on non-linear damage evolution model by Chaboche and Lemaitre
%Author: Shruti Motiwale
%Code Status - Works well as of 05/02/2016
%Please update the code status and specify modifications
%Derived from CurveFitting_Nf_2.m
%Confidence intervals may not be correct
%Modified the functions for alpha and damage rate
%********************************************************************
clear all;
tic
%% Define data
sdata_a=[0.861,1.134,1.134,1.239,1.239,1.323,1.344,1.344,1.491];
ndata_a=[10000,10000,10000,371,4245,960,215,2720,195];
sdata_p=[2.185,2.2325,2.2375,2.6125,2.66,2.8975,2.9925,3.0875,3.23,4.3325];
ndata_p=[10000,1835,10000,14,1860,5,15,27,700,5];
% sdata=[sdata_a sdata_p];
% ndata=[ndata_a ndata_p];
% sdata=[1.239,1.239,1.323,1.344,1.344,1.491,2.2325,2.66,2.9925,3.0875,3.23];
% ndata=[371,4245,960,215,2720,195,1835,1860,15,27,700];
sdata=[1.239,1.239,1.323,1.344,1.344,1.491,2.2325,2.6125,2.8975,2.9925,3.0875,4.3325];
ndata=[371,4245,960,215,2720,195,1835,14,5,15,27,5];
% sdata1=sdata_a;
% ndata1=ndata_a;
UTS_a=3.9;
UTS_p=8.6;
UTS=6.25;%average UTS
sf_a=0;%fatigue strength
sf_p=0;

%% Log curve fitting
% 1/(1-alpha)=(UTS-Smax)/(a*Samp); Nf=(UTS-Smax)/(a*Samp)*(Samp/M0)^(-beta);
objfcn = @(v)log10((UTS-sdata)./(v(1)*sdata).*((sdata/v(2)).^(-v(3)))) - log10(ndata);
% objfcn = @(v)(-1/v(2))*log10(sdata/v(1))-log10(ndata); %SN Curve
% objfcn = @(v)log10((UTS-sdata)./(v(1)*sdata.^v(2))) - log10(ndata); % Reduced Parameter Chaboche Equation
% modelfun = @(v,sdata)log10((UTS-sdata)/(v(1)*sdata)*(sdata/v(2)).^(-v(3)));
% modelfun = @(v,sdata)log10(((UTS-sdata)./(v(1)*sdata)).*(v(2)./sdata).^(v(3)));
% modelfun = @(v,sdata)log10((UTS-sdata)./(v(1)*sdata.^v(2))); %Reduced Parameter Chaboche equation
% modelfun = @(v,sdata)(-1/v(2))*log10(sdata/v(1)); %SN Curve
opts = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
v0 = [5,50,2.5]; %Arbitrary initial guess, a,Mo,beta
% [v_de,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(sdata, log10(ndata), modelfun, v0);
[vde, resnorm, residuals, exitflag,output] = lsqnonlin(objfcn,v0,[0,0,0],[10,Inf,Inf],opts);
vde, resnorm, exitflag,output.firstorderopt %#ok<*NOPTS>
problem=createOptimProblem('lsqnonlin','x0',v0,'objective',objfcn,'lb',[0,0,0],'ub',[10,Inf,Inf],'xdata',sdata,'ydata',log10(ndata));
ms = MultiStart('PlotFcns',@gsplotbestf)
[vm_de, errormulti, exitflag, output, solutions]=run(ms,problem,100)
[v_de, resnorm, residuals, exitflag, output, lambda, jacobian] = lsqnonlin(objfcn,vm_de,[0,0,0],[10,Inf,Inf],opts);
v_de, resnorm, exitflag,output.firstorderopt %#ok<*NOPTS>
%%
% ci=nlparci(v_de,R,J)
ci=nlparci(v_de, residuals, jacobian)
for i=1:length(sdata)
%     nf(i)=(UTS-sdata(i))/(v_de(1)*sdata(i)^v_de(2));
%     nf(i)=(sdata(i)/v_de(1))^(-1/v_de(2));
    nf(i)=(UTS-sdata(i))/(v_de(1)*sdata(i))*(sdata(i)/v_de(2))^(-v_de(3));
end

%% Plot data points
figure(10)
hold all

% Plot damage evolution
h1 = plot(log10(nf),sdata,'r'); set(h1, 'Linewidth' , 3);
% h1 = plot(log10(nf),sdata,'r'); set(h1, 'Linewidth' , 3);
%Plot data
% plot(log(ndata_a)/log(10),sdata_a,'k*')
h1 = plot(log10(ndata),sdata,'o'); set(h1, 'Linewidth' , 3);
set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'XMinorTick', 'on','YMinorTick','on');
box on
legend('Damage evolution','Green''s data')
xlabel('log(N_f)')
ylabel('Stress (MPa)')
% plot(log(sort(ndata1))/log(10), 4.0264*(sort(ndata1)).^(-0.204))

% figure(13)
% hold all
% %Plot SN curve
% % h2 = plot(dm_sn,sdata_p,'g');
% % set(h2, 'Linewidth' , 3);
% % set(gca, 'FontSize', 20);
% xlabel('N')
% ylabel('Amplitude of Stress (Mpa)')
% % Plot damage evolution
% h1 = plot(nf,sdata,'r');
% set(h1, 'Linewidth' , 3);%Plot data
% % plot(ndata_a,sdata_a,'k*')
% h1 = plot(ndata,sdata,'ko')
% set(h1, 'Linewidth' , 3);
% legend('SN Curve','Damage evolution','Green''s data Posterior')
% plot(sort(ndata1), 4.0264*(sort(ndata1)).^(-0.204))
toc
save('curvefit.mat')

