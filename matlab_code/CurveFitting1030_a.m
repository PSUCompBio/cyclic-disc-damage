clear all;
clc;

%%
sdata_a=[0.861,1.134,1.134,1.239,1.239,1.323,1.344,1.344,1.491];
ndata_a=[10000,10000,10000,371,4245,960,215,2720,195];
sdata_p=[2.185,2.2325,2.2375,2.6125,2.66,2.8975,2.9925,3.0875,3.23,4.3325];
ndata_p=[10000,1835,10000,14,1860,5,15,27,700,5];
% sdata1=[sdata_a sdata_p];
% ndata1=[ndata_a ndata_p];
% sdata1=sdata_a;
% ndata1=ndata_a;
UTS_a=3.9;
UTS_p=8.6;
sf_a=0;%fatigue strength
sf_p=0;

%% Fit coefficients to SN curve for power law
objfcn = @(v1)(sdata_a/v1(1)).^(-v1(2)) - ndata_a; 
opts = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective')
v0 = [4,6.1451]; %Arbitrary initial guess
[vsn, resnorm, residuals, exitflag,output] = lsqnonlin(objfcn,v0,[0,0],[Inf,Inf],opts);
vsn, resnorm, exitflag,output.firstorderopt %#ok<*NOPTS>
sn = (sdata_a/vsn(1)).^(-vsn(2));
problem=createOptimProblem('lsqnonlin','x0',v0,'objective',objfcn,'lb',[0,0],'ub',[Inf,Inf],'xdata',sdata_a,'ydata',ndata_a);
ms = MultiStart('PlotFcns',@gsplotbestf)
[vm_sn,errormulti]=run(ms,problem,1500)
for i=1:length(sdata_a)
    dm_sn(i)=(sdata_a(i)/vm_sn(1)).^(-vm_sn(2));
end

%% Fit coefficients for simplified damage equation
%v2(1)=a, v2(2)=M0, v2(3)=beta, sf=0.45*UTS
objfcn = @(v2)1/v2(1)*(UTS_a-sdata_a)/max((sdata_a-sf_a),0)*(sdata_a/v2(2)).^(-v2(3)) - ndata_a;
opts = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective')
v0 = [0.5,4,0.5]; %Arbitrary initial guess
[vde, resnorm, residuals, exitflag,output] = lsqnonlin(objfcn,v0,[0,0,0],[],opts);
vde, resnorm, exitflag,output.firstorderopt %#ok<*NOPTS>
n_de = 1/vde(1)*(UTS_a-sdata_a)/max((sdata_a-0.45*UTS_a),0)*(sdata_a/vde(2)).^(-vde(3));
problem=createOptimProblem('lsqnonlin','x0',v0,'objective',objfcn,'lb',[0.5,0,0],'ub',[1,+Inf,+Inf],'xdata',sdata_a,'ydata',ndata_a);
ms = MultiStart('PlotFcns',@gsplotbestf)
[vm_de,errormulti]=run(ms,problem,1500)
for i=1:length(sdata_a)
    dm_de(i)=1/vm_de(1)*(UTS_a-sdata_a)/max((sdata_a-sf_a),0)*(sdata_a(i)/vm_de(2))^(-vm_de(3));
end

% Plot data points
% figure(8)
% hold all
% %Plot SN curve
% h2 = plot(log(dm_sn)/log(10),sdata_a,'g');
% set(h2, 'Linewidth' , 3);
% set(gca, 'FontSize', 20);
% xlabel('log(N)')
% ylabel('Amplitude of Stress (Mpa)')
% % Plot damage evolution
% h1 = plot(log(dm_de)/log(10),sdata_a,'r-.');
% set(h1, 'Linewidth' , 3);
% %Plot data
% plot(log(ndata_a)/log(10),sdata_a,'k*')
% % plot(log(ndata_p)/log(10),sdata_p,'b*')
% legend('SN Curve','Damage evolution','Green''s data Anterior')
% % plot(log(sort(ndata1))/log(10), 4.0264*(sort(ndata1)).^(-0.204))
% hold off

figure(13)
hold all
%Plot SN curve
h2 = plot(dm_sn,sdata_a,'g');
set(h2, 'Linewidth' , 3);
set(gca, 'FontSize', 20);
xlabel('N')
ylabel('Amplitude of Stress (Mpa)')
% Plot damage evolution
h1 = plot(dm_de,sdata_a,'r-.');
set(h1, 'Linewidth' , 3);%Plot data
plot(ndata_a,sdata_a,'k*')
% plot(ndata_p,sdata_p,'b*')
legend('SN Curve','Damage evolution','Green''s data Anterior')
% plot(sort(ndata1), 4.0264*(sort(ndata1)).^(-0.204))
hold off
save('curvefit_a.mat')

