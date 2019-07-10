function [Smean, Smax, Samp, cycles] = RainflowCounting(fname,disc_no)

load(fname)
s=C(:,disc_no);
%Cycle Extraction from input data
[tp,te]=sig2ext(s,Time);
rf=rainflow(tp,te);
% rfm=rfmatrix(rf); %not sure yet what the bins are
% figure(14)
% plot(rfm);

Famp=rf(1,:).';
Fmean=rf(2,:).';

area=312; %mm^2;
Smean=Fmean/area; %MPa
Samp=Famp/area;   %MPa
Smax=Smean+Samp;
cycles=rf(3,:).';