function [Smean, Smax, Samp, cycles] = RainflowCounting2(von,disc_no, time)

%load(fname)
s=von(:, 1);
%Cycle Extraction from input data
[tp,te]=sig2ext(s,time);
rf=rainflow(tp,te);
% rfm=rfmatrix(rf); %not sure yet what the bins are
% figure(14)
% plot(rfm);

Samp=rf(1,:).';
Smean=rf(2,:).';

%area=312; %mm^2;
%Smean=Fmean/area; %MPa
%Samp=Famp/area;   %MPa
Smax=Smean+Samp;
cycles=rf(3,:).';