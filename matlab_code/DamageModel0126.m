%********************************************************************
%Damage Evolution code based on non-linear damage evolution model by
%Chaboche and Lemaitre
%Author: Shruti motiwale
%********************************************************************

%Load input data
clear all;
clc;
load Run_noHelmet.mat
disc_no=2;
s=C(:,disc_no);
%disc=2; 

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

%Damage Model
a=0.9711;
beta=3.0216;
b=1;
UTS=3.9;%MPa
% sf=0.45*UTS; %Fatigue Limit
sf=0;
M0=14.7358;%MPa
% Smax=-Smax;
% Smean=-Smean;
% Samp=-Samp;
p=365*60; %Number of days
n=9000; %Number of steps in a day
%Rearranging cycles to arrange according to descending order of Smax
matrix=[Smax Samp Smean cycles];
mat=sortrows(matrix,-1);
Smax=mat(:,1);
Samp=mat(:,2);
Smean=mat(:,3);
cycles=n*mat(:,4); %each cycle is applied 6000 times at once
m=length(Samp);
r=0; %recovery measure

alpha(1)=1-a*max((Smax(1)-sf)/(UTS-Smax(1)),0);
if alpha(1)<1
    Nf(1)=(1/(1-alpha(1)))*abs((Samp(1)/(M0*(1-b*Smean(1)))))^(-beta);
%     Nf(1)=(1+r)*Nf(1);
%     recovery
    dam(1)=(cycles(1)/Nf(1))^(1/(1-alpha(1)));
    dam(1)=dam(1)*(1-r);
    Nf(1)=cycles(1)/(dam(1)^(1-alpha(1)));
    alpha_mod(1)=alpha(1);
    Nf_mod(1)=Nf(1);
    Neff(1)=cycles(1);
else
    dam(1)=0;
    alpha_mod(1)=0.9999;
    Nf_mod(1)=(1/(1-alpha_mod(1)))*(abs(Samp(1)/M0))^(-beta);
    Neff(1)=cycles(1);
end


for k=2:m
        alpha(k)=1-a*max((Smax(k)-sf)/(UTS-Smax(k)),0);
        if alpha(k)<1
            Nf(k)=(1/(1-alpha(k)))*(abs(Samp(k)/M0*(1-b*Smean(k))))^(-beta);
            Nf(k)=(1+r)*Nf(k);
            Neff(k)=cycles(k)+Nf(k)*(Neff(k-1)/Nf_mod(k-1))^((1-alpha(k))/(1-alpha_mod(k-1)));
            dam(k) = (Neff(k)/Nf(k))^(1/(1-alpha(k)));
            %recovery
            dam(k)=dam(k-1)+(dam(k)-dam(k-1))*(1-r);
            Nf(k)=Neff(k)/(dam(k)^(1-alpha(k)));
            alpha_mod(k)=alpha(k);
            Nf_mod(k)=Nf(k);
            Neff(k)=cycles(k);
        else
            dam(k)=dam(k-1)/((1+r)^(1/(1-alpha(k))));%damage reduced and Nf increased by a factor of (1+r) overnight
            alpha_mod(k)=alpha_mod(k-1);
            Nf_mod(k)=Nf_mod(k-1)*(1+r);
            Neff(k)=Neff(k-1);
        end
%         dam(k)=(1-r)*dam(k);
end

i=1;
for i=2:p
    for k=1:length(Samp)
        gi=m*(i-1)+k; %global index
        alpha(gi)=1-a*max((Smax(k)-sf)/(UTS-Smax(k)),0);
        if alpha(gi)<1
            Nf(gi)=(1/(1-alpha(gi)))*(abs(Samp(k)/M0*(1-b*Smean(k))))^(-beta);
            Nf(gi)=(1+r)*Nf(gi);
            Neff(gi)=cycles(k)+Nf(gi)*(Neff(gi-1)/Nf_mod(gi-1))^((1-alpha(gi))/(1-alpha_mod(gi-1)));
            dam(gi)=(Neff(gi)/Nf(gi))^(1/(1-alpha(gi)));
            %recovery
            dam(gi)=dam(gi-1)+(dam(gi)-dam(gi-1))*(1-r);
%             dam(gi)=dam(gi-1)+(dam(gi)-dam(gi-1));
            Nf(gi)=Neff(gi)/(dam(gi)^(1-alpha(gi)));
            alpha_mod(gi)=alpha(gi);
            Nf_mod(gi)=Nf(gi);
        else
            dam(gi)=dam(gi-1)/((1+r)^(1/(1-alpha(k))));
            alpha_mod(gi)=alpha_mod(gi-1);
            Nf_mod(gi)=Nf_mod(gi-1)*(1+r);
            Neff(gi)=Neff(gi-1);
        end
%         dam(gi)=(1-r)*dam(gi);
    end
end

figure(10)
hold on;
reduced_dam=dam(1:m:end);
% steps=(1:m*i)/9*n;
steps=(1:i)*n;
days=steps/n;
h1=plot(days(1:400:end)/365,reduced_dam(1:400:end));
% h1=plot(reduced_dam,days);
set(h1, 'Linewidth' , 3);
set(gca, 'FontSize', 20);
ylabel('Damage')
xlabel('Number of years')
% legend('b=1','b=0.5','b=0.1')
