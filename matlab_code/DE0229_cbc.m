%********************************************************************
%Damage Evolution code based on non-linear damage evolution model by
%Chaboche and Lemaitre
%Author: Shruti Motiwale
%Code Status - 
%Please update the code status and specify modifications
%No Approximations
%********************************************************************
%Algorithm: If current alpha is 1, make damage due to current cycle zero.
%Also update current alpha, Nf and number of cycles to previous alpha and damage, so that you
%can use them to calculate damage for next iteration. let us store these
%modified vectors in alpha_mod, Nf_mod, cycles_mod etc.

clc;
clear all;%close all;

%load('Walk_ACH.mat');
%Define variables
a=4.8898; beta=2.9357; M0=21.9419; b=1;
UTS=6.25;%MPa

%disc_no=2;
%s=C(:,disc_no);
Smax=1; Samp=1; Smean=0; cycles=6;
%[Smean, Smax, Samp, cycles] = RainflowCounting('Walk_ACH.mat', Time);
m=length(Samp);


p=365*0.2; %number of days %try for 5 years
q=9000; %number of steps in a day
n=q*p; %Number of steps
m=length(Samp);
r=0; %recovery measure

alpha(1)=1-a*Smax(1)/(UTS-Smax(1));
if alpha(1)<1
    Nf(1)=(1/(1-alpha(1)))*abs((Samp(1)/(M0*(1-b*Smean(1)))))^(-beta);
    dam(1)=(cycles(1)/Nf(1))^(1/(1-alpha(1)));
    alpha_mod(1)=alpha(1);
    Nf_mod(1)=Nf(1);
    Neff(1)=cycles(1);
else
    dam(1)=0;
    alpha_mod(1)=0.9999;
    Nf_mod(1)=(1/(1-alpha_mod(1)))*(abs(Samp(1)/M0))^(-beta);
    Neff(1)=cycles(1);
end

if m>1
    for k=2:m
        alpha(k)=1-a*Smax(k)/(UTS-Smax(k));
        if alpha(k)<1
            Nf(k)=(1/(1-alpha(k)))*(abs(Samp(k)/M0*(1-b*Smean(k))))^(-beta);
            Neff(k)=cycles(k)+Nf(k)*(Neff(k-1)/Nf_mod(k-1))^((1-alpha(k))/(1-alpha_mod(k-1)));
            dam(k) = (Neff(k)/Nf(k))^(1/(1-alpha(k)));
            alpha_mod(k)=alpha(k);
            Nf_mod(k)=Nf(k);
            %             Neff(k)=cycles(k);
        else
            dam(k)=dam(k-1);
            alpha_mod(k)=alpha_mod(k-1);
            Nf_mod(k)=Nf_mod(k-1);
            Neff(k)=Neff(k-1);
        end
    end
else
    k=1;
end
dam(k)=(1-r)*dam(k);
Nf(k)=Neff(k)/(dam(k)^(1-alpha(k)));
Dtot(1)=dam(k);%No aging occurs in the first cycle as aging starts after the age of 20
Dmech(1)=dam(k);

for i=2:n
    for k=1:length(Samp)
        gi=m*(i-1)+k; %global index
        alpha(gi)=1-a*Smax(k)/(UTS-Smax(k));
        if alpha(gi)<1
            Nf(gi)=(1/(1-alpha(gi)))*(abs(Samp(k)/M0*(1-b*Smean(k))))^(-beta);
            Neff(gi)=cycles(k)+Nf(gi)*(Neff(gi-1)/Nf_mod(gi-1))^((1-alpha(gi))/(1-alpha_mod(gi-1)));
            dam(gi)=(Neff(gi)/Nf(gi))^(1/(1-alpha(gi)));
            alpha_mod(gi)=alpha(gi);
            Nf_mod(gi)=Nf(gi);
        else
            dam(gi)=dam(gi-1);
            alpha_mod(gi)=alpha_mod(gi-1);
            Nf_mod(gi)=Nf_mod(gi-1);
            Neff(gi)=Neff(gi-1);
        end
    end
    dam(gi)=(1-r)*dam(gi);
    Nf(gi)=Neff(gi)/(dam(gi)^(1-alpha(gi)));
    Dmech(i)=dam(gi);
%     if i<=365*thr_age
%         Dtot(i)=Dmech(i);
%     else
%         Dtot(i)=dam(gi)+Ca*(i-thr_age*365);
%     end
end

reduced_dam=dam(1:9:end);
steps=1:i;
days=steps/q;

h1=semilogx((1:n)*cycles,dam(1:end)); set(h1, 'Linewidth' , 5);
hold all
ylim([0,1])
%h1=plot(days/365, reduced_dam,'y');
set(gca, 'FontSize', 22);
ylabel('Damage')
%xlabel('Number of years')
xlabel('Number of cycles')