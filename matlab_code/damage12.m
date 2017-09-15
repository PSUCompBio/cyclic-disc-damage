%********************************************************************
%Damage Evolution function based on non-linear damage evolution model by
%Chaboche and Lemaitre
%Author: Shruti Motiwale
%Code Status - Works well as of 05/02/2016
%Please update the code status and specify modifications
%Derived from Damage 4.
%Aging function obtained from literature data calibration
%Changed the function for alpha
%********************************************************************

function [Dmech, Dtot, Da, d_out] = damage12(matrix, p, n, b, r)
% matrix=cycle data
% p = number of days
% n = number of steps in a day
% r = recovery

%Define variable
a=9.8865; beta=2.9737; M0=17.8154;
% b=0.9;
UTS=6.25;%MPa
% sf=0.45*UTS; %Fatigue Limit
% sf=0;


%Rearranging cycles to arrange according to descending order of Smax
% matrix=[Smax Samp Smean cycles];
mat=sortrows(matrix,-1); %the negative sign is for descending order
Smax=mat(:,1);
Samp=mat(:,2);
Smean=mat(:,3);
% Approximation strategy: all cycles of a certain amplitude 
% in a day are applied at once
cycles=n*mat(:,4); 
m=length(Samp);


%Damage Evolution
alpha(1)=1-a*Samp(1)/(UTS-Smax(1));
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
        alpha(k)=1-a*Samp(k)/(UTS-Smax(k));
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
Dtot(1)=dam(k);
Dmech(1)=dam(k);

for i=2:p
    for k=1:length(Samp)
        gi=m*(i-1)+k; %global index
        alpha(gi)=1-a*Samp(k)/(UTS-Smax(k));
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
    Da(i)=0.0117*i/365-0.2614;
    if Da(i)<0
        Da(i)=0;      
    end
    Dtot(i)=Dmech(i)+Da(i);
end
d_out=[Neff(end), Nf_mod(end), alpha_mod(end), i];
end