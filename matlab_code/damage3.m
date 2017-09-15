function [dam_eod, Dtot] = damage3(matrix, p, n, r, b)
% matrix=cycle data
% p = number of days
% n = number of steps in a day
% r = recovery

%Define variables
a=7;
beta=3.09;
% b=1;
UTS=6.25;%MPa
M0=14.8;%MPa
age=90*365;%days
%Assuming 0.6 damage due to aging in 90 years
Ca=0.6/age;%aging constant
% a=0.9711;
% beta=3.0216;
% b=0.9;
% UTS=3.9;%MPa
% M0=14.7358;%MPa


%Rearranging cycles to arrange according to descending order of Smax
mat=sortrows(matrix,-1); %the negative sign is for descending order
Smax=mat(:,1);
Samp=mat(:,2);
Smean=mat(:,3);
% Approximation strategy: all cycles of a certain amplitude 
% in a day are applied at once
cycles=n*mat(:,4); 
m=length(Samp);


%Damage Evolution
alpha(1)=1-a*Smax(1)/(UTS-Smax(1));
Nf(1)=(1/(1-alpha(1)))*abs((Samp(1)/(M0*(1-b*Smean(1)))))^(-beta);
dam(1)=(cycles(1)/Nf(1))^(1/(1-alpha(1)));
% alpha_mod(1)=alpha(1);
% Nf_mod(1)=Nf(1);
Neff(1)=cycles(1);

for k=2:m
    alpha(k)=1-a*Smax(k)/(UTS-Smax(k));
    Nf(k)=(1/(1-alpha(k)))*(abs(Samp(k)/M0*(1-b*Smean(k))))^(-beta);
    Neff(k)=cycles(k)+Nf(k)*(Neff(k-1)/Nf(k-1))^((1-alpha(k))/(1-alpha(k-1)));
    dam(k) = (Neff(k)/Nf(k))^(1/(1-alpha(k)));
end
dam(k)=(1-r)*dam(k);
Nf(k)=Nf(k)*(1+r);
Dtot(1)=dam(k)+Ca;

for i=2:p
    for k=1:length(Samp)
        gi=m*(i-1)+k; %global index
        alpha(gi)=1-a*Smax(k)/(UTS-Smax(k));
        Nf(gi)=(1/(1-alpha(gi)))*(abs(Samp(k)/M0*(1-b*Smean(k))))^(-beta);
        Neff(gi)=cycles(k)+Nf(gi)*(Neff(gi-1)/Nf(gi-1))^((1-alpha(gi))/(1-alpha(gi-1)));
        dam(gi)=(Neff(gi)/Nf(gi))^(1/(1-alpha(gi)));
    end
    dam(gi)=(1-r)*dam(gi);
    Nf(gi)=Nf(gi-1)*(1+r);
    Dtot(i)=dam(gi)+Ca;
end


%Post-processing
dam_eod=dam(1:m:end);
% % steps=(1:m*i)/9*n;
% steps=(1:i)*n_steps;
% days=steps/n;

end

%NOTES:
%1. No case considered for alpha=1, as in the present scenario, alpha<1
%always. If a positive fatigue limit is considered, alpha=1 case will have
%to be considered.