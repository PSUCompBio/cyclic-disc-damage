function [dam_eod] = damage2(matrix, p, n, r, b)
% matrix=cycle data
% p = number of days
% n = number of steps in a day
% r = recovery

%Define variables
a=7;
beta=3.09;
% b=1;
UTS=6.25;%MPa
% sf=0.45*UTS; %Fatigue Limit
sf=0;
M0=14.8;%MPa
age=90*365;%days
%Assuming 0.6 damage due to aging in 90 years
Ca=0.6/age;%aging constant
% a=0.9711;
% beta=3.0216;
% b=0.9;
% UTS=3.9;%MPa
% % sf=0.45*UTS; %Fatigue Limit
% sf=0;
% M0=14.7358;%MPa
% Smax=-Smax;
% Smean=-Smean;
% Samp=-Samp;

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
alpha(1)=1-a*max((Smax(1)-sf)/(UTS-Smax(1)),0);
if alpha(1)<1
    Nf(1)=(1/(1-alpha(1)))*abs((Samp(1)/(M0*(1-b*Smean(1)))))^(-beta);
%     Nf(1)=(1+r)*Nf(1);
%     recovery
    dam(1)=(cycles(1)/Nf(1))^(1/(1-alpha(1)));
%     dam(1)=dam(1)*(1-r);
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
%             Nf(k)=(1+r)*Nf(k);
            Neff(k)=cycles(k)+Nf(k)*(Neff(k-1)/Nf_mod(k-1))^((1-alpha(k))/(1-alpha_mod(k-1)));
            dam(k) = (Neff(k)/Nf(k))^(1/(1-alpha(k)));
            %recovery should not be occuring after every block of loading
%             dam(k)=dam(k-1)+(dam(k)-dam(k-1))*(1-r);
%             Nf(k)=Neff(k)/(dam(k)^(1-alpha(k)));
            alpha_mod(k)=alpha(k);
            Nf_mod(k)=Nf(k);
            Neff(k)=cycles(k);
        else
%             dam(k)=dam(k-1)/((1+r)^(1/(1-alpha(k))));%damage reduced and Nf increased by a factor of (1+r) overnight
            dam(k)=dam(k-1);
            alpha_mod(k)=alpha_mod(k-1);
            Nf_mod(k)=Nf_mod(k-1);
            Neff(k)=Neff(k-1);
        end
        dam(k)=(1-r)*dam(k);
        Nf_mod(k)=Nf_mod(k)*(1+r);
end

i=1;
for i=2:p
    for k=1:length(Samp)
        gi=m*(i-1)+k; %global index
        alpha(gi)=1-a*max((Smax(k)-sf)/(UTS-Smax(k)),0);
        if alpha(gi)<1
            Nf(gi)=(1/(1-alpha(gi)))*(abs(Samp(k)/M0*(1-b*Smean(k))))^(-beta);
%             Nf(gi)=(1+r)*Nf(gi);
            Neff(gi)=cycles(k)+Nf(gi)*(Neff(gi-1)/Nf_mod(gi-1))^((1-alpha(gi))/(1-alpha_mod(gi-1)));
            dam(gi)=(Neff(gi)/Nf(gi))^(1/(1-alpha(gi)));
            %recovery
%             dam(gi)=dam(gi-1)+(dam(gi)-dam(gi-1))*(1-r);
%             Dtot(gi)=dam(gi)+Ca;
%             Nf(gi)=Neff(gi)/(dam(gi)^(1-alpha(gi)));
            alpha_mod(gi)=alpha(gi);
            Nf_mod(gi)=Nf(gi);
        else
%             dam(gi)=dam(gi-1)/((1+r)^(1/(1-alpha(k))));
%             Dtot(gi)=dam(gi)+Ca;
            alpha_mod(gi)=alpha_mod(gi-1);
            Nf_mod(gi)=Nf_mod(gi-1);
            Neff(gi)=Neff(gi-1);
        end
        dam(gi)=(1-r)*dam(gi);
        Nf_mod(gi)=Nf_mod(gi-1)*(1+r);
    end
end


%Post-processing
dam_eod=dam(1:m:end);
% Dtot_eod=Dtot(1:m:end);
% % steps=(1:m*i)/9*n;
% steps=(1:i)*n_steps;
% days=steps/n;

end