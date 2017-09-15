%********************************************************************
%Damage Evolution function based on non-linear damage evolution model by
%Chaboche and Lemaitre
%Author: Shruti Motiwale
%Code Status - Works well as of 04/23/2016
%Please update the code status and specify modifications
%Includes threshold for recovery parameter
%********************************************************************

function [Dmech, Dtot, d_out] = damage8(matrix, p, n, b, r, r1)
% matrix=cycle data
% p = number of days
% n = number of steps in a day
% r = recovery

%Define variables
age=90*365;%days
%Assuming 0.6 damage due to aging in 90 years
thr_age=20; %age threshold beyond which aging is assumed to begin
Ca=0.6/(age-thr_age*365);%aging constant
thr_dam=0.25;
thr_dam_day=thr_dam/age;
a=4.8898; beta=2.9357; M0=21.9419;
UTS=6.25;%MPa
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
alpha=zeros(m*p,1); Neff=zeros(m*p,1); alpha_mod=zeros(m*p,1); dam=zeros(m*p,1); Nf=zeros(m*p,1); Nf_mod=zeros(m*p,1); 
Dmech=zeros(p,1); Dtot=zeros(p,1); 
%Damage Evolution
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
if dam_day<thr_dam_day
    r=r1;
end
dam(k)=(1-r)*dam(k);
Nf(k)=Neff(k)/(dam(k)^(1-alpha(k)));
Dtot(1)=dam(k)+Ca;
Dmech(1)=dam(k);

for i=2:p
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
    Dmech(i)=dam(gi);
    dam_day=Dmech(i)-Dmech(i-1);
    if dam_day<thr_dam_day
        r=r1;
    end
    Dmech(i)=Dmech(i-1)+dam_day*10^(-r);
    Nf(gi)=Neff(gi)/(dam(gi)^(1-alpha(gi)));
    if i<=365*thr_age
        Dtot(i)=Dmech(i);
    else
        Dtot(i)=Dmech(i)+Ca*(i-thr_age*365);
    end
end
d_out=[Neff(end), Nf_mod(end), alpha_mod(end), i];
end