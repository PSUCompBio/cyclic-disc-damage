%********************************************************************
%Damage Evolution function based on non-linear damage evolution model by
%Chaboche and Lemaitre
%Author: Shruti Motiwale
%Code Status - Works well as of 04/23/2016
%Please update the code status and specify modifications
%Use when initial damage is not zero. Specify initial state of damage using
%the argument d_in.
%********************************************************************

function [Dmech, Dtot, d_out] = damage6(d_in,matrix, p, n, b, r)
% matrix=cycle data
% p = number of days
% n = number of steps in a day
% r = recovery

%Define variables
age=90*365;%days
%Assuming 0.6 damage due to aging in 90 years
thr=20; %age threshold beyind which aging is assumed to begin
Ca=0.6/(age-thr*365);%aging constant
a=4.8898; beta=2.9357; M0=21.9419;
UTS=6.25;%MPa
% sf=0.45*UTS; %Fatigue Limit
sf=0;

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

Neff_p=d_in(1);
Nf_mod_p=d_in(2);
alpha_mod_p=d_in(3);
p_i=d_in(4);

%Damage Evolution
alpha(1)=1-a*max((Smax(1)-sf)/(UTS-Smax(1)),0);
if alpha(1)<1
    Nf(1)=(1/(1-alpha(1)))*abs((Samp(1)/(M0*(1-b*Smean(1)))))^(-beta);
    Neff(1)=cycles(1)+Nf(1)*(Neff_p/Nf_mod_p)^((1-alpha(1))/(1-alpha_mod_p));
    dam(1)=(Neff(1)/Nf(1))^(1/(1-alpha(1)));
    alpha_mod(1)=alpha(1);
    Nf_mod(1)=Nf(1);
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
dam(k)=(1-r)*dam(k);
Nf(k)=Nf(k)*(1+r);
Dmech(1)=dam(k);
if (p_i+1)<=365*thr
    Dtot(1)=Dmech(1);
else
    Dtot(1)=Dmech(1)+Ca*(p_i+1-thr*365);
end

for i=2:p
    for k=1:length(Samp)
        gi=m*(i-1)+k; %global index
        alpha(gi)=1-a*max((Smax(k)-sf)/(UTS-Smax(k)),0);
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
    Nf(gi)=Nf(gi-1)*(1+r);
    Dmech(i)=dam(gi);
    if (p_i+i)<=365*thr
        Dtot(i)=Dmech(i);
    else
        Dtot(i)=Dmech(i)+Ca*(p_i+i-thr*365);
    end
end

d_out=[Neff(end), Nf_mod(end), alpha_mod(end), p_i+i];
end