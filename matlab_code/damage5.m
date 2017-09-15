%********************************************************************
%Damage Evolution function based on non-linear damage evolution model by
%Marko and Starkey
%Author: Shruti Motiwale
%Code Status - Results not verified
%Please update the code status and specify modifications
%********************************************************************

function [Dmech, Dtot] = damage5(matrix, p, n, b)
% matrix=cycle data
% p = number of days
% n = number of steps in a day
% b = exponent of cycle ratio

%Define variables
age=90*365;%days
%Assuming 0.6 damage due to aging in 90 years
r=0.1;
Ca=0.6/(age-20*365);%aging constant
a=3; beta=2.9357;
% sf=0.45*UTS; %Fatigue Limit
M0=14.1672;%MPa

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

Nf=(Samp/M0).^(-beta);
%Damage Evolution
for i=1:p
    for k=1:length(Samp)
        gi=m*(i-1)+k; %global index
        if gi==1
            dam(gi)=(cycles(k)/Nf(k))^b;
        else
            dam(gi)=dam(gi-1)+(cycles(k)/Nf(k))^b;
        end
    end
    Dmech(i)=dam(gi);
    if i<365*20
        Dtot(i)=Dmech(i);
    else
        Dtot(i)=dam(gi)+Ca*(i-20*365);
    end
end

