function [reduced_dam] = damage1(matrix, p, n, r)
% matrix=cycle data
% p = number of days
% n = number of steps
% r = recovery

%Define variables
a=0.9711;
beta=3.0126;
b=1;
UTS=3.9;%MPa
% sf=0.45*UTS; %Fatigue Limit
sf=0;
M0=14.7358;%MPa
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
mat=sortrows(matrix,-1);
Smax=mat(:,1);
Samp=mat(:,2);
Smean=mat(:,3);
% Approximation strategy: all cycles of a certain amplitude 
% in a day are applied at once
cycles=n*mat(:,4); 
m=length(Samp);


%Damage Evolution
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
            alpha_mod(gi)=alpha_mod(gi-1);
            Nf_mod(gi)=Nf_mod(gi-1)*(1+r);
            Neff(gi)=Neff(gi-1);
        end
%         dam(gi)=(1-r)*dam(gi);
    end
end


%Post-processing
reduced_dam=dam(1:m:end);
% % steps=(1:m*i)/9*n;
% steps=(1:i)*n_steps;
% days=steps/n;

end