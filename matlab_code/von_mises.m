close all
von_noach=zeros(247, 1);
von_ach=zeros(247, 1);
for i=1:247
    von_noach(i)=(1/sqrt(2))*(sqrt((2*((walk_noACH_3dStress(i, 2))^2))+6*(((walk_noACH_3dStress(i, 1))^2)+((walk_noACH_3dStress(i, 3))^2))));
    von_ach(i)=(1/sqrt(2))*(sqrt((2*((walk_ACH_3dStress(i, 2))^2))+6*(((walk_ACH_3dStress(i, 1))^2)+((walk_ACH_3dStress(i, 3))^2))));
end

figure
plot(time, von_noach);
hold on
plot(time, walk_noACH_3dStress(:, 2));
figure
plot(time, von_ach);
hold on
plot(time, walk_ACH_3dStress(:, 2));
hold off


