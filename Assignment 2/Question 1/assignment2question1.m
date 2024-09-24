%%Assignment 2 Question 1%%
%Variables
E1=135;
E2=10.5;
G12=6.6;
v12=0.24;
x=(-90:1:90);

%Stiffness Matrix
S11=1/E1;
S22=1/E2;
S12=(-v12)/E1;
S66=1/G12;

%Qmatrix
Q11=S22/((S11*S22)-(S12.^2));
Q22=S11/((S11*S22)-(S12.^2));
Q12=(-S12)/((S11*S22)-(S12.^2));
Q66=G12;

%Rotation
Qbar11=(Q11*((cosd(x)).^4))+(Q22*((sind(x).^4)))+ (2*(Q12+(2*Q66)))*((sind(x)).^2).*((cosd(x)).^2);

%Plot
hold on
plot(x,Qbar11)
title ('Qbar11 terms from theta -90deg to 90deg')
xlabel 'Theta (deg)'
ylabel 'Qbar11 (GPa)'
hold off







