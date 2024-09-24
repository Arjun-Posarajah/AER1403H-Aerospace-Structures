%%Question 3%%
t=0.25;%mm
d=5;%um
vf=0.66;
E1=125;
E2=7.5;
G12=6.2;
v12=0.26;
x=30;
%%%Part 1 stiffness matrix in the material coordinate system
S11=1/E1;
S22=1/E2;
S12=(-v12)/E1;
S66=1/G12;

Q11=S22/((S11*S22)-(S12.^2));
Q22=S11/((S11*S22)-(S12.^2));
Q12=(-S12)/((S11*S22)-(S12.^2));
Q66=G12;

Q=[Q11 Q12 0; Q12 Q22 0; 0 0 Q66]

%%%Part2 stiffness matrix for the same lamina rotated +30o from the global coordinate system
Qbar11= (Q11*((cosd(x)).^4))+(Q22*((sind(x).^4)))+ (2*(Q12+(2*Q66)))*((sind(x)).^2)*((cosd(x)).^2);
Qbar12= (Q11+Q22-(4*Q66))*((cosd(x)).^2)*((sind(x)).^2) + (Q12)*(((cosd(x)).^4)+((sind(x)).^4));
Qbar22=(Q11)*((sind(x)).^4)+(Q22)*((cosd(x)).^4)+(2*(Q12+(2*Q66)))*(((sind(x)).^2)*((cosd(x)).^2));
Qbar16=(Q11-Q12-(2*Q66))*((cosd(x)).^3)*(sind(x))-((Q22-Q12-(2*Q66))*(cosd(x))*(((sind(x)).^3)));
Qbar26=(Q11-Q12-(2*Q66))*(cosd(x))*((sind(x)).^3)-(Q22-Q12-(2*Q66))*((cosd(x)).^3)*(sind(x));
Qbar66= ((Q11+Q22-(2*Q12)-(2*Q66))*((cosd(x)).^2)*((sind(x)).^2))+((Q66)*(((cosd(x)).^4)+((sind(x)).^4)));

Qbar=[Qbar11 Qbar12 Qbar16; Qbar12 Qbar22 Qbar26 ; Qbar16 Qbar26 Qbar66]








