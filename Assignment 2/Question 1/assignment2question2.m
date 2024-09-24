%%Question 2%%
%Student Number 1004881737->Last 8 digits: 04881737%
%Ply Angles
x1=(04-50)*1.8; %-82.8deg
x2=(88-50)*1.8; %68.4deg
x3=(17-50)*1.8; %-59.4deg
x4=(37-50)*1.8; %-23.4deg
t=0.125;%mm
E1=131;
E2=9.8;
G12=5.8;
v12=0.22;

%%%Part 1 stiffness matrix in the material coordinate system
S11=1/E1;
S22=1/E2;
S12=(-v12)/E1;
S66=1/G12;

Q11=S22/((S11*S22)-(S12.^2));
Q22=S11/((S11*S22)-(S12.^2));
Q12=(-S12)/((S11*S22)-(S12.^2));
Q66=G12;

%%%Part2 stiffness matrix for the same lamina rotated angles from the global coordinate system

%%x1
x1Qbar11= (Q11*((cosd(x1)).^4))+(Q22*((sind(x1).^4)))+ (2*(Q12+(2*Q66)))*((sind(x1)).^2)*((cosd(x1)).^2);
x1Qbar12= (Q11+Q22-(4*Q66))*((cosd(x1)).^2)*((sind(x1)).^2) + (Q12)*(((cosd(x1)).^4)+((sind(x1)).^4));
x1Qbar22=(Q11)*((sind(x1)).^4)+(Q22)*((cosd(x1)).^4)+(2*(Q12+(2*Q66)))*(((sind(x1)).^2)*((cosd(x1)).^2));
x1Qbar16=(Q11-Q12-(2*Q66))*((cosd(x1)).^3)*(sind(x1))-((Q22-Q12-(2*Q66))*(cosd(x1))*(((sind(x1)).^3)));
x1Qbar26=(Q11-Q12-(2*Q66))*(cosd(x1))*((sind(x1)).^3)-(Q22-Q12-(2*Q66))*((cosd(x1)).^3)*(sind(x1));
x1Qbar66= ((Q11+Q22-(2*Q12)-(2*Q66))*((cosd(x1)).^2)*((sind(x1)).^2))+((Q66)*(((cosd(x1)).^4)+((sind(x1)).^4)));

%%x2
x2Qbar11= (Q11*((cosd(x2)).^4))+(Q22*((sind(x2).^4)))+ (2*(Q12+(2*Q66)))*((sind(x2)).^2)*((cosd(x2)).^2);
x2Qbar12= (Q11+Q22-(4*Q66))*((cosd(x2)).^2)*((sind(x2)).^2) + (Q12)*(((cosd(x2)).^4)+((sind(x2)).^4));
x2Qbar22=(Q11)*((sind(x2)).^4)+(Q22)*((cosd(x2)).^4)+(2*(Q12+(2*Q66)))*(((sind(x2)).^2)*((cosd(x2)).^2));
x2Qbar16=(Q11-Q12-(2*Q66))*((cosd(x2)).^3)*(sind(x2))-((Q22-Q12-(2*Q66))*(cosd(x2))*(((sind(x2)).^3)));
x2Qbar26=(Q11-Q12-(2*Q66))*(cosd(x2))*((sind(x2)).^3)-(Q22-Q12-(2*Q66))*((cosd(x2)).^3)*(sind(x2));
x2Qbar66= ((Q11+Q22-(2*Q12)-(2*Q66))*((cosd(x2)).^2)*((sind(x2)).^2))+((Q66)*(((cosd(x2)).^4)+((sind(x2)).^4)));

%%x3
x3Qbar11= (Q11*((cosd(x3)).^4))+(Q22*((sind(x3).^4)))+ (2*(Q12+(2*Q66)))*((sind(x3)).^2)*((cosd(x3)).^2);
x3Qbar12= (Q11+Q22-(4*Q66))*((cosd(x3)).^2)*((sind(x3)).^2) + (Q12)*(((cosd(x3)).^4)+((sind(x3)).^4));
x3Qbar22=(Q11)*((sind(x3)).^4)+(Q22)*((cosd(x3)).^4)+(2*(Q12+(2*Q66)))*(((sind(x3)).^2)*((cosd(x3)).^2));
x3Qbar16=(Q11-Q12-(2*Q66))*((cosd(x3)).^3)*(sind(x3))-((Q22-Q12-(2*Q66))*(cosd(x3))*(((sind(x3)).^3)));
x3Qbar26=(Q11-Q12-(2*Q66))*(cosd(x3))*((sind(x3)).^3)-(Q22-Q12-(2*Q66))*((cosd(x3)).^3)*(sind(x3));
x3Qbar66= ((Q11+Q22-(2*Q12)-(2*Q66))*((cosd(x3)).^2)*((sind(x3)).^2))+((Q66)*(((cosd(x3)).^4)+((sind(x3)).^4)));

%%x4
x4Qbar11= (Q11*((cosd(x4)).^4))+(Q22*((sind(x4).^4)))+ (2*(Q12+(2*Q66)))*((sind(x4)).^2)*((cosd(x4)).^2);
x4Qbar12= (Q11+Q22-(4*Q66))*((cosd(x4)).^2)*((sind(x4)).^2) + (Q12)*(((cosd(x4)).^4)+((sind(x4)).^4));
x4Qbar22=(Q11)*((sind(x4)).^4)+(Q22)*((cosd(x4)).^4)+(2*(Q12+(2*Q66)))*(((sind(x4)).^2)*((cosd(x4)).^2));
x4Qbar16=(Q11-Q12-(2*Q66))*((cosd(x4)).^3)*(sind(x4))-((Q22-Q12-(2*Q66))*(cosd(x4))*(((sind(x4)).^3)));
x4Qbar26=(Q11-Q12-(2*Q66))*(cosd(x4))*((sind(x4)).^3)-(Q22-Q12-(2*Q66))*((cosd(x4)).^3)*(sind(x4));
x4Qbar66= ((Q11+Q22-(2*Q12)-(2*Q66))*((cosd(x4)).^2)*((sind(x4)).^2))+((Q66)*(((cosd(x4)).^4)+((sind(x4)).^4)));


%%%Part3: Calculating A,B,D
A11=x1Qbar11*(t)+ x2Qbar11*(t)+x3Qbar11*(t)+ x4Qbar11*(t);
A12=x1Qbar12*(t)+ x2Qbar12*(t)+x3Qbar12*(t)+ x4Qbar12*(t);
A22=x1Qbar22*(t)+ x2Qbar22*(t)+x3Qbar22*(t)+ x4Qbar22*(t);
A16=x1Qbar16*(t)+ x2Qbar16*(t)+x3Qbar16*(t)+ x4Qbar16*(t);
A26=x1Qbar26*(t)+ x2Qbar26*(t)+x3Qbar26*(t)+ x4Qbar26*(t);
A66=x1Qbar66*(t)+ x2Qbar66*(t)+x3Qbar66*(t)+ x4Qbar66*(t);

B11= 0.5*((x1Qbar11*(3/64))+(x2Qbar11*(1/64))+(x3Qbar11*(1/64))+(x4Qbar11*(3/64)));
B12= 0.5*((x1Qbar12*(3/64))+(x2Qbar12*(1/64))+(x3Qbar12*(1/64))+(x4Qbar12*(3/64)));
B22= 0.5*((x1Qbar22*(3/64))+(x2Qbar22*(1/64))+(x3Qbar22*(1/64))+(x4Qbar22*(3/64)));
B16= 0.5*((x1Qbar16*(3/64))+(x2Qbar16*(1/64))+(x3Qbar16*(1/64))+(x4Qbar16*(3/64)));
B26= 0.5*((x1Qbar26*(3/64))+(x2Qbar26*(1/64))+(x3Qbar26*(1/64))+(x4Qbar26*(3/64)));
B66= 0.5*((x1Qbar66*(3/64))+(x2Qbar66*(1/64))+(x3Qbar66*(1/64))+(x4Qbar66*(3/64)));

D11=(1/3)*((x1Qbar11*(7/512))+(x2Qbar11*(1/512))+(x3Qbar11*(1/512))+(x4Qbar11*(7/512)));
D12=(1/3)*((x1Qbar12*(7/512))+(x2Qbar12*(1/512))+(x3Qbar12*(1/512))+(x4Qbar12*(7/512)));
D22=(1/3)*((x1Qbar22*(7/512))+(x2Qbar22*(1/512))+(x3Qbar22*(1/512))+(x4Qbar22*(7/512)));
D16=(1/3)*((x1Qbar16*(7/512))+(x2Qbar16*(1/512))+(x3Qbar16*(1/512))+(x4Qbar16*(7/512)));
D26=(1/3)*((x1Qbar26*(7/512))+(x2Qbar26*(1/512))+(x3Qbar26*(1/512))+(x4Qbar26*(7/512)));
D66=(1/3)*((x1Qbar66*(7/512))+(x2Qbar66*(1/512))+(x3Qbar66*(1/512))+(x4Qbar66*(7/512)));

%%%Part4 Laminate Stiffness Matrix
A= [A11 A12 A16; A12 A22 A26; A16 A26 A66]
B= [B11 B12 B16; B12 B22 B26; B16 B26 B66]
D= [D11 D12 D16; D12 D22 D26; D16 D26 D66]

Qij=[A11 A12 A16 B11 B12 B16; A12 A22 A26 B12 B22 B26;A16 A26 A66 B16 B26 B66; B11 B12 B16 D11 D12 D16 ;B12 B22 B26 D12 D22 D26; B16 B26 B66 D16 D26 D66]




















