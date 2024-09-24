%%Question 1%%
%Ply Angles
x1=0;
x2=45;
x3=60;
t=0.125;%mm
E1=125;
E2=8.8;
G12=6.8;
v12=0.24;

%%%Part 1 stiffness matrix in the material coordinate system
S11=1/E1;
S22=1/E2;
S12=(-v12)/E1;
S66=1/G12;

Q11=S22/((S11*S22)-(S12.^2));
Q22=S11/((S11*S22)-(S12.^2));
Q12=(-S12)/((S11*S22)-(S12.^2));
Q66=G12;

Q= [Q11 Q12 0; Q12 Q22 0; 0 0 Q66];

%%%Part2 stiffness matrix for the same lamina rotated angles from the global coordinate system

%%x1
x1Qbar11= (Q11*((cosd(x1)).^4))+(Q22*((sind(x1).^4)))+ (2*(Q12+(2*Q66)))*((sind(x1)).^2)*((cosd(x1)).^2);
x1Qbar12= (Q11+Q22-(4*Q66))*((cosd(x1)).^2)*((sind(x1)).^2) + (Q12)*(((cosd(x1)).^4)+((sind(x1)).^4));
x1Qbar22=(Q11)*((sind(x1)).^4)+(Q22)*((cosd(x1)).^4)+(2*(Q12+(2*Q66)))*(((sind(x1)).^2)*((cosd(x1)).^2));
x1Qbar16=(Q11-Q12-(2*Q66))*((cosd(x1)).^3)*(sind(x1))-((Q22-Q12-(2*Q66))*(cosd(x1))*(((sind(x1)).^3)));
x1Qbar26=(Q11-Q12-(2*Q66))*(cosd(x1))*((sind(x1)).^3)-(Q22-Q12-(2*Q66))*((cosd(x1)).^3)*(sind(x1));
x1Qbar66= ((Q11+Q22-(2*Q12)-(2*Q66))*((cosd(x1)).^2)*((sind(x1)).^2))+((Q66)*(((cosd(x1)).^4)+((sind(x1)).^4)));

x1Qbar= [x1Qbar11 x1Qbar12 x1Qbar16; x1Qbar12 x1Qbar22 x1Qbar26; x1Qbar16 x1Qbar26 x1Qbar66];

%%x2
x2Qbar11= (Q11*((cosd(x2)).^4))+(Q22*((sind(x2).^4)))+ (2*(Q12+(2*Q66)))*((sind(x2)).^2)*((cosd(x2)).^2);
x2Qbar12= (Q11+Q22-(4*Q66))*((cosd(x2)).^2)*((sind(x2)).^2) + (Q12)*(((cosd(x2)).^4)+((sind(x2)).^4));
x2Qbar22=(Q11)*((sind(x2)).^4)+(Q22)*((cosd(x2)).^4)+(2*(Q12+(2*Q66)))*(((sind(x2)).^2)*((cosd(x2)).^2));
x2Qbar16=(Q11-Q12-(2*Q66))*((cosd(x2)).^3)*(sind(x2))-((Q22-Q12-(2*Q66))*(cosd(x2))*(((sind(x2)).^3)));
x2Qbar26=(Q11-Q12-(2*Q66))*(cosd(x2))*((sind(x2)).^3)-(Q22-Q12-(2*Q66))*((cosd(x2)).^3)*(sind(x2));
x2Qbar66= ((Q11+Q22-(2*Q12)-(2*Q66))*((cosd(x2)).^2)*((sind(x2)).^2))+((Q66)*(((cosd(x2)).^4)+((sind(x2)).^4)));

x2Qbar= [x2Qbar11 x2Qbar12 x2Qbar16; x2Qbar12 x2Qbar22 x2Qbar26; x2Qbar16 x2Qbar26 x2Qbar66];

%%x3
x3Qbar11= (Q11*((cosd(x3)).^4))+(Q22*((sind(x3).^4)))+ (2*(Q12+(2*Q66)))*((sind(x3)).^2)*((cosd(x3)).^2);
x3Qbar12= (Q11+Q22-(4*Q66))*((cosd(x3)).^2)*((sind(x3)).^2) + (Q12)*(((cosd(x3)).^4)+((sind(x3)).^4));
x3Qbar22=(Q11)*((sind(x3)).^4)+(Q22)*((cosd(x3)).^4)+(2*(Q12+(2*Q66)))*(((sind(x3)).^2)*((cosd(x3)).^2));
x3Qbar16=(Q11-Q12-(2*Q66))*((cosd(x3)).^3)*(sind(x3))-((Q22-Q12-(2*Q66))*(cosd(x3))*(((sind(x3)).^3)));
x3Qbar26=(Q11-Q12-(2*Q66))*(cosd(x3))*((sind(x3)).^3)-(Q22-Q12-(2*Q66))*((cosd(x3)).^3)*(sind(x3));
x3Qbar66= ((Q11+Q22-(2*Q12)-(2*Q66))*((cosd(x3)).^2)*((sind(x3)).^2))+((Q66)*(((cosd(x3)).^4)+((sind(x3)).^4)));

x3Qbar= [x3Qbar11 x3Qbar12 x3Qbar16; x3Qbar12 x3Qbar22 x3Qbar26; x3Qbar16 x3Qbar26 x3Qbar66];

%%%Part3: Calculating A,B,D
%Need to input z manually for spacing changes/ply changes

A11=x1Qbar11*(t)+ x2Qbar11*(t)+x3Qbar11*(t);
A12=x1Qbar12*(t)+ x2Qbar12*(t)+x3Qbar12*(t);
A22=x1Qbar22*(t)+ x2Qbar22*(t)+x3Qbar22*(t);
A16=x1Qbar16*(t)+ x2Qbar16*(t)+x3Qbar16*(t);
A26=x1Qbar26*(t)+ x2Qbar26*(t)+x3Qbar26*(t);
A66=x1Qbar66*(t)+ x2Qbar66*(t)+x3Qbar66*(t);

B11= 0.5*((x1Qbar11*(-1/32))+(x2Qbar11*(0))+(x3Qbar11*(1/32)));
B12= 0.5*((x1Qbar12*(-1/32))+(x2Qbar12*(0))+(x3Qbar12*(1/32)));
B22= 0.5*((x1Qbar22*(-1/32))+(x2Qbar22*(0))+(x3Qbar22*(1/32)));
B16= 0.5*((x1Qbar16*(-1/32))+(x2Qbar16*(0))+(x3Qbar16*(1/32)));
B26= 0.5*((x1Qbar26*(-1/32))+(x2Qbar26*(0))+(x3Qbar26*(1/32)));
B66= 0.5*((x1Qbar66*(-1/32))+(x2Qbar66*(0))+(x3Qbar66*(1/32)));

D11=(1/3)*((x1Qbar11*(13/2048))+(x2Qbar11*(1/2048))+(x3Qbar11*(13/2048)));
D12=(1/3)*((x1Qbar12*(13/2048))+(x2Qbar12*(1/2048))+(x3Qbar12*(13/2048)));
D22=(1/3)*((x1Qbar22*(13/2048))+(x2Qbar22*(1/2048))+(x3Qbar22*(13/2048)));
D16=(1/3)*((x1Qbar16*(13/2048))+(x2Qbar16*(1/2048))+(x3Qbar16*(13/2048)));
D26=(1/3)*((x1Qbar26*(13/2048))+(x2Qbar26*(1/2048))+(x3Qbar26*(13/2048)));
D66=(1/3)*((x1Qbar66*(13/2048))+(x2Qbar66*(1/2048))+(x3Qbar66*(13/2048)));

%%%Part4 Laminate Stiffness Matrix
A= [A11 A12 A16; A12 A22 A26; A16 A26 A66];
B= [B11 B12 B16; B12 B22 B26; B16 B26 B66];
D= [D11 D12 D16; D12 D22 D26; D16 D26 D66];

Qij=[A11 A12 A16 B11 B12 B16; A12 A22 A26 B12 B22 B26;A16 A26 A66 B16 B26 B66; B11 B12 B16 D11 D12 D16 ;B12 B22 B26 D12 D22 D26; B16 B26 B66 D16 D26 D66];

%Part5 Solve Laminate Stiffness Matrix for Midplane Strains
F=[6 ;-2; 0 ;2 ; 2 ;0];
Eo=linsolve(Qij,F);

%Part6 Calculate Strains at Lamina Mid-Planes in Global Coordinates
y1=-0.125;
y2=0;
y3=0.125;

E1= [Eo(1,1)+(y1*(Eo(4,1)));Eo(2,1)+(y1*(Eo(5,1))); Eo(3,1)+(y1*(Eo(6,1)))];
E2= [Eo(1,1)+(y2*(Eo(4,1)));Eo(2,1)+(y2*(Eo(5,1))); Eo(3,1)+(y2*(Eo(6,1)))];
E3= [Eo(1,1)+(y3*(Eo(4,1)));Eo(2,1)+(y3*(Eo(5,1))); Eo(3,1)+(y3*(Eo(6,1)))];

%Part7 Calculate Stress
Ss1=x1Qbar*E1
Ss2=x2Qbar*E2
Ss3=x3Qbar*E3

%Part8 Calculate Lamina Stresses in Material Coordinates
Tk1= [(cosd(x1))^2 (sind(x1))^2 2*(cosd(x1))*(sind(x1));(sind(x1))^2 (cosd(x1))^2 -2*(cosd(x1))*(sind(x1)); -1*(cosd(x1))*(sind(x1)) (cosd(x1))*(sind(x1)) ((cosd(x1)^2)-((sind(x1))^2))];
Tk2= [(cosd(x2))^2 (sind(x2))^2 2*(cosd(x2))*(sind(x2));(sind(x2))^2 (cosd(x2))^2 -2*(cosd(x2))*(sind(x2)); -1*(cosd(x2))*(sind(x2)) (cosd(x2))*(sind(x2)) ((cosd(x2)^2)-((sind(x2))^2))];
Tk3= [(cosd(x3))^2 (sind(x3))^2 2*(cosd(x3))*(sind(x3));(sind(x3))^2 (cosd(x3))^2 -2*(cosd(x3))*(sind(x3)); -1*(cosd(x3))*(sind(x3)) (cosd(x3))*(sind(x3)) ((cosd(x3)^2)-((sind(x3))^2))];

S1=Tk1*Ss1
S2=Tk2*Ss2
S3= Tk3*Ss3

%Part9 Tsai Hill
sigmaLplus=900;
sigmaLminus=750;
sigmaTplus=50;
sigmaTminus=180;
TLT= 80;

ply1Tsaihill=((S1(1,1).^2)/(sigmaLminus .^2))-(((S1(1,1))*(S1(2,1)))/(sigmaLminus .^2))+(((S1(2,1)).^2)/(sigmaTminus .^2))+(((S1(3,1)).^2)/(TLT.^2))
ply2Tsaihill=((S2(1,1).^2)/(sigmaLminus .^2))-(((S2(1,1))*(S2(2,1)))/(sigmaLminus .^2))+(((S2(2,1)).^2)/(sigmaTplus .^2))+(((S2(3,1)).^2)/(TLT.^2))
ply3Tsaihill=((S3(1,1).^2)/(sigmaLplus .^2))-(((S3(1,1))*(S3(2,1)))/(sigmaLplus .^2))+(((S3(2,1)).^2)/(sigmaTplus .^2))+(((S3(3,1)).^2)/(TLT.^2))

%Part10 Tsai Wu
F11= 1/((sigmaLplus)*(sigmaLminus))
F22= 1/((sigmaTplus)*(sigmaTminus))
F1= (1/(sigmaLplus))-(1/(sigmaLminus))
F2= (1/(sigmaTplus))-(1/(sigmaTminus))
F66= 1/(TLT.^2)
F12= (-1/2)*(sqrt(F11*F22))

ply1TsaiWu= (F11*(S1(1,1).^2))+(F22*(S1(2,1).^2))+(F66*(S1(3,1).^2))+(F1*S1(1,1))+(F2*S1(2,1))+ 2*(F12*S1(1,1)*S1(2,1))
ply2TsaiWu= (F11*(S2(1,1).^2))+(F22*(S2(2,1).^2))+(F66*(S2(3,1).^2))+(F1*S2(1,1))+(F2*S2(2,1))+ 2*(F12*S2(1,1)*S2(2,1))
ply3TsaiWu= (F11*(S3(1,1).^2))+(F22*(S3(2,1).^2))+(F66*(S3(3,1).^2))+(F1*S3(1,1))+(F2*S3(2,1))+ 2*(F12*S3(1,1)*S3(2,1))
