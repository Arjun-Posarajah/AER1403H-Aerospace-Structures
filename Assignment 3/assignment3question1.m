%%Question 1%%
%Ply Angles
x1=60;
x2=45;
x3=0;
x4=30;
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

Q11=S22/((S11*S22)-(S12.^2))
Q22=S11/((S11*S22)-(S12.^2))
Q12=(-S12)/((S11*S22)-(S12.^2))
Q66=G12

%%%Part2 stiffness matrix for the same lamina rotated angles from the global coordinate system

%%x1
x1Qbar11= (Q11*((cosd(x1)).^4))+(Q22*((sind(x1).^4)))+ (2*(Q12+(2*Q66)))*((sind(x1)).^2)*((cosd(x1)).^2);
x1Qbar12= (Q11+Q22-(4*Q66))*((cosd(x1)).^2)*((sind(x1)).^2) + (Q12)*(((cosd(x1)).^4)+((sind(x1)).^4));
x1Qbar22=(Q11)*((sind(x1)).^4)+(Q22)*((cosd(x1)).^4)+(2*(Q12+(2*Q66)))*(((sind(x1)).^2)*((cosd(x1)).^2));
x1Qbar16=(Q11-Q12-(2*Q66))*((cosd(x1)).^3)*(sind(x1))-((Q22-Q12-(2*Q66))*(cosd(x1))*(((sind(x1)).^3)));
x1Qbar26=(Q11-Q12-(2*Q66))*(cosd(x1))*((sind(x1)).^3)-(Q22-Q12-(2*Q66))*((cosd(x1)).^3)*(sind(x1));
x1Qbar66= ((Q11+Q22-(2*Q12)-(2*Q66))*((cosd(x1)).^2)*((sind(x1)).^2))+((Q66)*(((cosd(x1)).^4)+((sind(x1)).^4)));

x1Qbar= [x1Qbar11 x1Qbar12 x1Qbar16; x1Qbar12 x1Qbar22 x1Qbar26; x1Qbar16 x1Qbar26 x1Qbar66]


%%x2
x2Qbar11= (Q11*((cosd(x2)).^4))+(Q22*((sind(x2).^4)))+ (2*(Q12+(2*Q66)))*((sind(x2)).^2)*((cosd(x2)).^2);
x2Qbar12= (Q11+Q22-(4*Q66))*((cosd(x2)).^2)*((sind(x2)).^2) + (Q12)*(((cosd(x2)).^4)+((sind(x2)).^4));
x2Qbar22=(Q11)*((sind(x2)).^4)+(Q22)*((cosd(x2)).^4)+(2*(Q12+(2*Q66)))*(((sind(x2)).^2)*((cosd(x2)).^2));
x2Qbar16=(Q11-Q12-(2*Q66))*((cosd(x2)).^3)*(sind(x2))-((Q22-Q12-(2*Q66))*(cosd(x2))*(((sind(x2)).^3)));
x2Qbar26=(Q11-Q12-(2*Q66))*(cosd(x2))*((sind(x2)).^3)-(Q22-Q12-(2*Q66))*((cosd(x2)).^3)*(sind(x2));
x2Qbar66= ((Q11+Q22-(2*Q12)-(2*Q66))*((cosd(x2)).^2)*((sind(x2)).^2))+((Q66)*(((cosd(x2)).^4)+((sind(x2)).^4)));

x2Qbar= [x2Qbar11 x2Qbar12 x2Qbar16; x2Qbar12 x2Qbar22 x2Qbar26; x2Qbar16 x2Qbar26 x2Qbar66]

%%x3
x3Qbar11= (Q11*((cosd(x3)).^4))+(Q22*((sind(x3).^4)))+ (2*(Q12+(2*Q66)))*((sind(x3)).^2)*((cosd(x3)).^2);
x3Qbar12= (Q11+Q22-(4*Q66))*((cosd(x3)).^2)*((sind(x3)).^2) + (Q12)*(((cosd(x3)).^4)+((sind(x3)).^4));
x3Qbar22=(Q11)*((sind(x3)).^4)+(Q22)*((cosd(x3)).^4)+(2*(Q12+(2*Q66)))*(((sind(x3)).^2)*((cosd(x3)).^2));
x3Qbar16=(Q11-Q12-(2*Q66))*((cosd(x3)).^3)*(sind(x3))-((Q22-Q12-(2*Q66))*(cosd(x3))*(((sind(x3)).^3)));
x3Qbar26=(Q11-Q12-(2*Q66))*(cosd(x3))*((sind(x3)).^3)-(Q22-Q12-(2*Q66))*((cosd(x3)).^3)*(sind(x3));
x3Qbar66= ((Q11+Q22-(2*Q12)-(2*Q66))*((cosd(x3)).^2)*((sind(x3)).^2))+((Q66)*(((cosd(x3)).^4)+((sind(x3)).^4)));

x3Qbar= [x3Qbar11 x3Qbar12 x3Qbar16; x3Qbar12 x3Qbar22 x3Qbar26; x3Qbar16 x3Qbar26 x3Qbar66]

%%x4
x4Qbar11= (Q11*((cosd(x4)).^4))+(Q22*((sind(x4).^4)))+ (2*(Q12+(2*Q66)))*((sind(x4)).^2)*((cosd(x4)).^2);
x4Qbar12= (Q11+Q22-(4*Q66))*((cosd(x4)).^2)*((sind(x4)).^2) + (Q12)*(((cosd(x4)).^4)+((sind(x4)).^4));
x4Qbar22=(Q11)*((sind(x4)).^4)+(Q22)*((cosd(x4)).^4)+(2*(Q12+(2*Q66)))*(((sind(x4)).^2)*((cosd(x4)).^2));
x4Qbar16=(Q11-Q12-(2*Q66))*((cosd(x4)).^3)*(sind(x4))-((Q22-Q12-(2*Q66))*(cosd(x4))*(((sind(x4)).^3)));
x4Qbar26=(Q11-Q12-(2*Q66))*(cosd(x4))*((sind(x4)).^3)-(Q22-Q12-(2*Q66))*((cosd(x4)).^3)*(sind(x4));
x4Qbar66= ((Q11+Q22-(2*Q12)-(2*Q66))*((cosd(x4)).^2)*((sind(x4)).^2))+((Q66)*(((cosd(x4)).^4)+((sind(x4)).^4)));

x4Qbar= [x4Qbar11 x4Qbar12 x4Qbar16; x4Qbar12 x4Qbar22 x4Qbar26; x4Qbar16 x4Qbar26 x4Qbar66]

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

%Part5 Solve Laminate Stiffness Matrix for Midplane Strains
F=[15; 20; 0; 10; 5; 0];
Eo1=linsolve(Qij,F);
Eo= [.001; .001; 0.001; 1; 1; 1].* Eo1


%Part6 Calculate Strains at Lamina Mid-Planes in Global Coordinates
%%%%midplane is 0.0625mm

y3= 0.0625;

E3=[-7.1+(y3*85.234) ; -3.5+(y3*-45.0943); 9.3+(y3*-120.3244)]


%Part7 Calculate Stress
Ss3= x3Qbar*E3

%Part8 Calculate Lamina Stresses in Material Coordinates
Tk3= [(cosd(x3))^2 (sind(x3))^2 2*(cosd(x3))*(sind(x3));(sind(x3))^2 (cosd(x3))^2 -2*(cosd(x3))*(sind(x3)); -1*(cosd(x3))*(sind(x3)) (cosd(x3))*(sind(x3)) ((sind(x3))^2)-((cosd(x3))^2)];
S3= Tk3*Ss3


%Part9 Maximum Stress Failure Criterion 

%%Compare actual values no math needed; on paper

%Part10 Tsai Hill
sigmaLplus=850;
sigmaLminus=700;
sigmaTplus=40;
sigmaTminus=100;
TLT= 75;

Tsaihill=((S3(1,1).^2)/(sigmaLplus .^2))-(((S3(1,1))*(S3(2,1)))/(sigmaLplus .^2))+(((S3(2,1)).^2)/(sigmaTminus .^2))+(((S3(3,1)).^2)/(TLT.^2))

%Part11 Tsai Wu
F11= 1/((sigmaLplus)*(sigmaLminus));
F22= 1/((sigmaTplus)*(sigmaTminus));
F1= (1/(sigmaLplus))-(1/(sigmaLminus));
F2= (1/(sigmaTplus))-(1/(sigmaTminus));
F66= 1/(TLT.^2);
F12= (-1/2)*(sqrt(1/(sigmaLplus*sigmaLminus*sigmaTplus*sigmaTminus)));

TsaiWu= (F11*(S3(1,1).^2))+(F22*(S3(2,1).^2))+(F66*(S3(3,1).^2))+(F1*S3(1,1))+(F2*S3(2,1))+ 2*(F12*S3(1,1)*S3(2,1))












