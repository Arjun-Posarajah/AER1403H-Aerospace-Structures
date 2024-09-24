%%Question 1%%
%Ply Angles
x=-90:1:90;
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

%%x
Qbar11= (Q11*((cosd(x)).^4))+(Q22*((sind(x).^4)))+ (2*(Q12+(2*Q66))).*((sind(x)).^2).*((cosd(x)).^2);
Qbar12= (Q11+Q22-(4*Q66))*((cosd(x)).^2).*((sind(x)).^2) + (Q12)*(((cosd(x)).^4)+((sind(x)).^4));
Qbar22=(Q11)*((sind(x)).^4)+(Q22)*((cosd(x)).^4)+(2*(Q12+(2*Q66)))*(((sind(x)).^2).*((cosd(x)).^2));
Qbar16=(Q11-Q12-(2*Q66))*((cosd(x)).^3).*(sind(x))-((Q22-Q12-(2*Q66))*(cosd(x)).*(((sind(x)).^3)));
Qbar26=(Q11-Q12-(2*Q66))*(cosd(x)).*((sind(x)).^3)-(Q22-Q12-(2*Q66))*((cosd(x)).^3).*(sind(x));
Qbar66= ((Q11+Q22-(2*Q12)-(2*Q66))*((cosd(x)).^2).*((sind(x)).^2))+((Q66)*(((cosd(x)).^4)+((sind(x)).^4)));

%%%Part3: Calculating A,B,D
A11=Qbar11*(t)+ Qbar11*(t)+Qbar11*(t)+ Qbar11*(t);
A12=Qbar12*(t)+ Qbar12*(t)+Qbar12*(t)+ Qbar12*(t);
A22=Qbar22*(t)+ Qbar22*(t)+Qbar22*(t)+ Qbar22*(t);
A16=Qbar16*(t)+ Qbar16*(t)+Qbar16*(t)+ Qbar16*(t);
A26=Qbar26*(t)+ Qbar26*(t)+Qbar26*(t)+ Qbar26*(t);
A66=Qbar66*(t)+ Qbar66*(t)+Qbar66*(t)+ Qbar66*(t);

B11= 0.5*((Qbar11*(3/64))+(Qbar11*(1/64))+(Qbar11*(1/64))+(Qbar11*(3/64)));
B12= 0.5*((Qbar12*(3/64))+(Qbar12*(1/64))+(Qbar12*(1/64))+(Qbar12*(3/64)));
B22= 0.5*((Qbar22*(3/64))+(Qbar22*(1/64))+(Qbar22*(1/64))+(Qbar22*(3/64)));
B16= 0.5*((Qbar16*(3/64))+(Qbar16*(1/64))+(Qbar16*(1/64))+(Qbar16*(3/64)));
B26= 0.5*((Qbar26*(3/64))+(Qbar26*(1/64))+(Qbar26*(1/64))+(Qbar26*(3/64)));
B66= 0.5*((Qbar66*(3/64))+(Qbar66*(1/64))+(Qbar66*(1/64))+(Qbar66*(3/64)));

D11=(1/3)*((Qbar11*(7/512))+(Qbar11*(1/512))+(Qbar11*(1/512))+(Qbar11*(7/512)));
D12=(1/3)*((Qbar12*(7/512))+(Qbar12*(1/512))+(Qbar12*(1/512))+(Qbar12*(7/512)));
D22=(1/3)*((Qbar22*(7/512))+(Qbar22*(1/512))+(Qbar22*(1/512))+(Qbar22*(7/512)));
D16=(1/3)*((Qbar16*(7/512))+(Qbar16*(1/512))+(Qbar16*(1/512))+(Qbar16*(7/512)));
D26=(1/3)*((Qbar26*(7/512))+(Qbar26*(1/512))+(Qbar26*(1/512))+(Qbar26*(7/512)));
D66=(1/3)*((Qbar66*(7/512))+(Qbar66*(1/512))+(Qbar66*(1/512))+(Qbar66*(7/512)));

%%%Part4 Laminate Stiffness Matrix
A= [A11 A12 A16; A12 A22 A26; A16 A26 A66];
B= [B11 B12 B16; B12 B22 B26; B16 B26 B66];
D= [D11 D12 D16; D12 D22 D26; D16 D26 D66];

Qij=[A11 A12 A16 B11 B12 B16; A12 A22 A26 B12 B22 B26;A16 A26 A66 B16 B26 B66; B11 B12 B16 D11 D12 D16 ;B12 B22 B26 D12 D22 D26; B16 B26 B66 D16 D26 D66];

%Part5 Solve Laminate Stiffness Matrix for Midplane Strains
F=[10 ;5; 0 ;10; 0 ;0];
Eo1=linsolve(Qij,F);
%Eo = [.001; .001; 0.001; 1; 1; 1].* Eo1 %convert the KN/m to N/m

%Part6 Calculate Strains at Lamina Mid-Planes in Global Coordinates
%%%%midplane is 0.0625mm

y3= 0.0625;

E3=[5.9+(y3*41.6103) ; (y3*-1.5344); -3.7+(y3*-3.0309)]


%Part7 Calculate Stress
Ss3= x3Qbar*E3

%Part8 Calculate Lamina Stresses in Material Coordinates
Tk3= [(cosd(x3))^2 (sind(x3))^2 2*(cosd(x3))*(sind(x3));(sind(x3))^2 (cosd(x3))^2 -2*(cosd(x3))*(sind(x3)); -1*(cosd(x3))*(sind(x3)) (cosd(x3))*(sind(x3)) ((sind(x3))^2)-((cosd(x3))^2)];
S3= Tk3*Ss3


%Part9 Maximum Stress Failure Criterion 

%%Compare actual values no math needed on paper

%Part10 Tsai Hill
sigmaLplus=850;
sigmaLminus=700;
sigmaTplus=40;
sigmaTminus=100;
TLT= 75;

Tsaihill=((S3(1,1).^2)/(sigmaLplus .^2))-(((S3(1,1))*(S3(2,1)))/(sigmaLplus .^2))+(((S3(2,1)).^2)/(sigmaTminus .^2))+(((S3(3,1)).^2)/(TLT.^2))

%Part11 Tsai Wu
F11= 1/((sigmaLplus)*(sigmaLminus))
F22= 1/((sigmaTplus)*(sigmaTminus))
F1= (1/(sigmaLplus))-(1/(sigmaLminus))
F2= (1/(sigmaTplus))-(1/(sigmaTminus))
F66= 1/(TLT.^2);
F12= (-1/2)*(sqrt(1/(sigmaLplus*sigmaLminus*sigmaTplus*sigmaTminus)))

TsaiWu= (F11*(S3(1,1).^2))+(F22*(S3(2,1).^2))+(F66*(S3(3,1).^2))+(F1*S3(1,1))+(F2*S3(2,1))+ 2*(F12*S3(1,1)*S3(2,1))









