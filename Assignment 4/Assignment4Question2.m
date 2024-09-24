L=1500;%mm
b=35;%mm
t=3.5;
c=30;%mm
d=t+c;%mm
E=75000; %MPa
SigmaF=-600; %MPa Composite Compressive Strength
SigmaC= -4.5; %MPa Polymer Compressive Strength
TaoC= 2.6; %MPa Polymer Shear Strength


%Elastic Indentation
Pelasticindentation= b*t*((((pi^2)*d*E*(SigmaC^2))/(3*L))^(1/3))

%Core Shear
Pcoreshear= 2*b*d*TaoC

%Face Microbukling
Pmicrobuckling= (4*b*d*t*SigmaF)/L

