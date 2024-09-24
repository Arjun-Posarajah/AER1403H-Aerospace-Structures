%Lamina Stresses in material Coordinates
S=800
%Angle range
x=[0:1:90]
%Transform Matrix
T=[(cosd(x)).^2 (sind(x)).^2 2.*(cosd(x)).*(sind(x));(sind(x)).^2 (cosd(x)).^2 -2.*(cosd(x)).*(sind(x)); -1*(cosd(x)).*(sind(x)) (cosd(x)).*(sind(x)) ((sind(x)).^2)-((cosd(x)).^2)];

%Solving for Sxx
%Sxx=linsolve(T,S)
 Sxx=S./T

 plot(x,Sxx)

%Tsai Wu Criterion 
sigmaLplus=1500;
sigmaLminus=-1200;
sigmaTplus=40;
sigmaTminus=-190;
TLT=75;

F11= 1/((sigmaLplus)*(sigmaLminus));
F22= 1/((sigmaTplus)*(sigmaTminus));
F1= (1/(sigmaLplus))-(1/(sigmaLminus));
F2= (1/(sigmaTplus))-(1/(sigmaTminus));
F66= 1/(TLT.^2);
F12= (-1/2)*(sqrt(1/(sigmaLplus*sigmaLminus*sigmaTplus*sigmaTminus)));

TsaiWu= (F11*(S.^2))+(F22*(0.^2))+(F66*(0.^2))+(F1*S)+(F2*0)+ 2*(F12*S*0)