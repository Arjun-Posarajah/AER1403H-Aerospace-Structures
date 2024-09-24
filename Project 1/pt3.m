pf= 1600+(30*48)
Ef= (40+81)*10^9
SigmaF= (200+(100*7))*10^6
pc= 20+(5*37)
SigmaC= (0.5+(6.5*((37/100)^1.5)))*10^6
Tc= (0.5+(4.5*((37/100)^1.5)))*10^6

%M=(tbar).*(cbar)+ (pc/pf).*(cbar);

%%Meshgrid system

% Create a 100x100 meshgrid
cbar = linspace(0, 0.3, 100);
tbar = linspace(0, 0.3, 100);
[X, Y] = meshgrid(cbar, tbar);

% Define three equations
cs= ((2*Tc)/SigmaF).*(X)+((2*Tc)/SigmaF).*X.*Y;
mb= 4*(X.^2).*Y+4.*(Y.^2).*(X.^2);
ei=(1/SigmaF).*(Y.*X).*(((((pi^2)*Ef*(SigmaC^2))/3).*X.*(1+Y)).^(1/3));
M=(Y).*(X)+ (pc/pf).*(X);

ind = zeros(100);

for i=1:100
    for j=1:100
        if cs(i,j)<mb(i,j) && cs(i,j)<ei(i,j)
            ind(i,j)=cs(i,j);
        elseif mb(i,j)<cs(i,j) && mb(i,j)<ei(i,j)
                ind(i,j)=mb(i,j);
        else 
                    ind(i,j)=ei(i,j);
                
          end 
     end
end 

plot(ind,M)