%%%Project 1 Arjun Posarajah 
%%Optimal Design of Composite Sandwich Beams
%Student Number: 1004881737
%Last 7 digits: 4881737
%               tuvwxyz
pf= 1600+(30*48)
Ef= (40+81)*10^9
SigmaF= (200+(100*7))*10^6
pc= 20+(5*37)
SigmaC= (0.5+(6.5*((37/100)^1.5)))*10^6
Tc= (0.5+(4.5*((37/100)^1.5)))*10^6

%tbar=0:0.001:0.5;
cbar=0:0.001:0.5;


%CoreShear and Microbuckling 
tbar= Tc./(2*SigmaF.*cbar);
plot(cbar,tbar)
hold on

%%Elastic Indentation and Microbuckling 
fimplicit(@(cbar,tbar)((tbar.^2)+tbar+1)-(((pi^2)*Ef*(SigmaC^2))/(192*(SigmaF^3)*(cbar.^2))))

%Coreshear and Elastic Indentation 
fimplicit(@(cbar,tbar)((((pi^2)*Ef*(SigmaC^2))/(24*(Tc^3))).*(tbar.^3).*cbar)-(tbar.^2)-(2.*tbar)-1)
axis([0 0.03 0 0.3])
xlabel c/l
ylabel t/c
%hold off

%%Meshgrid system
% Create a 100x100 meshgrid
cbar = linspace(0, 0.3, 5000);
tbar = linspace(0, 0.3, 5000);
[X, Y] = meshgrid(cbar, tbar);

% Define three equations
cs= ((2*Tc)/SigmaF).*(X)+((2*Tc)/SigmaF).*X.*Y;
mb= 4*(X.^2).*Y+4.*(Y.^2).*(X.^2);
ei=(1/SigmaF).*(Y.*X).*(((((pi^2)*Ef*(SigmaC^2))/3).*X.*(1+Y)).^(1/3));

index = zeros(5000);

for i=1:5000
    for j=1:5000
        if cs(i,j)<mb(i,j) && cs(i,j)<ei(i,j)
            index(i,j)=1;
        elseif mb(i,j)<cs(i,j) && mb(i,j)<ei(i,j)
                index(i,j)=2;
        else 
                    index(i,j)=3;
                
          end 
     end
end 
contourf(X,Y,index)
        
%Trajectory of Optimal Design on the mechanism map    
yline([0.0359 0.116])

title('Failure Mechanism Map')
hold off


        
    