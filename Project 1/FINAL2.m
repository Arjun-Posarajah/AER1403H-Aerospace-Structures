pf= 1600+(30*48);
Ef= (40+81)*10^9;
SigmaF= (200+(100*7))*10^6;
pc= 20+(5*37);
SigmaC= (0.5+(6.5*((37/100)^1.5)))*10^6;
Tc= (0.5+(4.5*((37/100)^1.5)))*10^6;
pbar=pc/pf;

P=0.000001:0.000001:0.00004;

Mm =(0.1295.*P).^(1/2);
plot(P,Mm)
hold on

b=(((0.4839/(9*(pi^2)*((SigmaC/SigmaF)^2)*(Ef/SigmaF)))).^(1/4))

Me= 4*(((0.4839/(9*(pi^2)*((SigmaC/SigmaF)^2)*(Ef/SigmaF))))^(1/4)).*(P.^(3/4));
plot(P,Me,'--')
%axis([0 0.00001 0 0.001])
xlabel Pmin
ylabel Mmin
title ('Mmin vs Pmin')
legend('Micro Buckling','Elastic Indentation')
hold off
