x=[0 1 1 0];
y=[0 0 1 1];
coordinates=[0 0; 1 0; 1 1; 0 1];
xdisp=[0.02 0.05 0.03 -0.02];
ydisp=[0.01 -0.01 0.04 -0.01];

d=[0.02 0.01 0.05 -0.01 0.03 0.04 -0.02 -0.01];

E=7*10^10;
v=0.3;

D=(E/(1-(v^2))).*[1 v 0; v 1-v 0; 0 0 ((1-(2*v))/2)];

zetaVal=-1:0.1:1;
etaVal=-1:0.1:1;

xVal= 0:0.05:1;
yVal=0:0.05:1;

dispMatX= zeros(length(zetaVal),length(etaVal));
dispMatY= zeros(length(zetaVal),length(etaVal));

for i=1:length(zetaVal)
    for j=1:length(etaVal)
        zeta=zetaVal(i);
        eta=etaVal(j);

        N1 = 0.25*(1 - zeta)*(1 - eta);
        N2 = 0.25*(1 + zeta)*(1 - eta);
        N3 = 0.25*(1 + zeta)*(1 + eta);
        N4 = 0.25*(1 - zeta)*(1 + eta);

        interDisp=[N1 N2 N3 N4]*[xdisp' ydisp']
        interDispMatrix(:, j, i) = interDisp'
        dispMatX(j,i)=interDisp(1);
        dispMatY(j,i)=interDisp(2);
        
        J = (1/4)*[eta-1 1-eta 1+eta -1-eta; zeta-1 -1-zeta 1+zeta 1-zeta]*(coordinates);
        Hstar = (1/4)*inv(J)*[eta-1 1-eta 1+eta -1-eta; zeta-1 -1-zeta 1+zeta 1-zeta];
        H = [Hstar(1,1) 0 Hstar(1,2) 0 Hstar(1,3) 0 Hstar(1,4) 0; 0 Hstar(2,1) 0 Hstar(2,2) 0 Hstar(2,3) 0 Hstar(2,4); Hstar(2,1) Hstar(1,1) Hstar(2,2) Hstar(1,2) Hstar(2,3) Hstar(1,3) Hstar(2,4) Hstar(1,4)];

        strain = H*transpose(d);
        stress = E.*strain;

        sigmaColl(j,i,:) = stress;
        sigma11 = sigmaColl(:,:,1);
        sigma22 = sigmaColl(:,:,2);
        sigma12 = sigmaColl(:,:,3);
    end
end
[X,Y]= meshgrid(zetaVal,etaVal)

figure (1)
contourf(X,Y,sigma11)
hold on
contourf(X,Y,sigma22)
contourf(X,Y,sigma12)
hold off
xlabel('X-axis')
ylabel('Y-axis')
colorbar

figure (2)
[X,Y] = meshgrid(x_vals,y_vals);

contourf(X,Y,disp_mat_x,10,'-.','ShowText','on')
xlabel('X-Axis')
ylabel('Y-Axis')
colorbar
H