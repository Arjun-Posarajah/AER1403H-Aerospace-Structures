clc, clearvars, close all
% AER1403 Assignment 8 Q1
E = 70e9; %Pa
v = 0.3;


x = [0 1 1 0];
y = [0 0 1 1];

coordinates = [0 0; 1 0; 1 1; 0 1];

x_disp = [0.02 0.05 0.03 -0.02];
y_disp = [0.01 -0.01 0.04 -0.01];

d = [0.02 0.01 0.05 -0.01 0.03 0.04 -0.02 -0.01];
D=(E/(1-(v^2))).*[1 v 0; v 1-v 0; 0 0 ((1-(2*v))/2)];


zeta_vals = -1:0.1:1;
eta_vals = -1:0.1:1;

x_vals = 0:0.05:1;
y_vals = 0:0.05:1;

disp_mat_x = zeros(length(zeta_vals),length(eta_vals));
disp_mat_y = zeros(length(zeta_vals),length(eta_vals));

for i = 1:length(zeta_vals)
    for j = 1:length(eta_vals)

        zeta = zeta_vals(i);
        eta = eta_vals(j);
        
        N1 = 0.25*(1 - zeta)*(1 - eta);
        N2 = 0.25*(1 + zeta)*(1 - eta);
        N3 = 0.25*(1 + zeta)*(1 + eta);
        N4 = 0.25*(1 - zeta)*(1 + eta);

        inter_disp = [N1 N2 N3 N4]*[x_disp' y_disp'];
        disp_mat_x(j,i) = inter_disp(1);
        disp_mat_y(j,i) = inter_disp(2);

        J = (1/4)*[eta-1 1-eta 1+eta -1-eta; zeta-1 -1-zeta 1+zeta 1-zeta]*(coordinates);
        H_star = (1/4)*inv(J)*[eta-1 1-eta 1+eta -1-eta; zeta-1 -1-zeta 1+zeta 1-zeta];
        H = [H_star(1,1) 0 H_star(1,2) 0 H_star(1,3) 0 H_star(1,4) 0; 0 H_star(2,1) 0 H_star(2,2) 0 H_star(2,3) 0 H_star(2,4); H_star(2,1) H_star(1,1) H_star(2,2) H_star(1,2) H_star(2,3) H_star(1,3) H_star(2,4) H_star(1,4)];

        epsilon = H*transpose(d);
        sigma = E.*epsilon;

        sigma_collected(j,i,:) = sigma;
        sigma11 = sigma_collected(:,:,1);
        sigma22 = sigma_collected(:,:,2);
        sigma12 = sigma_collected(:,:,3);
    end
end
sigma_collected



[X,Y] = meshgrid(zeta_vals,eta_vals);

contourf(X,Y,sigma11)
hold on
contourf(X,Y,sigma22)
contourf(X,Y,sigma12)
hold off
xlabel('X-axis')
ylabel('Y-axis')
colorbar

H