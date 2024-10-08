% AER1403 Assignment 8 Q2

% AER1403 Assignment 8 Q2

E = 120*10^9;
nu = 0.3;
coords = [0.1 0.3 1.2 0.9; 1.3 0.4 0.1 1.1];

syms zeta eta;

N1 = 0.25*(1 - zeta)*(1 - eta);
N2 = 0.25*(1 + zeta)*(1 - eta);
N3 = 0.25*(1 + zeta)*(1 + eta);
N4 = 0.25*(1 - zeta)*(1 + eta);

J = 0.25*[eta-1 1-eta 1+eta -1-eta; zeta-1 -1-zeta 1+zeta 1-zeta]*[transpose(coords)];

D = (E/(1 - nu^2))*[1 nu 0; nu 1 0; 0 0 0.5*(1-nu)];

%Ae_1 = abs(0.5*det([0.1 1.3 1;0.3 0.4 1;0.9 1.1 1]));
%Ae_2 = abs(0.5*det([1.2 1.3 1;0.1 0.4 1;0.9 1.1 1]));
%Ae = Ae_1+Ae_2;

H_temp = inv(J)*[diff(N1,zeta) diff(N2,zeta) diff(N3,zeta) diff(N4,zeta); diff(N1,eta) diff(N2,eta) diff(N3,eta) diff(N4,eta)]

H = ([H_temp(1,1) 0 H_temp(1,2) 0 H_temp(1,3) 0 H_temp(1,4) 0; 0 H_temp(2,1) 0 H_temp(2,2) 0 H_temp(2,3) 0 H_temp(2,4); H_temp(2,1) H_temp(1,1) H_temp(2,2) H_temp(1,2) H_temp(2,3) H_temp(1,3) H_temp(2,4) H_temp(1,4)]);

K_temp = transpose(H)*D*H;

n_gp = [-0.7745966692 0 0.7745966692];
wts = [0.5555555556 0.8888888889 0.5555555556];

fin_K = zeros(size(K_temp));
for i =1:length(n_gp)
    for j = 1:length(n_gp)
        zeta = n_gp(i);
        wt_z = wts(i);

        eta = n_gp(j);
        wt_e = wts(j);        

        fin_K = fin_K + wt_z*wt_e*subs(K_temp*det(J));
    end
end

K = double(fin_K)
det(J)