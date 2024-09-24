CordMat=[0.2 0.8 0.6 0; 0 0 0.8 0.6];
E=7*10^10;
v=0.3;

syms zeta eta

%Plain Strain
m=E/((1+v)*(1-(2*v)));


%Shape Functions
N1= (0.25)*(1-zeta)*(1-eta);
N2=(0.25)*(1+zeta)*(1-eta);
N3= (0.25)*(1+zeta)*(1+eta);
N4= (0.25)*(1-zeta)*(1+eta);

J = 0.25*[eta-1 1-eta 1+eta -1-eta; zeta-1 -1-zeta 1+zeta 1-zeta]*[transpose(CordMat)];

D= m.*[1-v v 0 ; v 1-v 0; 0 0 (1-(2*v))/2];
te=0.01;

Hstar= inv(J)*[diff(N1,zeta) diff(N2,zeta) diff(N3,zeta) diff(N4,zeta); diff(N1,eta) diff(N2,eta) diff(N3,eta) diff(N4,eta)];

H= [Hstar(1,1) 0 Hstar(1,2) 0 Hstar(1,3) 0 Hstar(1,4) 0; 0 Hstar(2,1) 0 Hstar(2,2) 0 Hstar(2,3) 0 Hstar(2,4); Hstar(2,1) Hstar(1,1) Hstar(2,2) Hstar(1,2) Hstar(2,3) Hstar(1,3) Hstar(2,4) Hstar(1,4)];

Kstar= transpose(H)*D*H;

gp = [-0.577 0.577];
weight = [1 1];
%n_gp = [-0.7745966692 0 0.7745966692];
%wts = [0.5555555556 0.8888888889 0.5555555556];

fin_K = zeros(size(Kstar));
for i =1:length(gp)
    for j = 1:length(gp)
        zeta = gp(i);
        wt_z = weight(i);

        eta = gp(j);
        wt_e = weight(j);        

        fin_K = fin_K + wt_z*wt_e*subs(Kstar*det(J));
    end
end

K = te*double(fin_K)

F= transpose([0 0 0 0 250*10^3 0 0 0]);

kboundcond=K(5:8,5:8);
F=F(5:8);

displacement=kboundcond\F

ReactionForce=K*[0 0 0 0 displacement']'



