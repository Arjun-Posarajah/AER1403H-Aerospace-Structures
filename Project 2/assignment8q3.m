clc, clearvars, close all
syms e z R1x R1y R2x R2y d3x d3y d4x d4y
t = 0.01;
E = 70*10^9; %in GPa
v = 0.3;
%Matrix of Elastic Constants (Plane Strain)
D = E/(1+v)/(1-2*v).*[1-v,v,0; v,1-v,0; 0,0,(1-2*v)/2];
coord = [0.2,0; 0.8,0; 0.6,0.8; 0,0.6];
N1 = 0.25*(1 - z)*(1 - e);
N2 = 0.25*(1 + z)*(1 - e);
N3 = 0.25*(1 + z)*(1 + e);
N4 = 0.25*(1 - z)*(1 + e);
N1_z = diff(N1,z);
N1_e = diff(N1,e);
N2_z= diff(N2,z);
N2_e = diff(N2,e);
N3_z = diff(N3,z);
N3_e = diff(N3,e);
N4_z = diff(N4,z);
N4_e = diff(N4,e);
N =[N1_z,N2_z,N3_z,N4_z; N1_e,N2_e,N3_e,N4_e];
J = N*coord;
Hp = J\N;
H = [Hp(1,1),0,Hp(1,2),0,Hp(1,3),0,Hp(1,4),0; 0,Hp(2,1),0,Hp(2,2),0,Hp(2,3),0,Hp(2,4); Hp(2,1),Hp(1,1),Hp(2,2),Hp(1,2),Hp(2,3),Hp(1,3),Hp(2,4),Hp(1,4)];
n_gp = [1/sqrt(3),-1/sqrt(3)];
K = zeros(size(8,8));
for i=1:2
for j=1:2
z = n_gp(i);
e = n_gp(j);
K = K + subs(transpose(H)*D*det(J)*H);
end
end
Kt = t.*double(K);
disp(Kt)