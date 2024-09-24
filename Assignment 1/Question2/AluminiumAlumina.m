%%Aluminium and Alumina

vf=[0:0.01:1]; %vf=Volume Fraction of reinforcing material
pf=3.9;        %pf=density of reinforcing material 
pm=2.7;        %pm=density of the matrix material
Em=70;         %Em=youngs matrix material
Ef=340;        %Ef=youngs reinforcing material 

%Reuss Bound
Er= (((vf/Ef)+((1-vf)/Em)).^-1);

%Voigt Bound
Ev=(vf*Ef)+((1-vf)*Em);

%Density
p=(vf*pf)+((1-vf)*pm);

%Graph
hold on
plot(p,Er)
plot(p,Ev)
title('Modulus Bounds: Aluminum and Alumina');
xlabel('Density ,p (Mg/m^3)')
ylabel('Youngs Modulus, E (GPa)')
hold off

