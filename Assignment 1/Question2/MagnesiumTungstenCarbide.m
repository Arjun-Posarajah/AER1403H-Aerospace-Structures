%%Magnesium and Tungsten Carbide

vf=[0:0.01:1]; %vf=Volume Fraction of reinforcing material
pf=15.63;        %pf=density of reinforcing material 
pm=1.8;        %pm=density of the matrix material
Em=42;         %Em=youngs matrix material
Ef=696;        %Ef=youngs reinforcing material 

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
title('Modulus Bounds: Magnesium, Aluminum and Tungsten Carbide');
xlabel('Density ,p (Mg/m^3)')
ylabel('Youngs Modulus, E (GPa)')
hold off