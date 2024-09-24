t=0.125;%mm
d=6*10.^-6;%um
vf=0.65;
E1=131;
E2=10.3;
G12=6;
v12=0.22;
theta=30;

%Find the stiffneses matrix from the lamina in the material coordinate system
%Determine the angle theta by which the global coordinate system is rotatedthe material system
%Use the equation on page 41 to determine Qij components of stiffness matrix

S11=1/E1;
S22=1/E2;
S12=(-v12)/E1;
S66=1/G12;


Q11=S22/((S11*S22)-(S12.^2));
Q22=S11/((S11*S22)-(S12.^2));
Q12=(-S12)/((S11*S22)-(S12.^2));
Q66=G12;

Q=[Q11 Q12 0; Q12 Q22 0; 0 0 Q66];

disp(Q)