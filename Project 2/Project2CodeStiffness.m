%Project 2: Finite Element Solution for a Curved Bracket
%Name: Arjun Posarajah Student#:1004881737

% Parameters/Properties
inner_radius = 35E-3; 
outer_radius = 55E-3; 
width = 20E-3; 
thickness = 2E-3; 
num_theta = 20;                    % Circumferential Number of Elements
num_radial = 5;                   % Radial Number of Elements
elastic_modulus = 200e9; 
poisson_ratio = 0.3; 
applied_force = 4e3; % in Newtons

% Generate coordinates for nodes
theta_vals = linspace(-pi/2, pi/2, num_theta + 1);                 % Angle range from -pi/2 to pi/2
radial_vals = flip(linspace(inner_radius, outer_radius, num_radial + 1));     % Radial values from inner to outer radius

% Create nodes 
nodes = zeros((num_theta + 1) * (num_radial + 1), 2);
node_index = 1;
for i = 1:(num_radial + 1)
    for j = 1:(num_theta + 1)
        nodes(node_index, 1) = radial_vals(i) * cos(theta_vals(j));        % Convert polar to Cartesian coordinates
        nodes(node_index, 2) = radial_vals(i) * sin(theta_vals(j));
        node_index = node_index + 1;
    end
end

% The position of each node in the mesh to draw isoparametric elements
elements = zeros(num_theta * num_radial, 5);
element_index = 1;
for i = 1:num_radial
    for j = 1:num_theta
        node1 = (i - 1) * (num_theta+1) + j;      %Example: When dividing by 10 radial elements, 11 nodes will happen which is why we add by 1      
        node2 = node1 + 1;
        node3 = node1 + num_theta + 2;
        node4 = node1 + num_theta + 1;
        elements(element_index, :) = [node1, node2, node3, node4, node1];       % Store contour numbers element by element column wise.
        element_index = element_index + 1;
    end
end

% Plotting the mesh
figure;
hold on;
for i = 1:size(elements, 1)
    n = elements(i, :);
    x = nodes(n, 1);
    y = nodes(n, 2);
    plot(x, y, 'k'); % Plotting edges of each element
end
xlabel('X-axis');
ylabel('Y-axis');
title('Mesh of Semicircular Bracket');
axis equal;
hold off;

% Given geometry is in plane stress:

D = (elastic_modulus/1-((poisson_ratio)^2)) * [1, poisson_ratio, 0;poisson_ratio, 1, 0;0, 0, ((1-poisson_ratio)/2)];       % Pa

%Calculating global stiffness matrix
%Using 2-point Gauss quadrature to calculate each GN4Q value:

xi = [(1/sqrt(3)) (1/sqrt(3)) -(1/sqrt(3)) -(1/sqrt(3))];
n = [(1/sqrt(3)) -(1/sqrt(3)) (1/sqrt(3)) -(1/sqrt(3))];

GN4Q1 = (1/4) * [n(1)-1 1-n(1) 1+n(1) -n(1)-1; xi(1)-1 -xi(1)-1 1+xi(1) 1-xi(1)]; 
GN4Q2 = (1/4) * [n(2)-1 1-n(2) 1+n(2) -n(2)-1; xi(2)-1 -xi(2)-1 1+xi(2) 1-xi(2)];
GN4Q3 = (1/4) * [n(3)-1 1-n(3) 1+n(3) -n(3)-1; xi(3)-1 -xi(3)-1 1+xi(3) 1-xi(3)]; 
GN4Q4 = (1/4) * [n(4)-1 1-n(4) 1+n(4) -n(4)-1; xi(4)-1 -xi(4)-1 1+xi(4) 1-xi(4)]; 

%Allocating memory for the global stiffness matrix
TotalDOFs= (num_theta+1)*(num_radial+1)*2;
Totalelements=num_theta*num_radial;
k = sparse(TotalDOFs, TotalDOFs);
DOF=zeros(Totalelements,8);

XDOFs = (2 .* (elements(:,1:4))) - 1;
YDOFs = (2 .* (elements(:,1:4)));

for i = 1:Totalelements
        DOF(i,:) = [XDOFs(i,1); YDOFs(i,1); XDOFs(i,2); YDOFs(i,2); XDOFs(i,3); YDOFs(i,3); XDOFs(i,4); YDOFs(i,4);]; 
end

%Single Loop iteration calculates individual element stiffness matrix and
% reconciles that to the global stiffness matrix
xcoord= nodes(:,1);
ycoord= nodes(:,2);

for i = 1: Totalelements
    coord_ind = elements(i,:);
    Xval = xcoord(coord_ind);
    Yval = ycoord(coord_ind);
    coordmat = [Xval(1) Yval(1); Xval(2) Yval(2); Xval(3) Yval(3); Xval(4) Yval(4)];
    J1 = GN4Q1 * coordmat;
    J2 = GN4Q2 * coordmat;
    J3 = GN4Q3 * coordmat;
    J4 = GN4Q4 * coordmat;
    detJ1 = det(J1);
    invJ1 = inv(J1);
    detJ2 = det(J2);
    invJ2 = inv(J2);
    detJ3 = det(J3);
    invJ3 = inv(J3);
    detJ4 = det(J4);
    invJ4 = inv(J4);
    Hstar1 = invJ1\GN4Q1;
    Hstar2 = invJ2\GN4Q2;
    Hstar3 = invJ3\GN4Q3;
    Hstar4 = invJ4\GN4Q4;
    H1 = [Hstar1(1,1) 0 Hstar1(1,2) 0 Hstar1(1,3) 0 Hstar1(1,4) 0; 0 Hstar1(2,1) 0 Hstar1(2,2) 0 Hstar1(2,3) 0 Hstar1(2,4); Hstar1(2,1) Hstar1(1,1) Hstar1(2,2) Hstar1(1,2) Hstar1(2,3) Hstar1(1,3) Hstar1(2,4) Hstar1(1,4)];
    H2 = [Hstar2(1,1) 0 Hstar2(1,2) 0 Hstar2(1,3) 0 Hstar2(1,4) 0; 0 Hstar2(2,1) 0 Hstar2(2,2) 0 Hstar2(2,3) 0 Hstar2(2,4); Hstar2(2,1) Hstar2(1,1) Hstar2(2,2) Hstar2(1,2) Hstar2(2,3) Hstar2(1,3) Hstar2(2,4) Hstar2(1,4)];
    H3 = [Hstar3(1,1) 0 Hstar3(1,2) 0 Hstar3(1,3) 0 Hstar3(1,4) 0; 0 Hstar3(2,1) 0 Hstar3(2,2) 0 Hstar3(2,3) 0 Hstar3(2,4); Hstar3(2,1) Hstar3(1,1) Hstar3(2,2) Hstar3(1,2) Hstar3(2,3) Hstar3(1,3) Hstar3(2,4) Hstar3(1,4)];
    H4 = [Hstar4(1,1) 0 Hstar4(1,2) 0 Hstar4(1,3) 0 Hstar4(1,4) 0; 0 Hstar4(2,1) 0 Hstar4(2,2) 0 Hstar4(2,3) 0 Hstar4(2,4); Hstar4(2,1) Hstar4(1,1) Hstar4(2,2) Hstar4(1,2) Hstar4(2,3) Hstar4(1,3) Hstar4(2,4) Hstar4(1,4)];
    kstar = ((H1.')*D*H1*detJ1) + ((H2.')*D*H2*detJ2) + ((H3.')*D*H3*detJ3) + ((H4.')*D*H4*detJ4);    % Pa
    elementDOF = DOF(i,:);
    k(elementDOF,elementDOF) = k(elementDOF,elementDOF) + kstar;      % Pa
end

