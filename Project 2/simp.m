% Project 2: Finite Element Solution for a Curved Bracket
% Name: Arjun Posarajah Student#: 1004881737

% Parameters/Properties
inner_radius = 35E-3;
outer_radius = 55E-3;
width = 20E-3;
thickness = 2E-3;
numtheta = 4;
numradial = 4;
elasticmodulus = 200e9;
poissonratio = 0.3;
appliedforce = 4e3; % in Newtons

% Generate coordinates for nodes
thetavals = linspace(-pi/2, pi/2, numtheta + 1);
radialvals = flip(linspace(inner_radius, outer_radius, numradial + 1));

% Create nodes
nodes = zeros((numtheta + 1) * (numradial + 1), 2);
nodeindex = 1;
for i = 1:(numradial + 1)
    for j = 1:(numtheta + 1)
        nodes(nodeindex, :) = [radialvals(i) * cos(thetavals(j)), radialvals(i) * sin(thetavals(j))];
        nodeindex = nodeindex + 1;
    end
end

% Generate elements
elements = zeros(numtheta * numradial, 5);
element_index = 1;
for i = 1:numradial
    for j = 1:numtheta
        node1 = (i - 1) * (numtheta + 1) + j;
        node2 = node1 + 1;
        node3 = node1 + numtheta + 2;
        node4 = node1 + numtheta + 1;
        elements(element_index, :) = [node1, node2, node3, node4, node1];
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

% Given geometry is in plane stress
D = (elasticmodulus / (1 - poissonratio^2)) * [1, poissonratio, 0; poissonratio, 1, 0; 0, 0, (1 - poissonratio) / 2];

% Gauss quadrature points
xi = [1/sqrt(3), 1/sqrt(3), -1/sqrt(3), -1/sqrt(3)];
n = [1/sqrt(3), -1/sqrt(3), 1/sqrt(3), -1/sqrt(3)];

% Stiffness matrix assembly
k = sparse((numtheta + 1) * (numradial + 1) * 2, (numtheta + 1) * (numradial + 1) * 2);
DOF = reshape(1:(numtheta + 1) * (numradial + 1) * 2, [], 2);

for i = 1:numradial
    for j = 1:numtheta
        node_indices = elements((i - 1) * numtheta + j, :);
        x_vals = nodes(node_indices, 1);
        y_vals = nodes(node_indices, 2);
        
        J = zeros(2, 2);
        Hstar = cell(1, 4);
        
        for q = 1:4
          J = J + [xi(q); n(q)] * [x_vals; y_vals]';
            invJ = inv(J);
            Hstar{q} = invJ * [n(q) - 1, 1 - n(q), 1 + n(q), -n(q) - 1; xi(q) - 1, -xi(q) - 1, 1 + xi(q), 1 - xi(q)];
        end
        
        detJ = det(J);
        H = cellfun(@(Hstar_q) [Hstar_q(1, :), zeros(1, 4); zeros(1, 4), Hstar_q(2, :); Hstar_q(2, :), Hstar_q(1, :)], Hstar, 'UniformOutput', false);
        
        k_star = detJ * cellfun(@(H_q) (H_q' * D * H_q), H, 'UniformOutput', false);
        
        k(DOF(node_indices, :), DOF(node_indices, :)) = k(DOF(node_indices, :), DOF(node_indices, :)) + cell2mat(k_star);
    end
end

% Fixed Points
dof_fixed = reshape(1:(numtheta + 1) * 2:(numtheta + 1) * (numradial + 1) * 2, [], 1);
dof_fixed = [dof_fixed, dof_fixed + 1, dof_fixed + (numtheta + 1) * 2, dof_fixed + (numtheta + 1) * 2 + 1];

dof_all = 1:(numtheta + 1) * (numradial + 1) * 2;
free_dof = setdiff(dof_all, dof_fixed(:));

% Stiffness matrix for the given linear semi-circular beam
K = full(k(free_dof, free_dof));

% Force position finder function:
% Find the middle node of the middle thickness of the middle shape
fnode = (numradial - 1) * (numtheta + 1) + floor(numtheta / 2);
forcedof = [fnode * 2, fnode * 2 + 1];

% Apply the force at the closest node
force_magnitude = -4000; % 4 kN
applied_force_vector = zeros((numtheta + 1) * (numradial + 1) * 2, 1);
applied_force_vector(forcedof) = force_magnitude;

displacement = (1/2) * (K \ applied_force_vector(free_dof));

% Deformed shape
newnodes = nodes;
newnodes(free_dof) = newnodes(free_dof) + displacement;

% Plotting the original and deformed shapes
figure;
hold on;
for i = 1:size(elements, 1)
    n = elements(i, :);
    x = nodes(n, 1);
    y = nodes(n, 2);
    plot(x, y, 'k'); % Plotting edges of each element
end

for i = 1:size(elements, 1)
    n = elements(i, :);
    x = newnodes(n, 1);
    y = newnodes(n, 2);
    plot(x, y, 'r'); % Plotting deformed shape
end

xlabel('X-axis');
ylabel('Y-axis');
title('Deformed Shape of Semicircular Bracket');
axis equal;
hold off;