% Project 2: Finite Element Solution for a Curved Bracket
% Name: Arjun Posarajah Student#:1004881737

%%%% Meshing Section %%%%
% Parameters 
inner_radius = 35E-3; 
outer_radius = 55E-3; 
width = 20E-3; 
thickness = 2E-3; 
num_theta = 40;                    % Circumferential Number of Elements
num_radial = 10;                   % Radial Number of Elements

% Material Properties
elastic_modulus = 200e9; 
poisson_ratio = 0.3; 

% Applied load
applied_force = 4e3; % in Newtons

% Generate coordinates for nodes
theta_vals = linspace(-pi/2, pi/2, num_theta + 1);                 % Angle range from -pi/2 to pi/2
radial_vals = linspace(inner_radius, outer_radius, num_radial + 1);     % Radial values from inner to outer radius

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
        node1 = (i - 1) * (num_theta + 1) + j;      
        node2 = node1 + 1;
        node3 = node1 + num_theta + 2;
        node4 = node1 + num_theta + 1;
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

% Material Properties
elastic_modulus = 200e9;
poisson_ratio = 0.3;

% Initialize global stiffness matrix
num_nodes = size(nodes, 1);
K_global = zeros(2 * num_nodes, 2 * num_nodes);

% Shape functions for a bilinear quadrilateral element
shapeFunctions = @(xi, eta) [0.25 * (1 - xi) * (1 - eta), 0.25 * (1 + xi) * (1 - eta), 0.25 * (1 + xi) * (1 + eta), 0.25 * (1 - xi) * (1 + eta)];
shapeFunctionsDerivativeXi = @(xi, eta) [-0.25 * (1 - eta), 0.25 * (1 - eta), 0.25 * (1 + eta), -0.25 * (1 + eta)];
shapeFunctionsDerivativeEta = @(xi, eta) [-0.25 * (1 - xi), -0.25 * (1 + xi), 0.25 * (1 + xi), 0.25 * (1 - xi)];

% Loop over each element
for elem_idx = 1:size(elements, 1)
    % Extract nodes for the current element
    element_nodes = elements(elem_idx, :);
    
    % Initialize element stiffness matrix
    k_element = zeros(8, 8);
    
    % Loop over integration points (using a simple 2x2 Gauss quadrature)
    for xi = [-1, 1]
        for eta = [-1, 1]
            % Calculate shape functions and derivatives at the integration point
            N = shapeFunctions(xi, eta);
            dN_dxi = shapeFunctionsDerivativeXi(xi, eta);
            dN_deta = shapeFunctionsDerivativeEta(xi, eta);
            
            % Jacobian matrix
            J = [dN_dxi .* nodes(element_nodes, 1), dN_dxi .* nodes(element_nodes, 2);
                 dN_deta .* nodes(element_nodes, 1), dN_deta .* nodes(element_nodes, 2)];
            
            % Derivative of the shape functions with respect to physical coordinates
            dN_dx = J' \ [dN_dxi; dN_deta];
            
            % B matrix (strain-displacement matrix)
            B = zeros(3, 8);
            B(1, 1:2:end) = dN_dx(:, 1)';
            B(2, 2:2:end) = dN_dx(:, 2)';
            B(3, 1:2:end) = dN_dx(:, 2)';
            B(3, 2:2:end) = dN_dx(:, 1)';
            
            % Element stiffness matrix
