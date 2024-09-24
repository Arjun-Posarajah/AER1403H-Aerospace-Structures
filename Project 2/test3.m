%Project 2: Finite Element Solution for a Curved Bracket
%Name: Arjun Posarajah Student#:1004881737

%%%% Meshing Section %%%%
% Parameters/Properties
inner_radius = 35E-3; 
outer_radius = 55E-3; 
width = 20E-3; 
thickness = 2E-3; 
num_theta = 40;                    % Circumferential Number of Elements
num_radial = 10;                   % Radial Number of Elements
elastic_modulus = 200e9; 
poisson_ratio = 0.3; 
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

% Identify fixed nodes
fixed_nodes = [1, num_theta + 1, (num_theta + 1) * num_radial + 1, (num_theta + 1) * (num_radial + 1)];

% Initialize global stiffness matrix
num_nodes = (num_theta + 1) * (num_radial + 1);
K_global = zeros(2 * num_nodes, 2 * num_nodes);

% Loop over elements
for elem = 1:size(elements, 1)
    % Get nodal coordinates of the element
    nodes_elem = elements(elem, :);
    x_elem = nodes(nodes_elem, 1);
    y_elem = nodes(nodes_elem, 2);
    
    % Compute element stiffness matrix for quadrilateral element
    [K_elem] = compute_quadrilateral_element_stiffness(x_elem, y_elem, D);
    
    % Assemble element stiffness matrix into the global stiffness matrix
    K_global = assemble_global_stiffness(K_global, K_elem, nodes_elem);
end

% Apply fixed boundary conditions
K_global_reduced = remove_rows_and_columns(K_global, fixed_nodes);

% Function to remove rows and columns associated with fixed nodes
function [K_reduced] = remove_rows_and_columns(K, fixed_nodes)
    % Remove rows and columns associated with fixed nodes
    K_reduced = K;
    K_reduced(fixed_nodes, :) = [];
    K_reduced(:, fixed_nodes) = [];
end
% Function to compute element stiffness matrix for quadrilateral element
function [K_elem] = compute_quadrilateral_element_stiffness(x, y, D)
 % Assuming linear quadrilateral elements

 % Compute derivatives of shape functions with respect to xi and eta
dN_dxi = [-0.25*(1-y); 0.25*(1-y); 0.25*(1+y); -0.25*(1+y)];
dN_deta = [-0.25*(1-x), -0.25*(1+x), 0.25*(1+x), 0.25*(1-x)];

% Jacobian matrix
J = [dN_dxi * x, dN_deta * x; dN_dxi * y, dN_deta * y];
    
    % Inverse of Jacobian matrix
    invJ = inv(J);
    
    % Derivatives of shape functions with respect to x and y
    dN_dx = invJ(1, 1) * dN_dxi + invJ(1, 2) * dN_deta;
    dN_dy = invJ(2, 1) * dN_dxi + invJ(2, 2) * dN_deta;
    
    % B matrix (strain-displacement matrix)
    B = zeros(3, 8);
    B(1, 1:2:end) = dN_dx;
    B(2, 2:2:end) = dN_dy;
    B(3, 1:2:end) = dN_dy;
    B(3, 2:2:end) = dN_dx;
    
    % Compute element area
    A = polyarea(x, y);
    
    % Compute element stiffness matrix
    K_elem = A * B' * D * B;
end
