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


