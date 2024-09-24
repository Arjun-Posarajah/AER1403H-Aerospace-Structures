% -------------------------------------------------------------------------

% Matlab code for Project 2 - Meshing

% -------------------------------------------------------------------------

% Parameter of the given Semi-Circular Bracket.

radius_inner = 35E-3; 
radius_outer = 55E-3; 
beam_width = 20E-3; 
thickness = 2E-3; 
n_theta = 80;                    % Circumferential Number of Elements
n_radial = 20;                   % Radial Number of Elements

% Material Properties of the bracket.

E = 200e9; 
nu = 0.3; 

% Applied load

applied_load = 4e3; % in Newtons

% Generate coordinates for nodes

theta = linspace(-pi/2, pi/2, n_theta + 1);                 % Angle range from -pi/2 to pi/2
r = linspace(radius_inner, radius_outer, n_radial + 1);     % Radial values from inner to outer radius

% Create nodes 

nodes = zeros((n_theta + 1) * (n_radial + 1), 2);
n = 1;
for i = 1:(n_radial + 1)
    for j = 1:(n_theta + 1)
        nodes(n, 1) = r(i) * cos(theta(j));        % Convert polar to cartesian co-ordinates
        nodes(n, 2) = r(i) * sin(theta(j));
        n = n + 1;
        
    end
end

% The position of each node in the mesh to draw isoparametric elements

elements = zeros(n_theta * n_radial, 5);
element_index = 1;
for i = 1:n_radial
    for j = 1:n_theta
        node1 = (i - 1) * (n_theta+1) + j;      %When dividing by 10 radial elements, 11 will happen which is why we add by 1      
        node2 = node1 + 1;
        node3 = node1 + n_theta + 2;
        node4 = node1 + n_theta + 1;
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