nx = 25;
ny = 80;
inner_radius = 35e-3;   % 35 mm converted to meters
bracket_width = 20e-3;  % 20 mm converted to meters

% Calculate outer radius
outer_radius = inner_radius + bracket_width;

% Parametric Equation for a Semicircle (from top to bottom)
theta = linspace(-pi/2, pi/2, ny);
inner_x = inner_radius * cos(theta);
inner_y = inner_radius * sin(theta) + 0.055;
outer_x = outer_radius * cos(theta);
outer_y = outer_radius * sin(theta) + 0.055;

% Generate quadrilateral mesh
X = zeros(ny, nx);
Y = zeros(ny, nx);

% The position of each node in the mesh to draw isoparametric elements

elements = zeros(ny * nx, 5);
x_values = zeros(ny * nx, 5);
y_values = zeros(ny * nx, 5);
element_index = 1;

for i = 1:nx
    for j = 1:ny
        node1 = (i - 1) * (ny + 1) + j;
        node2 = node1 + 1;
        node3 = node1 + ny + 2;
        node4 = node1 + ny + 1;

        % Check that indices are within bounds
        if max([node1, node2, node3, node4]) <= numel(x)
            % Store node indices
            elements(element_index, :) = [node1, node2, node3, node4, node1];

            % Store x and y values for each node
            x_values(element_index, :) = [x(node1), x(node2), x(node3), x(node4), x(node1)];
            y_values(element_index, :) = [y(node1), y(node2), y(node3), y(node4), y(node1)];

            element_index = element_index + 1;
        else
        end
    end
end


for i = 1:ny
    for j = 1:nx
        X(i, j) = inner_x(i) + (outer_x(i) - inner_x(i)) * (j - 1) / (nx - 1);
        Y(i, j) = inner_y(i) + (outer_y(i) - inner_y(i)) * (j - 1) / (nx - 1);
    end
end

% Plotting the Semicircles and Mesh
figure;
hold on;
plot(inner_x, inner_y, 'b', 'LineWidth', 2);
plot(outer_x, outer_y, 'r', 'LineWidth', 2);

% Plot quadrilateral elements
for i = 1:ny-1
    for j = 1:nx-1
        plot([X(i, j), X(i, j+1), X(i+1, j+1), X(i+1, j), X(i, j)],[Y(i, j), Y(i, j+1), Y(i+1, j+1), Y(i+1, j), Y(i, j)], 'k-', 'LineWidth', 0.5);
    end
end

axis equal;
title('Top-to-Bottom Inner and Outer Semicircles with Quadrilateral Mesh');
legend('Inner Semicircle', 'Outer Semicircle', 'Mesh Lines');
xlim([0 0.055]);
ylim([0 0.11]);
hold off;


