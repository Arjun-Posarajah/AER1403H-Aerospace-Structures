% Define three equations
cs= ((2*Tc)/SigmaF).*(cbar)+((2*Tc)/SigmaF).*cbar.*tbar;
mb= 4*(cbar.^2).*tbar+4.*(tbar.^2).*(cbar.^2);
ei=(1/SigmaF).*(tbar.*cbar).*(((((pi^2)*Ef*(SigmaC^2))/3).*cbar.*(1+tbar)).^(1/3));

% Initialize arrays to store the results
cs1 = zeros(100, 100);
mb2 = zeros(100, 100);
ei3 = zeros(100, 100);

% Loop through each point in the meshgrid and evaluate the equations
for i = 1:100
    for j = 1:100
        cs1(i, j) = cs(X(i, j), Y(i, j));
        mb2(i, j) = mb(X(i, j), Y(i, j));
        ei3(i, j) = ei(X(i, j), Y(i, j));
    end
end

% Initialize an array to store the lowest equation index (1, 2, or 3) at each point
min_results = zeros(100, 100);

% Loop through each point and determine the lowest equation
for i = 1:100
    for j = 1:100
        values = [cs1(i, j), mb2(i, j), ei3(i, j)];
        [~, min_idx] = min(values);
        min_results(i, j) = min_idx;
    end
end

% Now you have min_results containing the lowest equation indices (1, 2, or 3) at each point.

% Create the contour plot for the min_results array
contour(cbar, tbar, min_results, [1, 2, 3]);
xlabel('c/l');
ylabel('t/c');
legend('cs', 'mb', 'ei');

% Add labels for the contour lines
clabel(C, 'FontSize', 12);

% Specify line colors for each equation
hfail = gca;
hfail = get(hfail, 'failure');
set(hfailure(1), 'Color', 'red');
set(hfailure(2), 'Color', 'blue');
set(hfailure(3), 'Color', 'green');

% Customize the legend
legend('cs', 'mb', 'ei', 'Location', 'NorthEast');

% Display the colorbar
%colorbar('Ticks', [1, 2, 3], 'TickLabels', {'Equation 1', 'Equation 2', 'Equation 3'});
hold off
