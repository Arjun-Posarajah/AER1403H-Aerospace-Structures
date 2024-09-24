% Define the symbolic variables
syms x y lambda

% Define the equation
eq1 = x*y + 0.067*x - lambda * (4*x^2*y + 4*x^2*y^2);

% Take the partial derivatives
d_eq_x = diff(eq1, x);
d_eq_y = diff(eq1, y);
d_eq_lambda = diff(eq1, lambda);

% Set up the system of equations
eqs = [d_eq_x == 0, d_eq_y == 0, d_eq_lambda == 0];

% Solve the system of equations
sol = solve(eqs, [x, y, lambda]);

% Extract the solutions
x_solution = sol.x;
y_solution = sol.y;
lambda_solution = sol.lambda;

disp(['x: ', char(x_solution)]);
disp(['y: ', char(y_solution)]);
disp(['lambda: ', char(lambda_solution)]);
