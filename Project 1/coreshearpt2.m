% Define the symbolic variables
pf= 1600+(30*48);
Ef= (40+81)*10^9;
SigmaF= (200+(100*7))*10^6;
pc= 20+(5*37);
SigmaC= (0.5+(6.5*((37/100)^1.5)))*10^6;
Tc= (0.5+(4.5*((37/100)^1.5)))*10^6;
syms x y lambda

% Define the equation
eq1 = x*y + (pc/pf)*x - lambda*(((2*Tc)/SigmaF)*(x+(x*y)));

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
