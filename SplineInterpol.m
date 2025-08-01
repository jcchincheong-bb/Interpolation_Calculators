% Clamped Cubic Spline with 5 data points

% Input: data points and end derivatives
x = [1 2 3 4 5];              % example x-values
y = [2 3 5 4 2];              % example y-values
fp0 = 1.0;                    % derivative at x1
fpn = -1.0;                   % derivative at x5

n = length(x) - 1;            % Number of spline intervals (4)
h = diff(x);                  % Step sizes

% Step 1: Setup the system
A = zeros(n+1);               % Coefficient matrix (5x5)
rhs = zeros(n+1,1);           % Right-hand side

% Clamped conditions
A(1,1) = 2*h(1);
A(1,2) = h(1);
rhs(1) = 3*( (y(2) - y(1))/h(1) - fp0 );

A(n+1,n) = h(n);
A(n+1,n+1) = 2*h(n);
rhs(n+1) = 3*( fpn - (y(end) - y(end-1))/h(n) );

% Interior equations
for i = 2:n
    A(i,i-1) = h(i-1);
    A(i,i)   = 2*(h(i-1) + h(i));
    A(i,i+1) = h(i);
    rhs(i) = 3*((y(i+1) - y(i))/h(i) - (y(i) - y(i-1))/h(i-1));
end

% Step 2: Solve for c coefficients
c = A \ rhs;

% Step 3: Compute a, b, d for each spline interval
a = y(1:n);
b = zeros(n,1);
d = zeros(n,1);

for i = 1:n
    b(i) = (y(i+1) - y(i))/h(i) - h(i)*(2*c(i) + c(i+1))/3;
    d(i) = (c(i+1) - c(i)) / (3*h(i));
end

% Output spline coefficients for each interval
% S_i(x) = a(i) + b(i)*(x - x(i)) + c(i)*(x - x(i))^2 + d(i)*(x - x(i))^3

fprintf('Interval spline coefficients:\n');
for i = 1:n
    s{i} = [a(i) b(i) c(i) d(i)];
    fprintf('Interval [%f, %f]:\n', x(i), x(i+1));
    fprintf('  a = %.4f\n', a(i));
    fprintf('  b = %.4f\n', b(i));
    fprintf('  c = %.4f\n', c(i));
    fprintf('  d = %.4f\n\n', d(i));
    
end

for i = 1:n
    x_span = linspace(x(i),x(i+1));
    plot(x_span, polyval(s{i},x_span));
    hold on
end
plot(x,y,'rx')