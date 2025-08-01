% Input: data points and end derivatives
x = [-1 0 1 2 3];                       % x-values
y = [-28 -10 6 26 56];                  % y-values
L = -21;                                % derivative at x1
R = -37;                                % derivative at xn
natural = 0;                            % flag variable for natural spline
n = length(x) - 1;                      % Number of spline intervals
slope = @(k)((y(k+1)-y(k))/(x(k+1)-x(k)));

% Setup the system
A = zeros(n+1);               % Coefficient matrix (5x5)
rhs = zeros(n+1,1);           % Right-hand side

if ~natural
    A(1,1) = (x(2)-x(1))/3; A(1,2) = (x(2)-x(1))/6;
    A(n+1,n) = (x(n+1)-x(n))/6; A(n+1,n+1) = (x(n+1)-x(n))/3;
    rhs(1) = slope(1) - L; 
    rhs(n+1) = R - (slope(n)); 
end


for i = 2:n
    j = i-1; A(i,j) = (x(j+1)-x(j))/6;
    j = i; A(i,j) = (x(j+1)-x(j-1))/3;
    j = i+1; A(i,j) = (x(j)-x(j-1))/6;
    rhs(i) = (slope(i)) - (slope(i-1));
end


% Solving z values
if natural
    A_n = A(2:n,2:n);
    rhs_n = rhs(2:n);
    z = A_n\rhs_n;
    z = [0;z;0];
else
    z = A\rhs;
end 



% Spline function
syms t
s = [];
for i = 1:n
    c = slope(i) - ( z(i+1)/6 + z(i)/3 )*( x(i+1)-x(i) );
    b = z(i)/2;
    a = (1/6)*( (z(i+1)-z(i)) / (x(i+1)-x(i)) );
    s = vertcat(s,y(i) + c*(t-x(i)) + b*(t-x(i))^2 + a*(t-x(i))^3);
end 

% Plotting
plot(x,y,'rx')
hold on
for i = 1:n
    fplot(s(i),[x(i),x(i+1)],'b-')
end

% Reference
poly = [1 -1 16 -10];
x_span = linspace(x(1),x(n+1));
plot(x_span,polyval(poly,x_span),'r-')