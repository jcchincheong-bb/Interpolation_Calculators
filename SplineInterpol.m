%% Input: data points and end derivatives
x = [-1 0 1 2 3 4 5 6 7];                       % x-values
y = [-28 -10 6 26 56 38 16 12 11];                  % y-values
L = -21;                                % derivative at x1
R = -37;                                % derivative at xn
natural = 0;                            % flag variable for natural spline
n = length(x) - 1;                      % Number of spline intervals
slope = @(k,x,y)((y(k+1)-y(k))/(x(k+1)-x(k)));

%% Setup the system
A = zeros(n+1);               % Coefficient matrix (5x5)
rhs = zeros(n+1,1);           % Right-hand side

% Clapped Spline set up
if ~natural
    A(1,1) = (x(2)-x(1))/3; A(1,2) = (x(2)-x(1))/6;
    A(n+1,n) = (x(n+1)-x(n))/6; A(n+1,n+1) = (x(n+1)-x(n))/3;
    rhs(1) = slope(1,x,y) - L; 
    rhs(n+1) = R - (slope(n,x,y)); 
end

% Natural Spline Set up
for i = 2:n
    j = i-1; A(i,j) = (x(j+1)-x(j))/6;
    j = i; A(i,j) = (x(j+1)-x(j-1))/3;
    j = i+1; A(i,j) = (x(j)-x(j-1))/6;
    rhs(i) = (slope(i,x,y)) - (slope(i-1,x,y));
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



%% Spline function
s = zeros(n,4);  % Array of polynomials (3rd order)
for i = 1:n
    a = (1/6)*( slope(i,x,z) );
    b = z(i)/2;
    c = slope(i,x,y) - ( z(i+1)/6 + z(i)/3 )*( x(i+1)-x(i) );
    d = y(i);
    ax = a*conv(conv([1 -x(i)],[1 -x(i)]),[1 -x(i)]);
    bx = b*conv([1 -x(i)],[1 -x(i)]);
    cx = c*[1 -x(i)];
    s(i,:) = ax + [0 bx] + [0 0 cx] + [0 0 0 d];
end 

%% Plotting
plot(x,y,'rx')
hold on
for i = 1:n
    x_span = linspace(x(i),x(i+1));
    plot(x_span,polyval(s(i,:),x_span),'b-')
end

%% Interpolated Polynomial - Newton
x_span = linspace(x(1),x(n+1));

% divided difference scheme
diffs = zeros(n+1,n+1);
diffs(:,1) = y(:);                             % First value is the y values

for j = 2:n+1
    for i = 1:n+2-j
        diffs(i,j) = (diffs(i+1,j-1) - diffs(i,j-1)) / (x(i+j-1) - x(i));
    end
end

coeffs = diffs(1,:);

% evaluating the newton polynomial
newpoly_val = zeros(1,length(x_span));
newpoly_val = newpoly_val + coeffs(1);                  % First value is zeroth order
for k = 1:n
    temp_poly = rec_conv(x(1:k));
    newpoly_val = newpoly_val + coeffs(k+1).*polyval(temp_poly,x_span);  % Super position
end

plot(x_span,newpoly_val,'r-')

function out = rec_conv(x)
    if length(x) == 1
        out = [1,-x];
    elseif length(x) == 2
        out = conv([1,-x(1)],[1,-x(2)]);
    else
        out = conv([1,-x(1)],rec_conv(x(2:end)));
    end
end