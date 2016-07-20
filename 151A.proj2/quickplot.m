% this script does the following things:

%(1)
%Calls proj2.m using both options for boundary conditions and stores the resulting
%piecewise cubic spline polynomials into a local function handle, and as well as their
%piecewise components into a cell/array.

%(2)
% creates an array/set of t-values for which we want to evaluate the functions at; we
% then evaluate the functions at those t-values. 

%(3) plots x vs t
%(4) plots y vs t
%(5) plots x vs y, with direction!

%(6) uses 3-pt midpoint and endpoint formulas to numerically compute derivative and plots 
% that against the derivative kind of evaluated through the piecewise cubic spline polynomial,
% which essentially plots two versions of the vector field.

%hence, the main purpose of this script is to plot things.


%Calls proj2.m using both options for boundary conditions and stores the resulting
%piecewise cubic spline polynomials into a local function handle, and as well as their
%piecewise components into a cell/array.
[spxi, spyi, sxi, syi] = proj2(1, 'data.mat');
[spxii, spyii, sxii, syii] = proj2(2, 'data.mat');

T = [];
xi = [];
yi = [];
xii = [];
yii = [];

for i = 0:0.05:6.2
	T(end + 1) = i;
	xi(end + 1) = spxi(i);
	yi(end + 1) = spyi(i);
	xii(end + 1) = spxii(i);
	yii(end + 1) = spyii(i);
end

%give direction to the x vs y graph we will plot later
dxi = gradient(xi);
dxii = gradient(xii);
dyi = gradient(yi);
dyii = gradient(yii);

%plot x and t, then y and t, then x vs y
figure(1);
plot(T, xi, 'color', [0 0 1]);
hold on;
plot(T, xii, 'color', [1 0 0]);

figure(2);
plot(T, yi, 'color', [0 0 1]);
hold on;
plot(T, yii, 'color', [1 0 0]);

figure(3);
quiver(xi, yi, dxi, dyi, 'color', [0 0 1]);
hold on;
quiver(xii, yii, dxii, dyii, 'color', [1 0 0]);


%experiments part with velocity
data = load('data.mat');
Data = data.ip;
[n, cols] = size(Data);

t = [];
x = [];
y = [];
for i = 1:n
	t(end + 1) = Data(i, 1);
	x(end + 1) = Data(i, 2);
	y(end + 1) = Data(i, 3);
end

h = t(2) - t(1);

%use midpoint formula to generate velocity vectors v
vx = [];
vy = [];
vx(end + 1) = (1/(2*h))*(-3*x(1) + 4*x(2) - x(3)); %handle first point
vy(end + 1) = (1/(2*h))*(-3*y(1) + 4*y(2) - y(3));
for i = 2:(n-1)
	vx(end + 1) = (1/(2*h))*(x(i+1) - x(i-1)); %handle middle points 
	vy(end + 1) = (1/(2*h))*(y(i+1) - y(i-1));
end
vx(end + 1) = (1/(2*h))*(-3*x(n) + 4*x(n-1) - x(n-2)); %handle last point
vy(end + 1) = (1/(2*h))*(-3*y(n) + 4*y(n-1) - y(n-2));

%use the cubic splne function to generate velocity vectors u
ux = [];
uy = [];
xplot = [];
yplot = [];
for i = 1:n
	xplot(end + 1) = spxi(t(i));
	yplot(end + 1) = spyi(t(i));
end
ux = gradient(xplot);
uy = gradient(yplot);

%plot x' and y' vector field
figure(4);
quiver(xplot, yplot, ux, uy);
hold on;
quiver(x, y, vx, vy, 'color', [1 0 0]);

