
%problem 3
'region of stability of multistep methods'
clf;
[X,Y] = meshgrid(linspace(-10,10), linspace(-10,10));
Z = X+Y*1i;
phi = ((-3*Z -4)) ./ (Z - 4);
contourf(X,Y,1-abs(phi), [0 0], 'LineWidth', 1);
set(gca, 'FontSize', 20, 'CLim', [0 1]);
colormap([1 .5 .8; 0 0 0]);
hold on;
plot([-3 2], [0 0], '--k', 'LineWidth', 1);
plot([0 0], [-2 2], '--k', 'LineWidth', 1);



%problem 4


%the following code was put as a function
'adam bashforth 2'
function [ x, y ] = ab2 ( f, xrange, y0, numSteps )
x(1) = xrange(1);											%first x value is starting interval
h = ( xrange(2) - xrange(1) ) / numSteps;					%step size
y(:,1) = y0;												%set first y value to given y0
i = 1;														%iterate starting at the first step i = 1
f_Val =  f( x(i), y(:,i) );									%get the f-value 
x_half = x(i) + 0.5 * h;									%get the next x value to be used
y_half = y(:,i) + 0.5 * h * f_Val;							%get the next y-value to be used
f_half = f( x_half, y_half );								%get the f-value of th enext value to be used
x(1,i+1) = x(1,i) + h;										%store this value in hte x-array
y(:,i+1) = y(:,i) + h * f_half;								%store this value in the y-array
									%we now have two initial values wi and wi+1, can do algorithm now
for i = 2 : numSteps				%iterate from the number of steps we will do
  fValueold=f_Val;					%update the f-value
  f_Val = f( x(i), y(:,i) );		%compute 
  x(1,i+1) = x(1,i) + h;
  y(:,i+1) = y(:,i) + h * ( 3 * f_Val - fValueold ) / 2;
end

%-----end code.

%[x, y] = ab2(f, [0 1], 0, 10) is passed to the function, where f is the desired y(t), xrange is the range [0,1],
% y0 is initial value, and numsteps is number of steps; had mesh size as 0.1 here, so numsteps is 10

%for the latter, we run

f = @(t,y) (t+2)^-2 * sin(y);
[x, y] = ab2(f, [0 1], 0, 10);
plot(x,y);




%problem 5
'eulers method'
function [ x, y ] = euler ( f, xrange, y0, numSteps )
x(1) = xrange(1);
h = ( xrange(2) - xrange(1) ) / numSteps;
y(:,1) = y0;
for i = 1 : numSteps
  x(1,i+1) = x(1,i) + h;
  y(:,i+1) = y(:,i) + h * f( x(i), y(:,i) );
end


