%theodore nguyen 704-156-701 math 151B Spr2016

%question 1C
'adam bashforth-2'
%the following code was put as a function --------------------------------------

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

%-----end code---------------------------------------------------------------.

% [X,Y] = ab2(f, [0, 2.5], alpha, 250) is passed to plot.

f = @(t, y) sin y + 2y;
alpha0 = 1;
alpha1 = -0.7;
[X, Y] = ab2(f, [0, 2.5], alpha0, 250);
figure 1;
plot(X,Y);
f0 = Y(251);

[X,Y] = ab2(f, [0, 2.5], alpha1, 250);
figure 2;
plot(X,Y);
f1 = Y(251);

alpha2 = alpha1 - (alpha1 - alpha0)*f1/(f0 - f1);

[X,Y] = ab2(f, [0, 2.5], alpha2, 250);
figure 3;
plot(X,Y);
f2 = Y(251);

alpha3 = alpha2 - (alpha2 - alpha1)*f2/(f1 - f2);

[X,Y] = ab2(f, [0, 2.5], alpha3, 250);
figure 4;
plot(X,Y);
f3 = Y(251);

alpha4 = alpha3 - (alpha3 - alpha2)*f3/(f2 - f3);

[X,Y] = ab2(f, [0, 2.5], alpha4, 250);
figure 5;
plot(X,Y);

clear all



%question 2D
'finite difference method'
%y''=p(t)y'+q(t)y+r(t), a<=t<=b, y(a)=alp, y(b)=beta

A = pi/2;
B = (3*pi)/2;
%we want h = 0.01. This gives us N+1=314 or N+1=315
%N+1=314 gives h > 0.01, N+1=315 gives h < 0.01.
%so, we make N+1=315. Therefore, N = 314
N = 314;
alpha = 1;
beta = 1;
p = @(t) sin(t);
q = @(t) -1*cos(t);
r = @(t) -1*sin(t);

h = (B - A)/(N+1);
x = A + h;
a = []; 
b = [];
c = [];
d = [];
a(end + 1) = 2 + q(x)*h^2;
b(end + 1) = -1 + (h/2)*(p(x));
c(end + 1) = 0; %placeholder
d(end + 1) = (-1*h^2)*r(x) + (1 + (h/2)*p(x))*alpha;

for i = 2:N-1
	x = A + i*h;
	a(end + 1) = 2 + q(x)*h^2;
	b(end + 1) = -1 + (h/2)*p(x);
	c(end + 1) = -1 - (h/2)*p(x);
	d(end + 1) = r(x)*(-1*h^2);
end

x = B - h;
a(end + 1) = 2 + q(x)*h^2;
c(end + 1) = -1 - (h/2)*p(x);
d(end + 1) = r(x)*(-1*h^2) + (1 - (h/2)*p(x))*beta;

l = [];
u = [];
z = [];
l(end + 1) = a(1);
u(end + 1) = b(1)/a(1);
z(end + 1) = d(1)/l(1);

for i = 2:N-1
	l(end + 1) = a(i) - c(i)*u(i-1);
	u(end + 1) = b(i)/l(i);
	z(end + 1) = (d(i) - c(i)*z(i-1))/l(i);
end

l(end + 1) = a(N) - c(N)*u(N-1);
z(end + 1) = (d(N) - c(N)*z(N-1))/l(N);

w = zeros(1, N+2);
w(1) = alpha;
w(N+2) = beta;
w(N+1) = z(N);
for i = N-1:-1:1
	w(i+1) = z(i) - u(i)*w(i+2);
end

x = [];
for i = 0:N+1;
	x(end + 1) = A + i*h;
end

clear all




%question 3B
'Newtons iteration'
%I am HARD-CODING the values for the problem
xo = 1/5;
yo = 2/5;
zo = 3/5;
to = 4/5;
f1 = @(x, y, z, t) x + y + z + t - 1
f2 = @(x, y, z, t) x^2 + y^2 + z^2 + t^2 - 1
f3 = @(x, y, z, t) x^3 + y^3 + z^3 + t^3 - 1
f4 = @(x, y, z, t) x^4 + y^4 + z^4 + t^4 - 1

F = zeros(4, 1);
F(1) = f1(xo, yo, zo, to);
F(2) = f2(xo, yo, zo, to);
F(3) = f3(xo, yo, zo, to);
F(4) = f4(xo, yo, zo, to);

%i dont have the symbolic toolbox on my octave,
%weirdly enough, so I can't perform alot of 
%operations - hence I am hard-coding this problem
J = zeros(4, 4);
for i = 1:4 
	J(1,i) = 1;
	J(2,i) = (2*i)/5;
end
J(3,1) = 3/25;
J(3,2) = 12/25;
J(3,3) = 27/25;
J(3,4) = 48/25;
J(4,1) = 4/125;
J(4,2) = 32/125;
J(4,3) = 108/125;
J(4,4) = 256/125;


y = J\(-F);

x = zeros(4, 1);
x(1) = xo;
x(2) = yo;
x(3) = zo;
x(4) = to;

x = x + y;

clear all

%----------

alpha = zeros(1,5);
alpha(1) = 1;
alpha(2) = -0.7;

W = zeros(2, 5);
W(1,1) = alpha(1);
W(1,2) = alpha(2);

h = 0.01;

%plot for F alpha0 and alpha1, because we already have them
%use Adam-bashforth

ab2(f, [0 2.5], 1, 250)




