% theodore nguyen 704-156-701 math 151B Spr2016

% question 1e
printf("Question 1E - steepest descent method\n");

A = zeros(2,2);
A(1,1) = 1;
A(1,2) = 2;
A(2,1) = 2;
A(2,2) = 3;
b = zeros(2,1);
b(1) = 3;
b(2) = 5;
J = @(x) (1/2)*(norm((A*x) - b))^2;
x0 = [2 3]';
htop = norm(A'*A*x0 - A'*b)^2;
hbot = norm(A*A'*A*x0 - A*A'*b)^2;
h = htop/hbot;
x = x0;
TOL = 10^-9;
for i = 1:2
	g1 = J(x);
	z = x;
	z0 = norm(z);
	if(z0 == 0)
		printf("zero gradient\n");
		break;
	end
	z = z / z0;
	alpha1 = 0;
	alpha3 = 1;
	g3 = J(x - alpha3*z);
	while(g3 >= g1)
		alpha3 = alpha3/2;
		g3 = J(x - alpha3*z);
	end
	
	alpha2 = alpha3 / 2;
	g2 = J(x - alpha2*z);
	h1 = (g2 - g1) / alpha2;
	h2 = (g3 - g2) / (alpha3 - alpha2);
	h3 = (h2 - h1) / alpha3;
	alpha0 = 0.5*(alpha2 - h1/h3);
	g0 = J(x - alpha0*z);
	alphaArr = [alpha0 alpha1 alpha2 alpha3];
	gArr = [J(x - alpha0*z) J(x - alpha1*z) J(x - alpha2*z) J(x - alpha3*z) ];
	[M, I] = min(gArr);
	alpha = alphaArr(I);
	x = x - alpha*z;
end
g1
x

clear all;


% question 2c - power method

printf("Question 2C - Power Method\n");
A = zeros(2,2);			%matrix A
A(1,1) = 2;
A(2,2) = 2;
A(1,2) = 1;
A(2,1) = 1;
x0 = [1 0]';			%initial vector x
n = 2;					%dimensions
TOL = 10^-9;			%epsilon
N = 1000;				%max num iterations

Mu = 0;
x = x0;

[M,p] = max(abs(x));
x = x / x(p);
for k = 1:N
	y = A*x;
	Mu = y(p);
	[M,p] = max(abs(y));
	if y(p) == 0
		printf("A has eigenvalue 0\n");
		break;
	end
	ERR = max(abs(x - (y/y(p))));
	x = y / y(p);
	if ERR < TOL
		printf("Current error is less than tolerance\n");
		break;
	end
	if k == N
		printf("max num iterations exceeded\n");
		break;
	end
end
printf("eigenvalue is Mu\n");
Mu
printf("eigenvector is x\n");
x

clear all

% question 3b

printf("Question 3B - continuation method\n");
xo = 1/5;
yo = 2/5;
zo = 3/5;
to = 4/5;
f1 = @(x, y, z, t) x + y + z + t - 1;
f2 = @(x, y, z, t) x^2 + y^2 + z^2 + t^2 - 1;
f3 = @(x, y, z, t) x^3 + y^3 + z^3 + t^3 - 1;
f4 = @(x, y, z, t) x^4 + y^4 + z^4 + t^4 - 1;

Fp = zeros(4,1);
Fp(1) = f1(xo,yo,zo,to);
Fp(2) = f2(xo,yo,zo,to);
Fp(3) = f3(xo,yo,zo,to);
Fp(4) = f4(xo,yo,zo,to);

%hard-code initial jacobian

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

J1 = @(v) 1;
J2 = @(v) 2*v;
J3 = @(v) 3*v^2;
J4 = @(v) 4*v^3;

h = 0.01;
N = 100;
n = 4;
b = -1*h*Fp;
x = [ xo yo zo to]';
for i = 1:N
	A = J;
	k1 = A\b;
	for j = 1:4 %thru cols
		J(1,j) = 1;
		J(2,j) = J2(x(j) + 0.5*k1(j));
		J(3,j) = J3(x(j) + 0.5*k1(j));
		J(4,j) = J4(x(j) + 0.5*k1(j));
	end
	A = J;
	k2 = A\b;
	for j = 1:4 %thru cols
		J(1,j) = 1;
		J(2,j) = J2(x(j) + 0.5*k2(j));
		J(3,j) = J3(x(j) + 0.5*k2(j));
		J(4,j) = J4(x(j) + 0.5*k2(j));
	end
	A = J;
	k3 = A\b;
	for j = 1:4 %thru cols
		J(1,j) = 1;
		J(2,j) = J2(x(j) + k3(j));
		J(3,j) = J3(x(j) + k3(j));
		J(4,j) = J4(x(j) + k3(j));
	end
	A = J;
	k4 = A\b;
	x = x + (k1 + 2*k2 + 2*k3 + k4)/6;
end
printf("The continuation method approx with h = 0.01 is x\n");
x

clear all

