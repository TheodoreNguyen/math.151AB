% set parameters for Kepler's equation

T = 1;
e = 0.25;

% create an array t of 50 elements, where t = [0.01 0.03 0.05 .... 0.97 0.99]
t = [];
iterator = 0.01;
while (iterator < 1);
	t(end + 1) = iterator;
	iterator = iterator + 0.02;
end

%the spec wants 6 digits of accuracy, yeah?
accuracy = 0.000001;

%prepare an array to input corresponding E values to the t values.
E = [];

%for every value of t = 0.01, 0.03, ..., 0.97, 0.99; solve the equation
% f(E) = 2*pi*t/T - E + e*sin(E) = 0 for E, then add it to the E array.
for i = 1:numel(t);
	E(end+1) = KEPLER(T, e, t(i), accuracy);
end