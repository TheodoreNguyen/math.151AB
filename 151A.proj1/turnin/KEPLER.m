function [soln]  = KEPLER(T, e, t, accuracy)
	%INPUTS
	%T, e, t are all 'fixed' coefficients in kepler's equation
	%accuracy is the factor/number of digits for which we want E,
	%the solution to f(E) = 0, accurate to.
	
	%OUTPUTS
	%soln is the solution E' that satisfies f(E') = 0
	
	%create a dummy local variable for accuracy to pass to
	%bisection's btol
	btol = accuracy
	
	%define the Kepler function desired
	f = @(E) ((2*pi*t)/T) - E + e*sin(E)
	
	%define its derivative
	df = @(E) e*cos(E) - 1
	
	%calculate the interval endpoints based on function properties
	Ea = (2*pi*t/T - 1);
	Eb = (2*pi*t/T + 1);
	
	%calculate maximum number of iterations
	imax = (log(Eb - Ea) - log(accuracy))/2;
	imax = ceil(imax);
	
	%call the newton-bisection hybrid equation solver
	soln = HYBRID(f, df, Ea, Eb, btol, accuracy, imax)
end