function [p] = HYBRID (f, df, a, b, tolb, toln, nmax)
%INPUTS
%f is the desired function to find root of, df is the derivative of f
%a,b is for the initial interval [a,b] to perform the bisection method 
%tolb and toln are the tolerances for bisection and newton respectively
	%NOTE: tolb is overwritten immediately anyway, so I don't even know
	%why its passed to this function.
%nmax is the maximum number of iterations for newton's method

%OUTPUTS
%p is the approximated x-value such that f(p) = 0 through the method
	
	if(f(a) * f(b) > 0)
		printf "Error: endpoints are not valid\n"	%checks if feeding
		return									%valid endpts bisection
	end
	error = abs(b - a);		%initialize the interval for the bisection
	tolb = error/2;			%method, as well as its tolerance. 
	l = a;
	r = b;
	flag = 0;				%set the flag to zero, as no convergence yet
	while (flag == 0)
		[l, r, p] = BISECTION(f, l, r, tolb, 1)		%one step bisection
		
		%use the approximation from bisection as your initial approximation
		%for newtons for each iteration. This way, bisection will appropriately
		%run once per iteration of this while loop to surely converge, but
		%only does so if the approximations fed to Newton's doesn't already 
		%cause Newtons to converge
		[p, flag] = NEWTONS(f, df, p, toln, nmax)

		if (flag == 1)		%Newtons either converged, or the input from				
			return;			%bisection into Newtons provided a success
		else
			tolb = tolb/2	%if no convergence on this iteration, halve/decrease
		end					%the tolerance fed to bisection: This way, each single
	end						%iteration of bisection must adhere to a more precise 
end							%result, so it should produce a more accurate approx 
							%for newtons, avoiding passing on the same approx