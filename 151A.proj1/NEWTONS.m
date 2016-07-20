function [p, flag] = NEWTONS (f, df, xo, TOL, imax)
%INPUTS
%f is the function, df is the derivative of f, xo is the initial guess
%TOL is the tolerance, imax is the maximum number of iterations

%OUTPUTS
%p is the approximated x such that f(x) = 0
%flag is an indicator = 1 if newton's method converges, 0 otherwise
	
	i = 1;		%initialize the iterator, flag=0 cuz hasnt converged
	po = xo;	%let po be the local version of the approx x
	flag = 0	%initialize p to be the initialized approx
	p = po;
	while (i <= imax)
		
		if(df(po) == 0)		%if f'(po)=0, newtons method diverges b/c cant divide 0
			printf "f'(po) = 0; cannot divide by zero in this iteration. Exiting.\n"
			return
		end
		p = po - (f(po)/df(po));	%compute the x-intercept of the tangent line
		
		if (abs(p - po) < TOL) 		%if the distance between this x-approx and the
			flag = 1;				%last x-approx is less than the tolerance
			return;					%then the soln is acceptable for this tolerance
		end
		i = i + 1;					%otherwise, prepare for next iteration
		po = p;						%incrementing i and setting up new po
	end
	printf "NEWTONS: Max number of iterations reached. Exiting with last computed solutions.\n"
end