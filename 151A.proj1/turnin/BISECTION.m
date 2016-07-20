function [ A, B, p ] = BISECTION( f, a, b, TOL, imax )
%INPUTS
%f is the function we want to find a zero of
%a is the starting position of the initial interval we perform BISECTION
%b is the ending position of the initial interval we perfrorm BISECTION
%TOL is the acceptable tolerance for which our solution can be off by 
%imax is the maximum number of iterations allowed 

%OUTPUTS
%A is the starting position of the final interval, while B is the ending
%the final interval is that which we computed the zero with/ended up with
%p is the x-value of the position of the zero of f in the interval
    i = 1;		%we initialize the iterator and the function value of the left endpoint
    FA = f(a);				
	p = a;		%give return values a predefined value
    A = a;
    B = b;
	
	if( (f(a) * f(b)) > 0 ) %check if initial endpoints satisfy bisection method reqs
		printf "The function has the same sign at both inputted endpoints.\n"
		return 
	end
	
    while (i <= imax)
        p = (A + (B - A)/2)		%set p to be the midpoint of interval [A,B]
        FP = f(p);				%evaluate the function at p
        if (FP == 0) 			%if f(p) = 0, then p is a zero and we are done
            return
		else if( ((B - A)/2) < TOL )	%if half the current interval is smaller than the 
			return						%tolerance, p is close enough to the x such that 			
        end								%f(x) = 0 for this tolerance
		
        i = i + 1;			%ready for next iteration
		
        if((FA * FP) > 0)	%if fcn value at left endpoint is the same sign as the the 
            A = p;			%fcn value at the current estimate p, then the true zero is 
            FA = FP;		%right of p; set the new left endpoint as p
        else
            B = p;			%otherwise, from ^^, then true zero is left of p, set the 
        end					%new right endpoint as p
    end
    printf "BISECTION: Max number of iterations reached. Exiting with last computed solutions.\n"
end

