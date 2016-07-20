function[B, C, D] = CUBICSPSLINE(option, T, X)
%INPUT 
% -->option: an integer 1 or 2; if anything else the function will exit. If option = 1, then the
% interpolation will use the Natural boundary condition. if option = 2, then it will use the 
% clamped boundary condition
% --> T: an n x 1 array (aka a vector) of data points with the values of t (the horizontal axis)
% --> X: an n x 1 array of data points with the values of x(t) (the vertical axis)

%OUTPUT
% B, C and D : each of these are (n-1) x 1 arrays that contain the b'ith, c'ith, and d'ith
% coefficient for the cubic spline interpolation at each of its i'th entries.

	%give feedback on the option selection
	if option == 1
		fprintf('You have chosen to use the Natural Boundary Condition in interpolation.\n');
	elseif option == 2
		fprintf('You have chosen to use the Clamped Boundary Condition in interpolation.\n');
	else
		fprintf('Invalid condition input. Please Enter as argument the integer 1 or 2.\n');
		return;
	end
	
	%initialize the return variables and the local form of the input variables
	B = [];
	C = [];
	D = [];
	
	t = T;
	x = X;
	
	%check if the vector sizes are equal
	[n, cols1] = size(x);
	[m, cols2] = size(t);
	if (n ~= m) | (cols1 ~= cols2)
		fprintf('The two vectors containing the data points have different dimensions.\n');
		return;
	end
	
	%fill up the return vectors
	for i = 1:(n-1)
		B(end + 1) = 0;
		C(end + 1) = 0;
		D(end + 1) = 0;
	end
	
	%H is h_i = t_i+1 - t_i, for i = 0... n-1  -> we initialize and give values to array of h-values used later
	H = [];
	for i = 1:(n-1)
		H(end + 1) = t(i + 1) - t(i);
	end
	
	% M is the solution, A is the tridiagonal matrix, R is the vector values in AM = R, M is second derivative we want
	M = [];
	A = [];
	R = [];
	Mu = [];
	Lmda = [];
	
	%create Mu array 
	if option == 1			%NBC condition
		Mu(end + 1) = 0;
	elseif option == 2		%CBC condition
		Mu(end + 1) = 1;
	else
		fprintf('How did you get this far if the condition option was not 1 or 2?\n');
		return;
	end

	for i = 2:(n - 1)
		Mu(end + 1) = ( H(i) / ( H(i - 1) + H(i) ) );
	end
	
	%create Lmda array 
	for i = 1:(n - 2)
		Lmda(end + 1) = 1 - Mu(i);
	end

	if option == 1			%NBC condition
		Lmda(end + 1) = 0;
	elseif option == 2		%CBC condition
		Lmda(end + 1) = 1;
	else
		fprintf('How did you get this far if the condition option was not 1 or 2?\n');
		return;
	end
	
	%precompute f' using five-point endpoint formula in case of option 2 CBC
	h = (t(2) - t(1));
	dxzero = ((1/(12*h)) * (-25*x(1) + 48*x(2) - 36*x(3) + 16*x(4) - 3*x(5)));
	dxnth = ((1/(12*h)) * (-25*x(n) + 48*x(n - 1) - 36*x(n - 2) + 16*x(n - 3) - 3*x(n - 4)));
	
	%create R array for for x-values
	if option == 1			%NBC condition
		R(end + 1) = 0;
	else					%CBC condition
		R(end + 1) = ( ((6/(H(1)^2))*(x(2) - x(1))) - ((6/H(1))*dxzero) ); 
	end

	for i = 1:(n - 2)
		R(end + 1) = (6 / (H(i + 1) + H(i))) * ( ((x(i+2)-x(i+1))/H(i + 1)) - ((x(i+1) - x(i))/H(i)) );
	end

	if option == 1			%NBC condition
		R(end + 1) = 0;
	else					%CBC condition
		R(end + 1) = ((6/H(n - 1))*dxnth) - (6*(x(n) - x(n - 1))/(H(n -1)^2));
	end
	
	%create the tridiagonal matrix A
	for i = 1:n
		for j = 1:n
			if i == j
				A(i, j) = 2;
			elseif j == (i + 1)
				A(i, j) = Mu(i);
			elseif (i - j) == 1
				A(i, j) = Lmda(j);
			else
				A(i, j) = 0;
			end
		end
	end
	
	%call croutsolver to solve the linear system AM = R for M
	M = CroutSolver(A, R);
	if option == 1		%NBC condition; CBC does not alter M(1) or M(n)
		M(n) = 0;
		M(1) = 0;
	end
	
	%solve for coefficients 
	for i = 1:(n-1)
		B(i) = ((x(i + 1) - x(i))/H(i)) - ((H(i)/6) * (M(i + 1) + 2*M(i)));
		C(i) = (M(i)/2);
		D(i) = (M(i + 1) - M(i))/(6*H(i));
	end

	
end