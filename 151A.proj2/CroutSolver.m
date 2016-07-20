function [m] = CroutSolver(A, r)
	% INPUTS 
	% A is an n x n Tridiagonal matrix, r is an n x 1 vector
	
	% OUTPUTS
	% m is the solution m to the equation Am = r
	
	% obtain dimensions of the matrix A and then check if conditions hold
	[dim1, dim2] = size(A);
	if dim1 ~= dim2
		printf("Error: A is not a Square matrix.\n");
		return;
	elseif dim1 ~= size(r)
		printf("Error: The dimension of A is not equal to that of r.\n")
		return;
	elseif isbanded(A, 1, 1) ~= 1
		printf("Error: A is not Tridiagonal.\n")
		return;
	end
	n = dim1;

	% initialize the L and U matrices and the Z vector and m solution in correct dims.
	L = zeros(n, n);
	U = zeros(n, n);
	Z = [];
	m = [];
	for i = 1:n
		Z(end + 1) = 0;
		m(end + 1) = 0;
	end

	%Set up and solve Lz = b 
	L(1, 1) = A(1, 1);
	U(1, 2) = A(1, 2)/L(1, 1);
	Z(1) = r(1)/L(1, 1);

	for i = 2:(n-1)
		L(i, i - 1) = A(i, i - 1);
		L(i, i) = A(i, i) - (L(i, i - 1) * U(i - 1, i));
		U(i, i + 1) = A(i, i + 1)/L(i, i);
		Z(i) = (r(i) - (L(i, i - 1) * Z(i - 1)))/L(i, i);
	end
	
	L(n, n - 1) = A(n, n - 1);
	L(n, n) = A(n, n) - (L(n, n - 1) * U(n - 1, n));
	Z(n) = (r(n) - (L(n, n - 1) * Z(n - 1)) )/L(n, n);
	
	%Solve Um = Z
	m(n) = Z(n);
	for i = (n-1):-1:1
		m(i) = Z(i) - ( U(i, i + 1) * m(i + 1) );
	end

end