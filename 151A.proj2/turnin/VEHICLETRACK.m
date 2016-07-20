function[Bxi, Cxi, Dxi, Byi, Cyi, Dyi, Bxii, Cxii, Dxii, Byii, Cyii, Dyii] = VEHICLETRACK(file)
%INPUT
% ---> file: this is a type string that is the name of the data file containing the n x 3 data 
% matrix to be loaded

%OUTPUT
% 12 arguments, grouped into FOUR BCD trios. 
% each B, C, and D trio are each vectors containing the coefficients to the cubic spline 
% interpolation.
% if the trio has an "x" in the name, then it is interpolating between x and t 
% if the trio has a "y" in the name, then it is interpolating between y and t
% if the trio has one "i" in the name, then it used natural boundary condition in interpolating
% if the trio has two "ii" in the name, then it used clamped boundary condition in interpolating

	structdata = load(file);			%extract the data from the file into the structure
	Data = file.ip;					%extract the n x 3 matrix from the structure into a matrix
	[n, cols] = size(Data);		%extract the number of rows (and cols) of the matrix
	
	%initialize data vectors to consist of the raw data points corresponding to t, x(t), y(t)
	t = []			
	x = []
	y = []
	
	%fill up those data vectors with the data points.
	for i = 1:n
		t(end + 1) = Data(i, 1);
		x(end + 1) = Data(i, 2);
		y(end + 1) = Data(i, 3);
	end
	
	%We now call CUBICSPLINE to get all the coefficients for the splines... there's alot
	
	%interpolating X and T using NATURAL BOUNDARY CONDITION
	[Bxi, Cxi, Dxi] = CUBICSPLINE(1, t, x);
	
	%interpolating Y and T using NATURAL BOUNDARY CONDITION
	[Byi, Cyi, Dyi] = CUBICSPLINE(1, t, y);
	
	%interpolating X and T using CLAMPED BOUNDARY CONDITION
	[Bxii, Cxii, Dxii] = CUBICSPLINE(2, t, x);
	
	%interpolating Y and T using NATURAL boundary condition
	[Byii, Cyii, Dyii] = CUBICSPLINE(2, t, y);
	
	
end