% theodore nguyen 704-156-701 math 151B Spr2016 HW5

%problem 1C (NOT REQUIRED, but an illustration is nice, yeah?)
clear all
printf("Question 1C - Householder's Method\n");
A(4,4) = 10
A
A = zeros(4,4)
A(1,1) = 1;
A(1,2) = 2;
A(1,3) = 3;
A(1,4) = 4;
A(2,1) = 2;
A(3,1) = 3;
A(4,1) = 4;
A(2,2) = 4;
A(2,3) = 9;
A(2,4) = 7;
A(3,2) = 9;
A(4,2) = 7;
A(3,3) = 1;
A(3,4) = 1;
A(4,3) = 1;
A(4,4) = 10;
alpha = (29^.5);
r = ( 0.5 * alpha^2 - 0.5 * alpha * 2 ) ^ .5;
w = zeros(4,1);
w(2) = ( 2 - alpha) / (2 * r);
w(3) = 3 / (2 * r);
w(4) = 4 / (2 * r);
W = w(2:4);
P = eye(3,3) - 2*W*W'/(norm(W)^2);
Pw = eye(4,4);
Pw(2:4,2:4) = P;
PxA = Pw * A
PxAxP = Pw*A*Pw



%problem 3:
clear all
printf("Question 3 - QR algorithm\n");
A = zeros(2,2);
A(1,1) = 2;
A(1,2) = 1;
A(2,1) = 1;
A(2,2) = 0.5;

N = 100;

%type 1: store all previous iterations result in a cell.

T{1} = A;
U{1} = eye(2,2);

for k = 1:N
	Q{k} = zeros(2,2);
	R{k} = zeros(2,2);
	[Q{k},R{k}] = qr(T{k});
	T{k+1} = R{k}*Q{k};
	U{k+1} = U{k}*Q{k};
end
T{101};

%type 2: overwrite previous iterations

T = A;
U = eye(2,2);
Q = zeros(2,2);
R = zeros(2,2);

for k = 1:N
	[Q, R] = qr(T);
	T = R*Q;
	U = U * Q;
end
T