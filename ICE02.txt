% -----------------------------------------
% HW1
% LiXin
% 2022/2/28
% -----------------------------------------
close all; clear all; clc;

%% Input for Problem 1~5

% Problem 1
M = [1 2; 3 4];
even_index(M)

% Problem 2
v = [1 2 3 4];
flip_it(v)

% Problem 3
N = magic(5);
n = 2;
top_right(N, n);

% Problem 4
A = magic(5);
peri_sum(A)

% Problem 5
cal(pi/2)

%% Problem 6
h = 2;
v = sqrt(2 * 9.8 *h);
for i = 1:8 % rebound 8 times
    v = 0.85 * v;
end
h = v*v/(2*9.8)

%% Problem 7
n=1;
T=300;
R=0.08206;
a=1.39;
b=0.039;    
V=linspace(0.08,6,100);
P1=n*R*T./V;
P2=n*R*T./(V-n*b)-n*n*a./V.^2;
plot(V,P1)
hold on
plot(V,P2)

%% Problem 8

% a
x = 10 * rand(ceil(10*rand)+2,1)
% b
mysum=0;
for i = 1:size(x,1)
    mysum = mysum + x(i,1);
end
mysum
% c
if mysum == sum(x)
    disp('Congratulations!! you did it right')
    load handel;
    sound(y, Fs)
else
    fprintf('Sorry, %.2f ~= %.2f. Please try again.\n', mysum, sum(x))
end
% d
x = 10 * rand(ceil(10*rand)+2,1);
mysum=0;
i = 1;
while i <= size(x,1)
    mysum = mysum + x(i,1);
    i = i + 1;
end
mysum
if mysum == sum(x)
    disp('Congratulations!! you did it right')
    load handel;
    sound(y, Fs)
else
    fprintf('Sorry, %.2f ~= %.2f. Please try again.\n', mysum, sum(x))
end
    

%% Problem 1
% returns a matrix that contains only those elements of M that are in
% even rows and columns
function res = even_index(M) % M as a matrix
    res = M(2:2:end, 2:2:end);
end

%% Problem 2
% returns the opposite order of v
function w = flip_it(v) % row vector v, row vector w
    for i = size(v,2) : -1 : 1
        w(1,size(v,2) - i + 1) = v(1,i);
    end
end

%% Problem 3
% returns the n-by-n square subarray of N located at the top right corner
% of N
function M = top_right(N, n) % a matrix N, a scalar non-negative integer n
    M = N(1:n, end-n+1:end)
end

%% Problem 4
% add together the elements that are in the first and last rows and columns
function my_sum = peri_sum(A)
    my_sum = 0;
    my_sum = my_sum + sum(A(:,1)) + sum(A(:,end)) + sum(A(1,:)) + sum(A(end,:));
    my_sum = my_sum - A(1,1) - A(1,end) - A(end, end) - A(end,1);
end

%% Problem 5
% power series for sin(x), 5 terms are needed, when i exceeds 20, the loop 
% terminates
function sin_x = cal(x)
sin_x = 0;
    for i = [1:2:20]
        sin_x = sin_x + (-1)^floor(i/2)*x^i/factorial(i);
        sin(x)-sin_x
    end
end