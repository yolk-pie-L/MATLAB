% -----------------------------------------
% HW1
% LiXin
% 2022/2/18
% -----------------------------------------
clear all; close all; clc

%% Problem 1

% 1
r=0.15; % interest rate
L=50000; % loan amount
N=[0.5:1/12:20]; % number of years
PN=(r*L*(1+r/12).^(12*N))./(12*((1+r/12).^(12*N)-1)) % monthly payment

% 2
plot(N, PN) % plot the monthly payment and the number of years

% 3
text(10,1000,"LiXin")

%% Problem 2
clear all; close all; clc % erase all the workspace data, command window output, and close all figures
% a
A=[20 4 2 6; 6 37 2 3;8 5 9 9] % create the matrix
% b
x1=A(1,:) % assign the first row of A to a vector called x1
% c
y=A([end-1 end],:) % assign the last 2 rows of A to an array called y
% d
B=A(:,2:2:size(A,2)) % assign the even-numbered columns of A to an array called B
% e
C=A' % assign the transpose of A to C
% f
reciprocal = 1./A % compute the reciprocal of each element of A
% g
A(3, 2)=100 % change the number in column 2, row 3 of A to 100

%% Problem 3
clear all; close all; clc % erase all the workspace data, command window output, and close all figures
% a
figure % create an empty figure

% b
t=0:0.01:2*pi; % create a vector contains numbers range from 0 to 2*pi
y=sin(4*t).*cos(2*t); % calculate the result of the function according to the scope
ax1=subplot(1,2,1);
plot(t,y) % normal plot
ax2=subplot(1,2,2);
polarplot(t,y) % polar plot

% c
legend(ax1,'y^k_{max}=sin(4t)cos(2t)') % put legend on the first subgraph

% d
grid(ax1, 'on') % grid on the first subgraph

% e
text(0,0,'Lixin') % text on the second subgraph

%% Problem 4
clear all; close all; clc % erase all the workspace data, command window output, and close all figures
% a
num=input('Enter a number: ');

% b
num_cm=num * 2.54; % convert to cm
fprintf('%.2f inches is %.2f cm\n', num, num_cm)

% c
num_mm=num_cm * 10; % the number converted to mm
formatSpec='%.2f';
str = [num2str(num, formatSpec), ' is also ', num2str(num_mm, formatSpec), 'mm'];
disp(str)

%% Problem 5 
clear all; close all; clc % erase all the workspace data, command window output, and close all figures
Nr = logspace(4,8,100); % Reynolds number
for De = [20 100 1000 1000 100000] % D/epsilon
    f = 0.25./(log(1/(3.7*De)+5.74./Nr.^0.9)/log(10)).^2;
    loglog(Nr, f)
    hold on
end
grid on
grid minor
Nr=[1e2:0.2e4];
f = 64./Nr;
plot(Nr, f) % plot f = 64 / Nr
title('Moody''s Diagram')
ylabel('Friction Factor')
xlabel('Reynolds Number N_R')
legend('D/\epsilon = 20', 'D/\epsilon = 100', 'D/\epsilon = 1000', 'D/\epsilon = 10000', 'D/\epsilon = 100000','Laminar flow',...
'Location','southwest')
text(1e7,0.08,'Lixin') % print my name
xlim([0.8e3 1e8])
ylim([0.8e-2 1e-1])% adjust axis limits

%% Problem 6
clear all; close all; clc % erase all the workspace data, command window output, and close all figures
[x, y]=meshgrid(-2:0.1:2);
f=50*y.^2.*exp(-x.^2-0.5*y.^2);
C=x.*y;
surf(x,y,f,C)
xlabel('x');
ylabel('y');
zlabel('f(x,y) = 50y^2e^{-x^2-0.5y^2}');
title('My Plot Title')


