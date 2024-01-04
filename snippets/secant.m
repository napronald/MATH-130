clear all
clc
close all

% f = @(x) cos(x);
% p0 = 0.5;
% p1 = pi/4;


f = @(x) cos(x+sqrt(2))+x*((x/2)+sqrt(2));
% f = @(x) exp(6*x)+3*(log(2))^2*exp(2*x)-(log(8))*exp(4*x)-(log(2))^3;
% f = @(x) exp(6*x)+1.441*exp(2*x)-2.079*exp(4*x)-0.333;

p0 = -1;
p1 = -2;


tol = 10^-5;
max_iter = 50;
iter = 1;

q0 = f(p0)
q1 = f(p1)

while iter < max_iter
   p = p1-q1*(p1-p0)/(q1-q0)
   error(iter) = abs(p-p1);
   if error(iter) < tol
       root = p
       break
   else
       iter = iter + 1
       p0 = p1;
       q0 = q1
       p1 = p;
       q1 = f(p)
   end
end

x = error(1:end-1);
y = error(2:end);
logx = log(x)
logy = log(y)

figure(1)
plot(logx,logy,'r-o')
grid on
slope =  (logy(end) - logy(end-1)) / (logx(end) - logx(end-1))
