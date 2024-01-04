clear all
clc
close all

% f = @(x) x.^2-4;
% derf = @(x) 2*x;

f = @(x) cos(x+sqrt(2))+x*((x/2)+sqrt(2));
derf = @(x) x-sin(x+sqrt(2))+sqrt(2);

% f = @(x) exp(6*x)+3*(0.693147)^2*exp(2*x)-(2.079441)*exp(4*x)-(0.693147)^3
% derf = @(x) 6*exp(2*x)*(0.693147)^2-4*exp(4*x)*(2.079441)+6*exp(6*x)

% f = @(x) exp(6*x)+1.441*exp(2*x)-2.079*exp(4*x)-0.333;
% derf = @(x) 6*exp(6*x)+2.882*exp(2*x)-8.316*exp(4*x);

start = -1;
tol = 10^-5;
max_iter = 50;

[root, err_est, iter] = newton(f,derf,start,tol,max_iter)

x = err_est(1:end-1);
y = err_est(2:end);
logx = log(x);
logy = log(y);

figure(1)
plot(logx,logy,'r-o')
grid on

slope =  (logy(end) - logy(end-1)) / (logx(end) - logx(end-1))

function [root, err, iter] = newton(f,derf,start,tol,max_iter)
    iter = 1;
    while (iter < max_iter)
        p = start - (f(start)/derf(start))
        
        err(iter) = abs(p-start);

        if (err(iter) < tol)
            root = p;
            break
        else
            iter = iter + 1;
            start = p;
        end
    end
end 



% p0 = -2;
% p1 = -1;
% 
% tol = 10^-5;
% max_iter = 50;
% [root, err_est, iter] = secant(f,tol,max_iter,p0,p1)
% 
% x = err_est(1:end-1);
% y = err_est(2:end);
% logx = log(x)
% logy = log(y)
% 
% figure(1)
% plot(logx,logy,'r-o')
% grid on
% slope =  (logy(end) - logy(end-1)) / (logx(end) - logx(end-1))

% function [root, err, iter] = secant(f,tol,max_iter,p0,p1)
%     iter = 1;
%     q0 = f(p0);
%     q1 = f(p1);
%     while iter < max_iter
%        p = p1-q1*(p1-p0)/(q1-q0)
%        err(iter) = abs(p-p1);
%        if err(iter) < tol
%            root = p
%            break
%        else
%            iter = iter + 1
%            p0 = p1;
%            q0 = q1;
%            p1 = p;
%            q1 = f(p);
%        end
%     end
% end

