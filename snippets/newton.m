clear all
clc
close all

f = @(x) x.^3 + 4*x.^2 -10;
derf = @(x) 3*x.^2 +8*x; 

% f = @(x) x.^2-4;
% derf = @(x) 2*x;
% start = 1;

f = @(x) cos(x+sqrt(x+sqrt(2)))+x*((x/2)+sqrt(2));
derf = @(x) x-sin(x+sqrt(2))+sqrt(2);

start = -1.5;
tol = 10^-5;
max_iter = 50;

[output, err_est, iter] = newton(f,derf,start,tol,max_iter)

x = err_est(1:end-1);
y = err_est(2:end);
logx = log(x)
logy = log(y)

figure(1)
plot(logx,logy,'r-o')
grid on

figure(2)
num_iter = 1:1:iter;
plot(num_iter,err_est,'b-o')
xlabel('number of iterations')
ylabel('error')
grid on

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
            start = p
        end
    end
end 