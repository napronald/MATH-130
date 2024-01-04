clear all
clc
close all
input = [1 0 1 0 1 0 1 0 1]
if input(1) == 0
    sign = 1;
else 
    sign = -1;
end
c = input(2)*2^3+input(3)*2^2+input(4)*2^1+input(5)*2^0
f = input(6)*(1/2)+input(7)*(1/4)+input(8)*(1/8)+input(9)*(1/16)
number_8bit = sign*2^(c)*(1+f)
% f = @(x) exp(6*x)+1.441*exp(2*x)-2.079*exp(4*x)-0.333;
% derf = @(x) 6*exp(6*x)+2.882*exp(2*x)-8.316*exp(4*x);
% derf_2 = @(x) 36*exp(6*x)+5.764*exp(2*x)-33.264*exp(4*x);
% 
% max_iter = 50;
% iter = 1;
% tol = 10^-5;
% start = 0;
% 
% while (iter < max_iter)
%     p = start - (f(start)*derf(start))/(derf(start)^2-f(start)*derf_2(start));
% %     p = start - (2*f(start)*derf(start))/(2*derf(start)^2-f(start)*derf_2(start));
%     err(iter) = abs(p-start);
% 
%     if (err(iter) < tol)
%         root = p
%         break
%     else
%         iter = iter + 1
%         start = p;
%     end
% end
% 
% x = err(1:end-1);
% y = err(2:end);
% logx = log(x);
% logy = log(y);
% 
% figure(1)
% plot(logx,logy,'r-o')
% grid on
% 
% slope =  (logy(end) - logy(end-1)) / (logx(end) - logx(end-1))
