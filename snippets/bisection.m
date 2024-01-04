clear all
clc

f = @(x) x^3;
a = -3;
b = 1;

% f = @(x) sin(x);
% a = -pi/2;
% b = 3*pi/4;

% f = @(x) x^5+5*x^3-x^2+1;
% a = -1;
% b = 2;

% f =@(x) x.^2 - 4;
% a = -3;
% b= 3; 

tolerance = 10^-5;
if (f(a)*f(b))<0
    p = (a+b)/2; 
    while abs(f(p)) > tolerance
        p = (a+b)/2;
        if f(p) > 0
            b = p
            error = abs(p)
        elseif f(p) < 0 
            a = p
            error = abs(p)
        end
    end
    p
else
    fprintf('\n No root found %\n')
end

