clear all
clc
close all

%Discussion Topic 1: Finite precision (fp), round-off error,
%and ways to measure accuracy.

%Today we will be finding out which is the best method for
%evaluating the derivative of a function at a point. 

%1). the function f and the exact derivative
f =@(x) x.^2;
fprime =@(x) 2.*x;

%2). evaluate the derivative at point val
val = 2;

%3). below is an algorithm for evaluating the derivative 
%of a function at a point. The method depends on the 
%parameter h. 
h = 0.25;
fapprox1 = ( f(val+h) - f(val) ) ./ h

%4). Questions:
%   a). What is a good value of h for approximating f'(val)?
% The smaller the h the better
%   b). Is is more accurate to use smaller or larger values of h?
% smaller
%   c). Compute the error (abs & relative) for 5 different
%         h values. 
%   d). Plot the h vs. error (abs & relative, two seperate plots) using log scale. 
%   e). ** Is there a point where small h will lead to larger errorr?

%5). Here is another method for approximating the derivative.
%Note: this is using the same h as method 1.
fapprox2 = ( f(val+h) - f(val-h) ) ./ (2.*h)

%6). Questions:
%   a). For the same value of h, which method is more accurate?
%   b). Do question c from above, use the same 5 values of h.
%   c). Plot the h vs. error for method 2 in the same plot
%       as method 1. Remember to use log scale.
%       Note, you will have 2 plots. 1 for abs error with both methods,
%       and 1 for relative error with both methods.
%   d). How can you use the plot to determine which method is better?

expected = fprime(val);

% rel_1 = abs(fapprox1-expected)
% abs_1 = abs(fapprox1-expected) / abs(expected)
% 
% rel_2 = abs(fapprox2-expected)
% abs_2 = abs(fapprox2-expected) / abs(expected)

figure(1)
for h=1:5
    fapprox1 = ( f(val+h) - f(val) ) ./ h;
    rel_1 = abs(fapprox1-expected);
    abs_1 = abs(fapprox1-expected) / abs(expected);
    plot(h,rel_1,".-",'MarkerSize',10)
    plot(h,abs_1,".-",'MarkerSize',10)
    xlabel("h")
    ylabel("relative error")
    hold on
end
figure(2)
for h=1:5
    fapprox2 = ( f(val+h) - f(val) ) ./ h;
    rel_2 = abs(fapprox2-expected);
    abs_2 = abs(fapprox2-expected) / abs(expected);
    plot(h,rel_2,".-",'MarkerSize',10)
    plot(h,abs_1,".-",'MarkerSize',10)
    xlabel("h")
    ylabel("relative error")
    hold on
end







