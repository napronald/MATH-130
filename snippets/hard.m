clear all
clc
close all

f = @(x) sin(x);
x0 = 0;
xn = 2*pi;
N = 5;
data_x = linspace(x0,xn,N)
data_y = f(data_x)

data_z = linspace(x0,xn,N)

n = length(data_x);
N = zeros(n,n);
N(:,1)=data_y;
for i = 2:n
    for j = i:n
        N(j,i) = (N(j,i-1) - N(j-1,i-1)) / (data_x(j) - data_x(j-i+1))
    end
end
coefficents = diag(N);
pn = coefficents(1);

ai=1;
syms x 

% x = sym(sym(data_z))
x = data_z
for j=2:n
    for k=1:j-1
        ai=ai.*(x-data_x(k));
    end
    pn =pn+coefficents(j).*ai;
end
pn





