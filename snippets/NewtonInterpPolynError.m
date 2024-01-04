% fname = @(x) sin(x);
% x0 = 0;
% xn = 2*pi
% N = 5
% NewtonInterpPolynError(fname, x0, xn, N)
% 
% function dummy = NewtonInterpPolynError(fname, x0, xn, N) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab code to plot f(x), its interpolating polynomial P_n(x) and 
% the absolute error = |f(x) - P_n(x)|
% INPUT:  fname  =  the name of the function f(x) to interpolate 
%         x0     =  lower bound of the interval on which to interpolate 
%                   (also the first interpolating node) 
%         xn     =  upper bound of the interpolation interval 
%                   (also the last interpolating node) 
%         N      =  number of interpolation nodes (N = n+1)        
% OUTPUT: plots of f(z), P_n(z) and absolute error = |f(z) - P_n(z)|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a vector of Z values   
% Z = linspace(0,2*pi,5);  
   
% Compute f(Z)   
% fZ = feval(fname,Z);   
   
% Compute P_n(Z) on nodes X = [x_0,...,x_n]  
% X = linspace(x0,xn,N);
% pnZ = NewtonInterpPolyn(N,X,fname,Z);
% 
% % Plot f(z), P_n(z) and absolute error = |f(z) - P_n(z)|
% subplot(2,1,1)
% plot(Z,fZ,'-',Z,pnZ,'x');
% xlabel('z');
% ylabel('f(z), P_n(z)')
% titletxt = sprintf('f(z) = --,\t P_n(z) = xxx,\t n+1 = %d',N);
% title(titletxt)
% 
% subplot(2,1,2)
% plot(Z,abs(fZ-pnZ),'.');
% xlabel('z')
% ylabel('|f(z)-P_n(z)|')
% title('Absolute Error = |f(z)-P_n(z)|')
% 
% dummy = 1;
% return
% 
% end
% 


% p = @(x)  x.^2+1; 
% data_z = [1 1 1.5];

% NewtonInterpPolyn(p,data_z)

% N = 5;
% x0 = -2;
% xn =  2;
% 
% f = @(x) sin(x);
% % f = @(x) 1./(x.^2+1)
% Z = linspace(x0,xn,N)
% fZ = feval(f,Z) 
% 
% 
%  
% % Used NewtonDivDiff(Z,fZ)
% p_sin = @(x) 0.6369426751592356.*(x-0) + -0.4056959714390036.*(x-0).*(x-pi/2) + 0.08613502578322793.*(x-0).*(x-pi/2).*(x-pi) + 4.419677645800782e-18.*(x-0).*(x-pi/2).*(x-pi).*(x-3.*pi/2)
% pnZ = p_sin(Z)
% 
% % p_n2 = @(x) 0.2 + 0.3.*(x+2) + 0.1.*(x+2).*(x+1) + -0.19999999999999998.*(x+2).*(x+1).*(x-0) + 0.09999999999999999.*(x+2).*(x+1).*(x-0).*(x-1)
% % pnZ = p_n2(Z)
% 
% subplot(2,1,1)
% plot(Z,fZ,'b-',Z,pnZ,'rx-');
% xlabel('z');
% ylabel('f(z), P_n(z)')
% titletxt = sprintf('f(z) = --,\t P_n(z) = xxx,\t n+1 = %d',N);
% title(titletxt)
% 
% subplot(2,1,2)
% plot(Z,abs(fZ-pnZ),'.-');
% xlabel('z')
% ylabel('|f(z)-P_n(z)|')
% title('Absolute Error = |f(z)-P_n(z)|')

% clear all
% clc
% close all
% syms x;
% 
% data_x = [0 2 4];
% data_y = [1 5 17];
% n = length(data_x);
% N = zeros(n,n);
% N(:,1)=data_y;
% for i = 2:n
%     for j = i:n
%         N(j,i) = (N(j,i-1) - N(j-1,i-1)) / (data_x(j) - data_x(j-i+1))
%     end
% end
% coefficents = diag(N);
% pn = coefficents(1);
% 
% 
% ai=1;
% for j=2:n
%     for k=1:j-1
%         ai=ai*(x-data_x(k));
%     end
%     pn =pn+coefficents(j)*ai;
% end
% pn


clear all
clc
close all
syms x;

xi = [-1 0 1 2];
yi = [2 1 2 -7];
n = length(xi);
m=n+1;
M=zeros(n,m);
M(:,1)=xi;
M(:,2)=yi;
for j=3:m
    for i=j-1:n
        M(i,j) =(M(i,j-1)-M(i-1,j-1))/(M(i,1)-M(i-j+2,1))
    end
end
S=M(:,2:m)
d = diag(S)
pn=d(1)

for k=2:n
    ai=1;
    for i=1:k-1
        ai=ai*(x-xi(i));
    end
    pn =pn+d(k)*ai;
end
pn;


x = [-1 0 1];
y = [2 1 2];

NewtonDivDiff(x,y)

function NewtonDivDiff(x,y)
    n = length(x);
    M = zeros(n,n);
    M(:,1) = y;
    for i = 2:n
        for j = i:n
            M(j,i) = (M(j,i-1) - M(j-1,i-1)) / (x(j) - x(j-i+1))
        end
    end
end