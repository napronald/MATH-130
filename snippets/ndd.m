clear all 
clc
close all

data_x = [-1 0 1];
data_y = [2 1 2];

[a] = NewtonDivDiff(data_x,data_y)


N = 3;
x0 = -2;
xn =  2;
Z = linspace(x0,xn,N);

% polynz = NewtonInterpPolyn(N,data_x,Z)

p = @(x)  x.^2+1; 
data_z = [1 1 1.5];


N = 5;
x0 = 0;
xn =  2*pi;

f = @(x) sin(x);
% f = @(x) 1./(x.^2+1)
Z = linspace(x0,xn,N)
fZ = feval(f,Z) 


 
% Used NewtonDivDiff(Z,fZ)
p_sin = @(x) 0.6369426751592356.*(x-0) + -0.4056959714390036.*(x-0).*(x-pi/2) + 0.08613502578322793.*(x-0).*(x-pi/2).*(x-pi) + 4.419677645800782e-18.*(x-0).*(x-pi/2).*(x-pi).*(x-3.*pi/2)
pnZ = p_sin(Z)

% p_n2 = @(x) 0.2 + 0.3.*(x+2) + 0.1.*(x+2).*(x+1) + -0.19999999999999998.*(x+2).*(x+1).*(x-0) + 0.09999999999999999.*(x+2).*(x+1).*(x-0).*(x-1)
% pnZ = p_n2(Z)

subplot(2,1,1)
plot(Z,fZ,'b-',Z,pnZ,'rx-');
xlabel('z');
ylabel('f(z), P_n(z)')
titletxt = sprintf('f(z) = --,\t P_n(z) = xxx,\t n+1 = %d',N);
title(titletxt)

subplot(2,1,2)
plot(Z,abs(fZ-pnZ),'.-');
xlabel('z')
ylabel('|f(z)-P_n(z)|')
title('Absolute Error = |f(z)-P_n(z)|')

function a = NewtonDivDiff(data_x,data_y)
    n = length(data_x);
    N = zeros(n,n);
    N(:,1) = data_y;
    for i = 2:n
        for j = i:n
            N(j,i) = (N(j,i-1) - N(j-1,i-1)) / (data_x(j) - data_x(j-i+1));
        end
    end
    for i = 1:n
       a(i) = N(i,i);
    end
end


% 
% function polynz = NewtonInterpPolyn(N,data_x,Z)
%     N = 3;
%     poly_terms = ones(length(Z),length(N)-1)
%     for j = 1:N-1
%         for k = 1:j
%             poly_terms(:,j) = poly_terms(:,j)*(z-x(k))'
%         end
%     %multiply by coefficent
%     poly_terms(:,j) = a(1,j+1)*poly_terms(:,j)
%     end
% end

clear all
clc
close all

f = @(x) sin(x);
x0 = 0;
xn = 2*pi;

% f = @(x) 1./(x.^2+1);
% x0 =-2;
% xn=2;

N = 5;
num_predicted_points = 60

data_x = linspace(x0,xn,N)
data_y = f(data_x)

data_z = linspace(x0,xn,num_predicted_points) 

n = length(data_x);
N = zeros(n,n);
N(:,1)=data_y;
for i = 2:n
    for j = i:n
        N(j,i) = (N(j,i-1) - N(j-1,i-1)) / (data_x(j) - data_x(j-i+1));
    end
end
coefficents = diag(N);
pn = coefficents(1);

syms x 

% x = sym(sym(data_z))
x = data_z
for j=2:n
    ai=1;
    for k=1:j-1
        ai=ai.*(x-data_x(k));
    end
    pn =pn+coefficents(j).*ai;
end
pn

% f(data_z)
% pn(3)

subplot(2,1,1)
plot(data_x,f(data_x),'bo-',data_z,pn,'rx-');
xlabel('z');
ylabel('f(z), P_n(z)')
titletxt = sprintf('f(z) = --,\t P_n(z) = xxx,\t n+1 = %d',N);
title(titletxt)
legend('Predicted points','Real Function')


subplot(2,1,2)
plot(data_z,abs(f(data_z)-pn),'.-');
xlabel('z')
ylabel('|f(z)-P_n(z)|')
title('Absolute Error = |f(z)-P_n(z)|')
legend('Error')
