clear all
clc
close all 

%%No interval possible, |g'(x)| > 1
% g = @(x) x.^5+5*x.^3-x.^2+1+x 


%%No interval guaranteed by convergence theorem
% g = @(x) (-5*x^3+x^2-1)^1/5

%%No interval guaranteed by convergence theorem
% g = @(x) -x^7-6*x^5+x^4-5*x^3+x-1


%%Converges on [1.213,2] guaranteed by convergence theorem
% g = @(x) x^4-8*x^3+24*x^2-32*x+16+x 

%%guaranteed by convergence theorem CONVERGES [-1,2]
% g = @(x) ((x^4+24*x^2-32*x+16)/8)^(1/3) 

%%CONVERGES on [-1,2] guaranteed by convergence theorem
% g = @(x) -16/(x.^3-8*x.^2+24*x-32) 

g = @(x) x+0.5*x*(x.^3-9)*(x.^2+1);
x0 = 2.08;

% x0 = -0.966;
tol = 10^-6;
max_iter = 50;

fixed_pt_method(g,x0, tol, max_iter);

function root = fixed_pt_method(g, x0, tol, max_iter)
    
    iter = 0;
    while ( iter < max_iter);
        x_next = g(x0);
        err = abs(x_next-x0)
        iter = iter + 1
        if (err < tol)
            root = x_next
            break
        else
            x0 = x_next

        error(iter) = err;
        figure(2)
        num_iter = 1:1:iter;
        loglog(num_iter,error,'g-s')
        grid on
        xlabel('number of iterations')
        ylabel('error')
        end
    end
end
