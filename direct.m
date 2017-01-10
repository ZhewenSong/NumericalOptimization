function [inform, x] = direct(fun, x, directparams)
maxit = directparams.maxit;
tol = directparams.toler;
it = 0;
n = size(x.p, 1);
D = [eye(n), -eye(n)];
gamma = 1;  
global numf
numf = 0;

inform.status = 0;
x.f = feval(fun, x.p, 1);
phi = directparams.phi;
theta = directparams.theta;
while it < maxit
    if gamma <= tol
        inform.status = 1;
        break
    end
    found = 0;
    index = randperm(2*n);
    for p=index
        xnew = x.p + gamma * D(:,p);
        fnew = feval(fun, xnew, 1);
        if fnew < x.f - gamma^2
            found = 1;
            x.p = xnew;
            x.f = fnew;
            gamma = phi*gamma;
            break;
        end
    end
    if found == 0
        gamma = theta*gamma;
    end
    it = it + 1;
end
inform.iter = it;