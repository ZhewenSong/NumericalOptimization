function [inform, x] = BFGS(fun, x, qnparams)
n = size(x.p,1);
H = eye(n); 
maxit = qnparams.maxit;
it = 0;
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);
lsparams = struct('c1',1.0e-4,'c2',0.4,'maxit',20);
inform.status = 0;
global numf numg 
numf = 0;
numg = 0;
while it < maxit
    p = - H * x.g;
    xold = x.p;
    gold = x.g;
    [alpha, x] = StepSize(fun, x, p, 1, lsparams);
    s = x.p - xold;
    y = x.g - gold;
    if it == 0
        H = H * (y'*s)/(y'*y);
    end
    rho = 1/(y'*s);
    H = (eye(n) - rho*s*y') * H * (eye(n) - rho*y*s') + rho*s*s';
    it = it + 1;
    if norm(x.g) <= qnparams.toler * (1 + abs(x.f))
        inform.status = 1;
        break;
    end
end
inform.iter = it;