function [inform,x] = SteepDescent(fun,x,sdparams)
maxit = sdparams.maxit;
toler = sdparams.toler;
params = struct('c1',0.01, 'c2',0.6, 'maxit',100);
alfa = 1;
iter = 0;
inform.status = 0;
x.f = feval(fun,x.p,1);
x.g = feval(fun,x.p,2);
while iter < maxit
    [alfa, x] = StepSize(fun, x, -x.g, 1, params);
    iter = iter + 1;
    if max(abs(x.g)) <= toler * (1+abs(x.f))
        inform.status = 1;
        break;
    end
end
inform.iter = iter;
end

