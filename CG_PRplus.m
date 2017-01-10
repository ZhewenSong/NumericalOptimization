function [inform,x] = CG_PRplus(fun,x,sdparams)
maxit = sdparams.maxit;
toler = sdparams.toler;
params = struct('c1',0.01, 'c2',0.1, 'maxit',100);
alfa = 1;
it = 0;
inform.status = 0;
x.f = feval(fun,x.p,1);
x.g = feval(fun,x.p,2);
p = -x.g;
while it < maxit
    gold = x.g;
    [alfa, x] = StepSize(fun, x, p, 1, params);
    beta = max(x.g'*(x.g-gold)/(gold'*gold), 0);
    p = -x.g + beta*p;
    it = it + 1;
    if max(abs(x.g)) <= toler * (1+abs(x.f))
        inform.status = 1;
        break;
    end
end
inform.iter = it;

end

