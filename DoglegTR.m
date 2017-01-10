function [inform,x] = DoglegTR(fun, x, trparams)
    
    maxit = trparams.maxit;
    delta = trparams.delta;
    hatDelta = trparams.hatDelta;
    eta = trparams.eta;
    Delta0 = trparams.Delta0;
    toler = trparams.toler;
    x.f = feval(fun,x.p,1);
    x.g = feval(fun,x.p,2);
    x.h = feval(fun,x.p,4);
    Delta = Delta0;
    it = 0;
    %% main loop
    while norm(x.g) > toler && it < maxit
        %% modify the Hessian if need be
        [Q,Lambda] = eig(x.h);
        updated = 0;
        for i=1:size(Lambda)
            if Lambda(i,i) < delta
               Lambda(i,i) = delta;
               updated = 1;
            end
        end
        if updated
            B = Q*Lambda*Q';
        else
            B = x.h;
        end
        %% solve for p
        pB = -inv(B) * x.g;
        pU = -(x.g'*x.g)/(x.g'*B*x.g) * x.g;
        if norm(pB) <= Delta
            p = pB;
        else        
            a = (pB-pU)'*(pB-pU);
            b = 2*pU'*(pB-pU);
            c = pU'*pU - Delta^2;
            t1 = (-b + sqrt(b^2-4*a*c))/2/a;
            t2 = (-b - sqrt(b^2-4*a*c))/2/a;
            if t1 < -1 || t1 > 1
                tau = t2 + 1;
            else
                tau = t1 + 1;
            end
            if tau >= 0 && tau <= 1
                p = tau * pU;
            end
            if tau > 1 && tau <= 2
                p = pU + (tau - 1)*(pB - pU);
            end
        end
        
        %%
        m = @(p) x.f + x.g'*p + p'*B*p/2; 
        rho = (x.f - feval(fun,x.p+p,1))/(x.f - m(p));
        
        if rho < 1/4
            Delta = Delta/4;
        else
            if rho > 3/4 && norm(p) == Delta
                Delta = min(2*Delta, hatDelta);
            end
        end
        if rho > eta
            x.p = x.p + p;
            x.f = feval(fun,x.p,1);
            x.g = feval(fun,x.p,2);
            x.h = feval(fun,x.p,4);
        end
        
        
        it = it + 1;
        inform.iter = it;
        fprintf(1,'iter %3d: f=%12.5e, ||Df||=%12.5e, Delta=%7.2e\n',...
            inform.iter, x.f, norm(x.g), Delta);

    end
    %%
    
    if norm(x.g) <= toler
        inform.status = 1;
    else
        inform.status = 0;
    end
