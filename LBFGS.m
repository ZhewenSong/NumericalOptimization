function [inform, x] = LBFGS(fun, x, lbfgsparams)
n = size(x.p,1);
H = eye(n); 
maxit = lbfgsparams.maxit;
k = 0;
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);
lsparams = struct('c1',1.0e-4,'c2',0.4,'maxit',20);
inform.status = 0;
global numf numg 
numf = 0;
numg = 0;
m = lbfgsparams.m;
s = zeros(n, m);
y = zeros(n, m);
rho = zeros(m);
alpha = zeros(m);

while k < maxit
    if k < m
        p = - H * x.g;
    else
        %% Calculating p_k from Algo 7.4
        q = x.g;
        for i = k-1:-1:k-m
            id = mod(i, m) + 1; % because matlab array starts from 1, not 0
            alpha(id) = rho(id) * s(:, id)' * q;
            q = q - alpha(id) * y(:, id);
        end
        id = mod(k-1, m) + 1;
        H = (s(:, id)'*y(:, id))/(y(:, id)'*y(:, id)) * eye(n);
        r = H * q;
        for i = k-m:k-1
            id = mod(i, m) + 1;
            beta = rho(id) * y(:, id)' * r;
            r = r + s(:, id) * (alpha(id) - beta);
        end
        p = -r;
        %%
    end
    xold = x.p;
    gold = x.g;
    [a, x] = StepSize(fun, x, p, 1, lsparams);
    id = mod(k, m) + 1;
    s(:, id) = x.p - xold;
    y(:, id) = x.g - gold;
    rho(id) = 1/(y(:, id)'*s(:, id));    
    if k < m
        if k == 0
            H = H * (y(:, id)'*s(:, id))/(y(:, id)'*y(:, id));
        end
        H = (eye(n) - rho(id)*s(:, id)*y(:, id)') * H * ...
        (eye(n) - rho(id)*y(:, id)*s(:, id)') + rho(id)*s(:, id)*s(:, id)';
    end
    
    k = k + 1;
    if norm(x.g) <= lbfgsparams.toler * (1 + abs(x.f))
        inform.status = 1;
        break;
    end
end
inform.iter = k;