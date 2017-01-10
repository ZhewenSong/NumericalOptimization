mu=0.01; L=1; kappa=L/mu;
n=100;
A = randn(n,n); [Q,R]=qr(A);
D=rand(n,1); D=10.^D; Dmin=min(D); Dmax=max(D);
D=(D-Dmin)/(Dmax-Dmin);
D = mu + D*(L-mu);
A = Q'*diag(D)*Q;

x_opt = zeros(n,1);
f = @(x) x'*A*x/2;

%%%% Steepest descent with alpha_k = 1/L. %%%%%

av_sd = 0;
for i=1:10
    x0 = randn(n,1);
    epsilon = f(x0) - f(x_opt);
    x = x0;
    while epsilon > 1e-6 
        alpha = 1/L;
        x = x - alpha * A*x;
        epsilon = f(x) - f(x_opt);
        av_sd = av_sd + 1;
    end
end
av_sd = av_sd/10;

%%%% Steepest descent with exact line search. %%%%%

av_sde = 0;
for i=1:10
    x0 = randn(n,1);
    epsilon = f(x0) - f(x_opt);
    x = x0;
    while epsilon > 1e-6 
        alpha = (x'*A^3*x)\(x'*A^2*x);
        x = x - alpha * A*x;
        epsilon = f(x) - f(x_opt);
        av_sde = av_sde + 1;
    end
end
av_sde = av_sde/10;

%%%% Nesterov's optimal method %%%%%

av_nest = 0;
for i=1:10
    x0 = randn(n,1);
    epsilon = f(x0) - f(x_opt);
    x = x0;
    x_old = x0;
    while epsilon > 1e-6 
        alpha = 1/L;
        beta = (sqrt(L) - sqrt(mu))/(sqrt(L) + sqrt(mu)); 
        y = x + beta * (x - x_old);
        x_old = x;
        x = y - alpha * A*y;
        epsilon = f(x) - f(x_opt);
        av_nest = av_nest + 1;
    end
end
av_nest = av_nest/10;

%%%% CG method %%%%%

av_cg = 0;
for i=1:10
    x0 = randn(n,1);
    epsilon = f(x0) - f(x_opt);
    x = x0;
    r = A*x0;
    p = -r;
    while epsilon > 1e-6 
        alpha = -(p'*A*p)\(r'*p);
        x = x + alpha*p;
        r = A*x;
        beta = (p'*A*p)\(r'*A*p);
        p = -r + beta*p;
        epsilon = f(x) - f(x_opt);
        av_cg = av_cg + 1;
    end
end
av_cg = av_cg/10;

fprintf(1,' steepest descent - fixed steps : %7.1f\n', av_sd);
fprintf(1,' steepest descent - exact steps : %7.1f\n', av_sde);
fprintf(1,' Nesterov                       : %7.1f\n', av_nest);
fprintf(1,' conjugate gradient             : %7.1f\n', av_cg);
