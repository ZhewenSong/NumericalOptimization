function [f, g, G] = qinterp(x, y)
% size(x) = [(n+2)*(n+1)/2, n]
% size(y) = [(n+2)*(n+1)/2, 1]
% example:
% x = [0, 1, 2, 1, 0, 0;
%     0, 0, 0, 1, 2, 1]';
% y = [1, 2.0084, 7.0091, 1.0168, -0.9909, -0.9916]';
n = size(x,2);
q = (n+2)*(n+1)/2;
A = zeros(q);
index = 1;
A(:, index) = ones(q,1); %% corresponding to f
index = index + 1;
for i=1:n  %% corresponding to g
    A(:,index) = x(:,i);
    index = index + 1;
end
for i=1:n %% corresponding to G
    for j=i:n
        if i==j
            A(:, index) = 0.5*x(:,i).*x(:,j);
        else
            A(:, index) = x(:,i).*x(:,j);
        end
        index = index + 1;
    end
end
coeff = A\y;  %% solve for the coefficients
f = coeff(1);
g = coeff(2:n+1);
G = zeros(n);
index = n+2;
for i=1:n
    for j=i:n
        G(i,j) = coeff(index);
        G(j,i) = G(i,j);
        index = index + 1;
    end
end
