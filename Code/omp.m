function [theta, LL] = omp(n, X, s, alpha, GD_iter)
theta_prev = zeros(size(X,2), 1);

lambda =  exp(X*theta_prev)./(1+exp(X*theta_prev));
del = n' - lambda;
grad = X'*del;

[~, S] = max(abs(grad));

for r=1:s
   for itr=1:GD_iter
       theta_prev(S) = theta_prev(S)+...
           alpha*(X(:,S)'*(n' - exp(X(:,S)*theta_prev(S))./(1+exp(X(:,S)*theta_prev(S))) ));
   end
   
   lambda =  exp(X*theta_prev)./(1+exp(X*theta_prev));
   del = n' - lambda;
   grad = X'*del; grad_tmp = grad;
   grad_tmp(S) = 0;
   [~, j] = max(abs(grad_tmp));
   S = unique([S,j]);

end
theta = theta_prev;

LL = n*X*theta - sum(log(1+exp(X*theta)));

end