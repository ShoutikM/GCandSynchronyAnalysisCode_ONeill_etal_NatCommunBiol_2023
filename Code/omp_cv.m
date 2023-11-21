function [theta, LL] = omp_cv(n1, X1, n2, X2, s, alpha, GD_iter)
LL = zeros(1,s);
theta_prev = zeros(size(X1,2), 1);

lambda =  exp(X1*theta_prev)./(1+exp(X1*theta_prev));
del = n1' - lambda;
grad = X1'*del;

[~, S] = max(abs(grad));

for r=1:s
   for itr=1:GD_iter
       theta_prev(S) = theta_prev(S)+...
           alpha*(X1(:,S)'*(n1' - exp(X1(:,S)*theta_prev(S))./(1+exp(X1(:,S)*theta_prev(S))) ));
   end
   
   lambda =  exp(X1*theta_prev)./(1+exp(X1*theta_prev));
   del = n1' - lambda;
   grad = X1'*del; grad_tmp = grad;
   grad_tmp(S) = 0;
   [~, j] = max(abs(grad_tmp));
   S = unique([S,j]);
   
	LL(r) = n2*X2*theta_prev - sum(log(1+exp(X2*theta_prev)));

end
theta = theta_prev;

end