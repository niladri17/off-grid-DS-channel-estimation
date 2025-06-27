function [h_sparse, h_est, support] = omp(A,y,k)

r = y;
S = [];
[m,n]=size(A);
Phi = [];
colnorms = sqrt(diag(A'*A))';
for t = 1 : k
  lambda = (abs(r'*A)./colnorms).^2;
  [vals ind] = sort(lambda,'descend');
  S=[S ind(1)];
  Phi = [Phi A(:,ind(1))];
  xt = pinv(Phi)*y;
  at = Phi*xt;
  r = y - at;
end

sol = zeros(n,1);
sol(S) = xt;
h_sparse = sol;
h_est = xt;
support = S;
