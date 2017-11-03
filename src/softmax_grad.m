% calculate gradient 
function grad = softmax_grad(x, param)
% grad:  n*1
L = param.L;
n = param.n;
I = eye(n);

B = L - I;
A = inv(diag(x)*B + I);

grad = zeros(n,1);
for i = 1:n
   grad(i) =  B(i,:)*A(:,i);
end

end
