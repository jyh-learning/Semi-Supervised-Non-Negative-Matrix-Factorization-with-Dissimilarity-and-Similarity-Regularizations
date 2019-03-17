function [U,V,obj]=gssnmf(X,D,para,A,W)
% X in d \times N
lambda=para.lambda;  % the hyper papramter lambda

k=para.k;
maxiter=para.maxiter;
mu=para.mu;

d=size(X,1);
n=size(X,2);
L=A-W;

% init
U=rand(d,k);
V=rand(k,n);
obj(1)=sum(sum((X-U*V).^2))+lambda*sum(sum((D.*(V'*V)).^2))+mu*trace(V*L*V');
for iter=1:maxiter
    U=U.*((X*V')./(U*V*V'+eps));
    V=V.*((U'*X+mu*V*W)./(U'*U*V+lambda*V*D+mu*V*A+eps));
    Z=X-U*V;

    % normization
    V=V.*(repmat(sum(U,1)',1,n));
    U=U./(repmat(sum(U,1),d,1));
    obj(iter+1)=sum(sum(Z).^2)+lambda*sum(sum((D.*(V'*V)).^2))+mu*trace(V*L*V');
    disp(['the ', num2str(iter), ' obj is ', num2str(obj(iter))]);
    if (abs(obj(iter+1)-obj(iter))<10^-3)
        break;
    end
end