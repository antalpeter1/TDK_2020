function J = cost_func_simple(u,horizon,q,Q,R,A,B,xtrain,L,alpha,Kinv)

niter = horizon;
mean = zeros(4,niter+1);
mean(:,1) = q; % initial condition
sf = 0.1;

Ktest = zeros(1,length(xtrain(1,:)));
M = zeros(4,1);
S = M;
N = length(xtrain(1,:));

for iter = 1:niter
    for k = 1:4
        for j=1:N
            Ktest(1,j,k)=sf^2*kernel(mean(:,iter),xtrain(:,j),L(:,:,k));
        end
        M(k) = Ktest(:,:,k)*alpha(:,:,k);
        S(k) = sf^2 - Ktest(:,:,k)*Kinv(:,:,k)*Ktest(:,:,k)';
    end
    mean(:,iter+1) = A*mean(:,iter) + B*u(iter) + M;
end
J = Q(1)*(sum(mean(1,:).^2) + S(1)) +...
    Q(3)*(sum(mean(3,:).^2) + S(3)) + R*sum(u.^2);