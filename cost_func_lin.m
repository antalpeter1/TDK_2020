function J = cost_func_lin(u,horizon,q,Q,R,A,B)

niter = horizon;
q1 = zeros(4,niter+1);
q1(:,1) = q; % initial condition

for iter = 1:niter
    q1(:,iter+1) = A*q1(:,iter) + B*u(iter);
end

J = Q(1)*sum(q1(1,:).^2) + Q(3)*sum(q1(3,:).^2) + R*sum(u.^2);