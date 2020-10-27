%% LQR for the cart-pole system

load('cp_params.mat');
syms phi phid x xd I u
fr(1) = 0.2;              % friction magnitude
Ts = 0.01;
Duration = 2.5;
Q = [100 0 1 0];          % LQ weights
R = 0.0001;

q1 = [phi; phid; x; xd];
q1 = cp_dyn_v3(q1,u,par,fr);

Asym = jacobian(q1,[phi;phid;x;xd]);
A = double(subs(Asym,[phi;phid;x;xd;u],[0;0;0;0;0]));
Bsym = jacobian(q1,u);
B = double(subs(Bsym,[phi;phid;x;xd;u],[0;0;0;0;0]));
A = eye(4)+Ts*A;
B = Ts*B;
clear phi phid x xd I u k1 k2 k3 k4 q1

%% LQR
C = diag([1,1,1,1]);
D = zeros(4,1);
K = dlqr(A,B,diag(Q),R);

%% Simulation
qLQR = zeros(5,Duration/Ts+1);
q = [0.2;0;0;0];
qLQR(:,1) = [q; 0];
for i = 1 : Duration/Ts
    u = -K*q;
    [~,q] = ode45(@(t,q) cp_dyn_v3(q,u,par,fr),[0 Ts/2 Ts],q);
    q = q(3,:)';
    qLQR(:,i+1) = [q; u];
end
%% Plots
plotTraj(qLQR,Ts)
