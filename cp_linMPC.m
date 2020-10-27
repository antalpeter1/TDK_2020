%% MPC for the cart-pole system

load('cp_params.mat');
syms phi phid x xd I u
fr(1) = 0.2;              % friction magnitude
Ts = 0.01;
horizon = 25;
Duration = 2.5;
Q = [100 0 1 0];          % MPC weights
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

%% LQ feedback gain
C = diag([1,1,1,1]);
D = zeros(4,1);
K = dlqr(A,B,diag(Q),R);

%% Simulation
A2 = A-B*K;
u = 0;
qlinMPC = zeros(5,Duration/Ts+1);
q = [0.2; 0; 0; 0];        % initial condition
qlinMPC(:,1) = [q; u];
input = zeros(horizon,1);
options = optimoptions('fmincon','Algorithm','sqp');

for i = 1: Duration / Ts
    input = fmincon(@(u) cost_func_lin(u,horizon,q,Q,R,A2,B),input,[],...
        [],[],[],(K*q-40)*ones(horizon,1),(K*q+40)*ones(horizon,1),[],options);
    u = -K*q + input(1);
    [~,q] = ode45(@(t,q) cp_dyn_v3(q,u,par,fr),[0 Ts/2 Ts],q);
    q = q(3,:)';
    qlinMPC(:,i+1) = [q; u];
end
%% Plots
plotTraj(qlinMPC,Ts)

%% Training points
xtrain = zeros(4,15);
ytrain = xtrain;
for i = 1:length(ytrain(1,:))
    xtrain(:,i) = qlinMPC(1:4,2*i);
    ytrain(:,i) = qlinMPC(1:4,2*i+1)-A*qlinMPC(1:4,2*i)-B*qlinMPC(5,2*i);
end