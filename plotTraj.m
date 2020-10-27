function plotTraj(datas,Ts)
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',12);

tsimu = 0:Ts:(length(datas(1,:))-1)*Ts;
figure()
subplot(3,1,1)
plot(tsimu,datas(1,:),'k','LineWidth',1);
xlabel('$t$ (s)')
ylabel('$\varphi$ (rad)')
box off
subplot(3,1,2)
plot(tsimu,datas(3,:),'k','LineWidth',1);
xlabel('$t$ (s)')
ylabel('$x$ (m)')
box off
subplot(3,1,3)
plot(tsimu,datas(5,:),'k','LineWidth',1);
xlabel('$t$ (s)')
ylabel('$u$ (V)')
box off