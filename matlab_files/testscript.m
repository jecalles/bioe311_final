%% Declare and package parameters

%length and width
L=3;
W=4;

%initial conditions
initCond = zeros(L,W,6);
initCond(:,:,2) = ones(L,W);
initCond(:,:,4) = ones(L,W)*2;
initCond(:,:,6) = ones(L,W)* 5;

%diffusion constants
Dp1 = 0;
Dp2 = 0;
Dp3 = 0;
D_pi = [Dp1 Dp2 Dp3];

%diff eq parameters
alpha = 1000;
alpha0 = 1;
beta = 5;
param = [alpha alpha0 beta];

%noise parameters
mnoise = 0;
pnoise = 0; %mnoise/10;
noiseParam = [mnoise pnoise];

%time simulation parameters
t = 100;
dt = 0.01;

%periodicity
periodic_x = 0;
periodic_y = 0;
periodicity = [periodic_x periodic_y];

%% Simulate system

% instantiate object
sim = repressilator(L, W, initCond, D_pi, param, noiseParam, periodicity);

% simulate
output = sim.simulate(t, dt);

%% plot output
p2 = output.p2;
p1 = output.p1;

figure(1)
hold on;
plot(squeeze(p1(3,3,:)), '--r');
plot(squeeze(p2(3,3,:)), '--b');
hold off;
title('TetR species vs Time');
legend('p1', 'm1');

figure(2)
plot(squeeze(m1(3,3,:)),squeeze(p1(3,3,:))

%

