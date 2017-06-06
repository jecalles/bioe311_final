%% Declare and package parameters

%length and width
L=3;
W=4;

%initial conditions
initCond = zeros(L,W, 6);
initCond(:,:,4) = ones(L,W);

%diffusion constants
Dp1 = 0;
Dp2 = 0;
Dp3 = 0;
D_pi = [Dp1 Dp2 Dp3];

%diff eq parameters
alpha = 5*10^-1;
alpha0 = 5*10^-4;
beta = 300;
param = [alpha alpha0 beta];

%noise parameters
mnoise = 0;
pnoise = 0; %mnoise/10;
noiseParam = [mnoise pnoise];

%time simulation parameters
t = 100;
dt = 0.1;

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
p1 = output.p1;

figure(1)
plot(squeeze(p1(3,3,:)))

%%
figure(2)
subplot(2,2,1);
imagesc(p1(:,:,1));
colorbar;
title('T=1');

subplot(2,2,2);
imagesc(p1(:,:,30));
colorbar;
title('T=10');

subplot(2,2,3);
imagesc(p1(:,:,500));
colorbar;
title('T=30');

subplot(2,2,4);
imagesc(p1(:,:,end));
colorbar;
title('T=100');

