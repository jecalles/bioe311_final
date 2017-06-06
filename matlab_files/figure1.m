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
m1 = output.m1;
p1 = output.p1;
p2 = output.p2;
p3 = output.p3;

figure(1)
subplot(1,2,1);
hold on;
plot(squeeze(p1(3,3,:)), '--r');
plot(squeeze(p2(3,3,:)), '--g');
plot(squeeze(p3(3,3,:)), '--b');
hold off;
title('Protein Levels vs Time');
xlabel('Time');
ylabel('Conc (arb)');
legend('p1', 'p2', 'p3');

subplot(1,2,2);
plot(squeeze(m1(3,3,:)),squeeze(p1(3,3,:)));
title('P1 vs M1 Trajectory');
xlabel('m1 (arb)');
ylabel('p1 (arb)');

% make a gif!
figure(2)
obj = VideoWriter('Repressilator (No Noise, No Diffusion)'); 
open(obj);

for t = 1:length(p1(1,1,:))
        Z = sin(n*pi.*x/L) * ( A* sin(c*n*pi*t/L) * cos(c*n*pi*t/L) );
        plot(x,Z); ylim([-(A+B)/4,(A+B)/4])
        frame = getframe(gcf);
        writeVideo(obj, frame);
end
close(obj);


