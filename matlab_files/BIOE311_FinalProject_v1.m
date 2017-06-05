tic;
%define system parameters for system containing three diffusing proteins
%and three mRNA.

%initial concentrations
m10 = 1;
m20 = 1;
m30 = 1;

p10 = 1;
p20 = 1;
p30 = 1;

%diffusion constants
Dp1 = 0;
Dp2 = 0;
Dp3 = 0;

%other parameters
alpha = 5*10^-1;
alpha0 = 5*10^-4;
beta = 2*10^2;

%noise parameters
mnoise = 1
pnoise = 1*10^-1

%boundary conditions in time
t = 100;
numStepsT = 1000;
dt = t/numStepsT;

L=5;
x = 1:L;
y = 1:L;
dx = 1;
numStepsX = L/dx;

neumannm1 = (2*Dm1*dt)/(dx)^2;
neumannm2 = (2*Dm2*dt)/(dx)^2;
neumannm3 = (2*Dm3*dt)/(dx)^2;

neumannp1 = (2*Dp1*dt)/(dx)^2;
neumannp2 = (2*Dp2*dt)/(dx)^2;
neumannp3 = (2*Dp3*dt)/(dx)^2;

%define data structure to capture system evolution. each diffusion
%component will be represented by a three-dimensional matrix with
%dimensions x, y, numSteps

m1 = zeros(L,L,numStepsT);
m2 = zeros(L,L,numStepsT);
m3 = zeros(L,L,numStepsT);

p1 = zeros(L,L,numStepsT);
p2 = zeros(L,L,numStepsT);
p3 = zeros(L,L,numStepsT);

%%define coefficient matrix
set = ones(L,1);
cjn = diag(set);
for i=1:numStepsX
    for j=1:numStepsX
        if i == j
            cjn(i,j) = -2;
        end
        if i == j+1
            cjn(i,j) = 1;
        end
        if i == j-1
            cjn(i,j) = 1;
        end
    end
end
%apply no flow BCs
cjn(1,1) = -1;
cjn(end,end) = -1;
%apply sparse function
C = sparse(cjn);

%establish noise random population
m1noise = randomize2(m1, mnoise, L);
m2noise = randomize2(m2, mnoise, L);
m3noise = randomize2(m3, mnoise, L);
p1noise = randomize2(p1, pnoise, L);
p2noise = randomize2(p2, pnoise, L);
p3noise = randomize2(p3, pnoise, L);

%% evaluate PDEsusing diffp and diffm functions
for t = 1:numStepsT
m1 = diffm(m1, alpha, alpha0, p3, dx, t, dt, mnoise, C);
m2 = diffm(m2, alpha, alpha0, p2, dx, t, dt, mnoise, C);
m3 = diffm(m3, alpha, alpha0, p1, dx, t, dt, mnoise, C);

p1 = diffp(p1, Dp1, beta, m1, dx, t, dt, pnoise, C);
p2 = diffp(p2, Dp2, beta, m2, dx, t, dt, pnoise, C);
p3 = diffp(p3, Dp3, beta, m3, dx, t, dt, pnoise, C);
end

%%
figure(1)
T = 1:1:t+1;
plot(T,squeeze(p1(3,3,:)))

%%
figure
subplot(2,2,1);
imagesc(1:L,1:L,p1(:,:,1));
colorbar;
title('T=1');

subplot(2,2,2);
imagesc(1:L,1:L,p1(:,:,250));
colorbar;
title('T=10');

subplot(2,2,3);
imagesc(1:L,1:L,p1(:,:,500));
colorbar;
title('T=30');

subplot(2,2,4);
imagesc(1:L,1:L,p1(:,:,numStepsT));
colorbar;
title('T=100');

time = toc