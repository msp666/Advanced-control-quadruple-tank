%% Task3
s = tf('s')
%% T3.1 weighting function
close all
% dcgain = 1e5;
% hfgain = 0.1;
% omB1 = 0.5;
% omB2 = 0.5;
% Wp11 = makeweight(dcgain, [omB1, 2], hfgain);
% Wp22 = makeweight(dcgain, [omB1, 2], hfgain);
% Wp = blkdiag(Wp11, Wp22);

% Wp = 100*(s+0.05) * (3*s+1)/((s+0.00001) * (s+0.1));
Wp = (15*s+0.01)/((s+0.0001));
figure
bodemag(Wp);hold on
Wp = blkdiag(Wp, Wp);

dcgain_u = 0.1;
hfgain_u = 1e8;
omBu = 0.01;
Wu = makeweight(dcgain_u, [omBu, 0.8], hfgain_u);
Wu = blkdiag(Wu, Wu);

% dcgain_T = 0;
% hfgain_T = 1e3;
% omBT = 0.1;
% WT = makeweight(dcgain_T, [omBT, 0.2], hfgain_T);
% WT = 1000*s/(s+2) * (s+0.3)/(s+10);
WT = 10*s/(s+0.1);
bodemag(WT);
legend('Wp', 'WT')
WT = blkdiag(WT, WT);


%% T3.2 H-inf norm controller design
[K_hinf, hinfnorm_cl, gamma] = mixsyn(G, Wp, [], WT);

S_hinf = feedback(eye(2), G*K_hinf);
T_hinf = eye(2) - S_hinf;
KS_hinf = K_hinf * S_hinf;

figure
step(T_hinf);
title('Unit step response of H-inf norm controller');
t = 0:0.1:700;
u = ones(length(t), 2).*[12.1, 12.6];
figure;
h0 = zeros(1,length(T_hinf.A));
lsim(T_hinf, u, t, h0);
title('Simulation of H-inf controller')

figure
sigma(G*K_hinf, G, Wp);
legend('GK', 'G', 'Wp')
title('Open-loop function')
figure
sigma(S_hinf, gamma/Wp);
legend('S', 'gamma/Wp');
title('Sensitive function')
figure
sigma(T_hinf, gamma/WT);
legend('T', 'gamma/Wy')

% Gpdyn_hinf_cl = feedback(G*(eye(2) + WI*DeltaI)*K_hinf, [1,0;0,1]);
Gpdyn_hinf_cl = feedback(Gp*K_hinf, [1,0;0,1]);
% robstabopt = robOptions('VaryFrequency', 'on');
% [margin_hinf, mcu_hinf, info_hinf] = robstab(Gpdyn_hinf_cl, logspace(-10, 0, 100), robstabopt);
[margin_hinf, mcu_hinf] = robstab(Gpdyn_hinf_cl);
% figure
% semilogx(info_hinf.Frequency, info_hinf.Bounds);
% title('Stability margin over frequency')
% legend('Lower bound', 'Upper bound')