%% Task3
s = tf('s');

controller_select = 0;
noise_setting = 0;
%% T3.1 weighting function
% dcgain = 1e5;
% hfgain = 0.1;
% omB1 = 0.5;
% omB2 = 0.5;
% Wp11 = makeweight(dcgain, [omB1, 2], hfgain);
% Wp22 = makeweight(dcgain, [omB1, 2], hfgain);
% Wp = blkdiag(Wp11, Wp22);

% Wp = 100*(s+0.05) * (3*s+1)/((s+0.00001) * (s+0.1));
k=2;
ws = 10^(k-1)*0.014;
Wp = (0.5*s+0.5)/(s+0.5*1e-2);
Wp = blkdiag(Wp, Wp);

% dcgain_u = 0.2;
% hfgain_u = 1e8;
% omBu = 0.01;
% Wu = makeweight(dcgain_u, [omBu, 0.8], hfgain_u);
Wu = (s+0.0001)/(s+0.1);
Wu = blkdiag(Wu, Wu);

% dcgain_T = 0;
% hfgain_T = 1e3;
% omBT = 0.1;
% WT = makeweight(dcgain_T, [omBT, 0.2], hfgain_T);
WT = (1e4*s+6)/(s+3);
WT = blkdiag(WT,WT);
% bodemag(WT);
% legend('Wp', 'WT')

% [sv_Wo, w_Wo] = sigma(WI);
% rd_Wo = frd(1.1*sv_Wo(1,:), w_Wo);
% WT = fitmagfrd(rd_Wo, 2);
% WT = blkdiag(WT, WT);
figure(9)
sigma(Wp, WT, WI);


%% T3.2 H-inf norm controller design
[K_hinf, hinfnorm_cl, gamma] = mixsyn(G, Wp, [], WT);
K0 = blkdiag(1.5, 1.5);

S_hinf = feedback(eye(2), G*K_hinf);
T_hinf = minreal(eye(2) - S_hinf);
KS_hinf = K_hinf * S_hinf;

figure(10)
step(K0*T_hinf);
title('Unit step response of H-inf norm controller');
info_hinf = stepinfo(T_hinf);

t = 0:0.1:2000;
u = ones(length(t), 2).*[12.1, 12.6];
figure(11);
h0 = zeros(1,length(T_hinf.A));
% h0 = [h10, h20, h30, h40, 0, 0];
lsim(T_hinf, u, t, h0);
title('Simulation of H-inf controller')

figure(12)
sigma(G*K_hinf, G, Wp);
legend('GK', 'G', 'Wp')
title('Open-loop function'); hold on
figure(13)
sigma(S_hinf, gamma/Wp);
legend('S', 'gamma/Wp');
title('Sensitive function');hold on
figure(14)
sigma(T_hinf, gamma/WT);
legend('T', 'gamma/WT');hold on

% Gpdyn_hinf_cl = feedback(G*(eye(2) + WI*DeltaI)*K_hinf, [1,0;0,1]);
Gpdyn_hinf_cl = eye(2) - feedback(eye(2), Gp*K_hinf);

[margin_hinf, wcu_hinf] = robstab(Gpdyn_hinf_cl);


%% T3.3 Non-linear model
noise_setting = 1;
%% T3.4 Robust performance
[perfmarg_hinf, wcu_perf_hinf] = robgain(Gpdyn_hinf_cl, 1.5);
[perfmarg_PI, wcu_perf_PI] = robgain(Gp_cl, 2);

