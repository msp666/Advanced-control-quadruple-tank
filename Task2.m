%% Task2
NI = evalfr(Gtf, 0) / (evalfr(Gtf(1,1), 0) * evalfr(Gtf(2,2), 0));
NI_eig = eig(NI);
if all(NI_eig > 0)
    disp('The plant could be controlled by a decentralized controller');
end
%% T2.1 simulate the linearized model with decetralized controller
% E = (Gtf - [Gtf(1, 1) 0; 0 Gtf(2, 2)]) * inv([Gtf(1, 1) 0; 0 Gtf(2, 2)])
% figure
% bodemag(E)
s = tf('s');
% Define parameters
if isequal(mode, 1)
    K1 = 5;
    K2 = 5;
    T1 = 0.5;
    T2 = 0.5;
elseif isequal(mode, 2)
end
% build feedback loop
K = [K1*(1 + 1/(T1 * s)), 0;
    0, K2*(1 + 1/(T2 * s))];
Gtf_cl = feedback(G*K, [1,0; 0 1]);
figure;
bode(Gtf_cl);
t = 0:0.1:50;
u = ones(length(t), 2).*[12.1, 12.6];
figure;
h0 = [h10, h20, h30, h40, 0, 0];
lsim(Gtf_cl, u, t, h0);
title('Simulation of decentralized controller of linearised model')

%% T2.2 Model the parametric uncertainty
% Define parametric uncertainty
gp1 = ureal('gp1', g1, 'Range', [0.9, 1.1]*g1);
gp2 = ureal('gp2', g2, 'Range', [0.9, 1.1]*g2);
kp1 = ureal('kp1', k1, 'Range', [0.9, 1.1]*k1);
kp2 = ureal('kp2', k2, 'Range', [0.9, 1.1]*k2);
ap1 = ureal('ap1', a1, 'Range', [0.95, 1.05]*a1);
ap2 = ureal('ap2', a2, 'Range', [0.95, 1.05]*a2);
ap3 = ureal('ap3', a3, 'Range', [0.95, 1.05]*a3);
ap4 = ureal('ap4', a4, 'Range', [0.95, 1.05]*a4);
hp30 = ureal('hp30', h30, 'plusminus', 0.03);
hp40 = ureal('hp40', h40, 'plusminus', 0.04);

Ap = [-ap1/A1*sqrt(g/(2*h10)),          0,               ap3/A3*sqrt(g/(2*h30)),          0;
      0,                       -ap2/A2*sqrt(g/(2*h20)),              0,            ap4/A4*sqrt(g/(2*h40));
      0,                                0,               -ap3/A3*sqrt(g/(2*h30)),         0;
      0,                                0,                           0,            -ap4/A4*sqrt(g/(2*h40))];
Bp = [gp1*kp1/A1,     0;
          0,          gp2*kp2/A2;
          0,          (1-gp2)*kp2/A3;
      (1-gp1)*kp1/A4, 0];
Cp = G.C;
Dp = G.D;
Gp = uss(Ap, Bp, Cp, Dp);
Gp_cl = feedback(Gp*K, [1 0; 0 1]);


% fit the weighting function
w = logspace(-4, 2, 80);
G_mag = zeros(2, 2, length(w), 1000);
for k = 1:1000
    % select the parameter in the uncetain interval via uniform
    % distribution
    p1 = (rand(1,4) * 0.2 + 0.9) .* [g1, g2, k1, k2];
    p2 = (rand(1,4) * 0.1 + 0.95) .* [a1, a2, a3, a4];
    % parateter definition
    A = [-p2(1)/A1*sqrt(g/(2*h10)),          0,               p2(3)/A3*sqrt(g/(2*h30)),          0;
         0,                       -p2(2)/A2*sqrt(g/(2*h20)),              0,            p2(4)/A4*sqrt(g/(2*h40));
         0,                                0,                -p2(3)/A3*sqrt(g/(2*h30)),         0;
         0,                                0,                           0,            -p2(4)/A4*sqrt(g/(2*h40))];
    B = [p1(1)*p1(3)/A1,     0;
          0,          p1(2)*p1(4)/A2;
          0,          (1-p1(2))*p1(4)/A3;
      (1-p1(1))*p1(3)/A4, 0];
    C = G.C;
    D = G.D;
    G_tmp = ss(A, B, C, D);
    
    [Gp_mag(:, : ,:, k)] = squeeze(bode(inv(G) * (G_tmp - G), w)); 
    
end

%%
Gp_mag_max = max(Gp_mag, [], 4);
rd = frd(Gp_mag_max, w);
constraint11.LowerBound = rd(1,1,:); constraint11.UpperBound = [];
constraint12.LowerBound = rd(1,2,:); constraint12.UpperBound = [];
constraint21.LowerBound = rd(2,1,:); constraint21.UpperBound = [];
constraint22.LowerBound = rd(2,2,:); constraint22.UpperBound = [];
WI11 = fitmagfrd(rd(1, 1), 2, [], [], constraint11);
WI12 = fitmagfrd(rd(1, 2), 2, [], [], constraint12);
WI21 = fitmagfrd(rd(2, 1), 2, [], [], constraint21);
WI22 = fitmagfrd(rd(2, 2), 2, [], [], constraint22);
WI = [WI11, WI12; WI21, WI22];
% WI = [WI11, 0; 0, WI22];
figure;
bodemag(1/WI); hold on;
title('Weighting of dynamic uncertainty')

% model the dynamic uncetainty
DeltaI = ultidyn('DeltaI', [2, 2], 'Bound', 1);
Gpdyn_cl = feedback(G*(eye(2) + WI*DeltaI)*K, [1,0;0,1]);

%% T2.3 check the robust stability of dynamic uncertainty
[margin_dyn, wcu_dyn] = robstab(Gpdyn_cl);
if margin_dyn.LowerBound >= 1
    disp('The system is robustly stable');
else
    warning('The system is not robustly stable');
end

%% T2.4 robust stability of paramatric uncertainty

[margin_para, wcu_para] = robstab(Gp_cl);
if margin_para.LowerBound >= 1
    disp('The system is robustly stable');
else
    warning('The system is not robustly stable');
end
% [N, Delta] = lftdata(Gp_cl);
% n = length(Delta.NominalValue);
% Ntf = tf(minreal(N));
% M = Ntf(1:n, 1:n);
% M_hinfnorm = hinfnorm(M);
% if M_hinfnorm < 1
%     disp('This system is robustly stable');
% else
%     warning('This system is not robustly stable');
% end