%% Non minimum phase
%% Question 2
run('quadtank_nonminphase.m')
percentage1 = 10; % uncertainty in % 
percentage2 = 5; % uncertainty in % 
% define uncertain parameters
gp1 = ureal('g1', g1, 'Percentage', percentage1);
gp2 = ureal('g2', g2, 'Percentage', percentage1);
kp1 = ureal('k1', k1, 'Percentage', percentage1);
kp2 = ureal('k2', k2, 'Percentage', percentage1);

ap1 = ureal('a1', a1, 'Percentage', percentage2);
ap2 = ureal('a2', a2, 'Percentage', percentage2);
ap3 = ureal('a3', a3, 'Percentage', percentage2);
ap4 = ureal('a4', a4, 'Percentage', percentage2);

% set up uncertain TF based on symbolic calculation (Task1_Q2)
s = tf('s');
% from state space:
A = [ -ap1/A1*sqrt( (g/(2*h10)) ), 0, ap3/A1*sqrt( (g/(2*h30)) ), 0;...
         0, -ap2/A2*sqrt( (g/(2*h20)) ), 0, ap4/A2*sqrt( (g/(2*h40)) );...
         0, 0, -ap3/A3*sqrt( (g/(2*h30)) ), 0;...
         0, 0, 0, -ap4/A4*sqrt( (g/(2*h40)) ) ];
B = [ (gp1*kp1)/A1, 0;...   
         0, (gp2*kp2)/A2;...
         0, ((1-gp2)*kp2)/A3;...
         ((1-gp1)*kp1)/A4, 0 ];
C = [1, 0, 0, 0; 0, 1, 0, 0];
D = zeros(size(C,1),size(B,2));
Gp_np = uss(A,B,C,D);
% from symbolic transfer function:
% Gp_np = [gp1*kp1/(A1*s+ap1*sqrt(g/(2*h10))), -ap3*kp2*(gp2-1)*sqrt(g/(2*h30))/((A1*s+ap1*sqrt(g/(2*h10)))*(A3*s+ap3*sqrt(g/(2*h30))));...
%          -ap4*kp1*(gp1-1)*sqrt(g/(2*h40))/((A2*s+ap2*sqrt(g/(2*h20)))*(A4*s+ap4*sqrt(g/(2*h40)))), gp2*kp2/(A2*s+ap2*sqrt(g/(2*h30)))];
% G_np = Gp_np.NominalValue;
% G_np = zpk(Gp_np.nominal); % normal TF (without uncertainty)
G_np = tf(sys_np);
% bodemag ( usample ( Gp_np, 100 ) ) ;
%% Fit the weighting function
% construct closed loop with controller K_np (Task 2 Q1)
K3 = 0.7030; K4 = 0.6516;
T3 = 90.9288; T4 = 181.19; 
K_np = [0, K4*(1+1/(T4*s)); K3*(1+1/(T3*s)), 0];
% according to Skogestad p.294
w = logspace(-8, 4, 1000);
G_mag = zeros(2, 2, length(w), 1000);
SV = zeros(1000, length(w));
tic;
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
    C = [1, 0, 0, 0; 0, 1, 0, 0];
    D = zeros(size(C,1),size(B,2));
    G_temp = ss(A, B, C, D);
    % [Gp_mag(:, : ,:, k)] = squeeze(bode((G_temp - G_np)*inv(G_np), w));
    % for i = 1:length(w)
    %      Gf = freqresp((G_temp-G_np)/G_np, w(i));
    %      SV(k,i) = max(svd(Gf));
    % end   
    SV(k,:) = max( sigma((G_temp-G_np)/G_np, w), [], 1 );
    
end
toc
SV_max = max(SV, [], 1);
rd = frd (SV_max, w);
constraint.LowerBound = rd;
constraint.UpperBound = [];
WO1 = fitmagfrd(rd, 2, [], [], constraint);
WO2 = fitmagfrd(rd, 4, [], [], constraint);
figure;
bodemag(rd, WO1,'r', WO2, 'g', {0.0001,100})
title('Singular value weight')
legend('lo(\omega)','|WO|(j\omega) (order 2)','|WO|(j\omega) (order 4)')
WO = WO2;
% Gp_mag_max = max(Gp_mag, [], 4);
% rd = frd(Gp_mag_max, w);
% constraint11.LowerBound = rd(1,1,:); constraint11.UpperBound = [];
% constraint12.LowerBound = rd(1,2,:); constraint12.UpperBound = [];
% constraint21.LowerBound = rd(2,1,:); constraint21.UpperBound = [];
% constraint22.LowerBound = rd(2,2,:); constraint22.UpperBound = [];
% WO11 = fitmagfrd(rd(1, 1), 2, [], [], constraint11);
% WO12 = fitmagfrd(rd(1, 2), 2, [], [], constraint12);
% WO21 = fitmagfrd(rd(2, 1), 2, [], [], constraint21);
% WO22 = fitmagfrd(rd(2, 2), 2, [], [], constraint22);
% WO = [WO11, WO12; WO21, WO22];
% figure;
% bodemag(1/WO); hold on;
% title('Weighting of dynamic uncertainty')

% model the dynamic uncetainty
DeltaO = ultidyn('DeltaO', [2, 2], 'Bound', 1);
Tpdyn_np = feedback((eye(2) + WO*DeltaO)*G_np*K_np, eye(2));

figure;
bodemag(usample((Gp_np-G_np)/G_np, 100)); hold on;
bodemag(WO*DeltaO,'r'); 
title('Comparison: parametric vs. dynamic uncertainty')
legend('parametric', 'dynamic')

%% Question 3: check the robust stability of dynamic uncertainty
[margin_dyn, wcu_dyn] = robstab(Tpdyn_np);
if margin_dyn.LowerBound >= 1
    disp(['The system is robustly stable (dynamic uncertainty) - stability margin: ', num2str(margin_dyn.LowerBound)]);
else
    warning(['The system is not robustly stable (dynamic uncertainty) - stability margin: ', num2str(margin_dyn.LowerBound)]);
end
%% Question 4: robust stability of paramatric uncertainty 
Tppara_np = feedback(Gp_np*K_np, eye(2));
[margin_para, wcu_para] = robstab(Tppara_np);
if margin_para.LowerBound >= 1
    disp(['The system is robustly stable (parametric uncertainty) - stability margin: ', num2str(margin_para.LowerBound)]);
else
    warning(['The system is not robustly stable (parametric uncertainty) - stability margin: ', num2str(margin_para.LowerBound)]);
end
w = logspace(-8,1,200); 
[~, ~, info_para] = robstab(Tppara_np, w);
[~, ~, info_dyn] = robstab(Tpdyn_np,w);
figure;
semilogx(w, info_para.Bounds(:,1))
hold on;
semilogx(w, info_dyn.Bounds(:,1))
title('Stability Margin vs. Frequency')
ylabel('Margin (lower bound)')
xlabel('Frequency')
legend('parametric','dynamic')