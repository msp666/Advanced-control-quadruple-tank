%% Task 2: Decentralised control and uncertainty modeling(Non-minimum Phase System)
% Calculate Niederlinski Index(NI)
% For non-minimum phase system
G_tf_n = tf(minreal(G_n));
G_0_n = evalfr(G_tf_n,0);
G_0_n_det = det(G_0_n);
NI_n1 = G_0_n_det/(G_0_n(1,1)*G_0_n(2,2)); 
NI_n2 = -G_0_n_det/(G_0_n(1,2)*G_0_n(2,1));

%% Task 2.1 Question 1: Decentralised controller based on linear and nonlinear models
% Non-minimum phase system

s = tf('s');
Tao1 = A1/a1*sqrt(2*h10/g);
Tao2 = A2/a2*sqrt(2*h20/g);
Tao3 = A3/a3*sqrt(2*h30/g);
Tao4 = A4/a4*sqrt(2*h40/g);
c1 = Tao1*k1/A1;
c2 = Tao2*k2/A2;

% Boundary values for K and T
K1_B_n = (9.2*Tao1 - 45)/(45*g1*c1);
T1_B_n = 4*Tao1*K1_B_n*g1*c1*(log(10))^2/((K1_B_n*g1*c1+1)^2*(pi^2 + (log(10))^2));
K2_B_n = (9.2*Tao2 - 45)/(45*g2*c2);
T2_B_n = 4*Tao2*K2_B_n*g2*c2*(log(10))^2/((K2_B_n*g2*c2+1)^2*(pi^2 + (log(10))^2));

% Choose K1,T1,K2,T2
K1_n = 0.5;
K2_n = 0.6;
T1_n = 100;
T2_n = 110;

% Calculation of closed-loop transfer function, plot and display
Gtf_n = [g1*c1/(1+Tao1*s)                  (1-g2)*c1/((1+Tao3*s)*(1+Tao1*s));
      (1-g1)*c2/((1+Tao4*s)*(1+Tao2*s))     g2*c2/(1+Tao2*s)               ];
K_n = [0, K1_n*(1+1/(T1_n*s)); K2_n*(1+1/(T2_n*s)), 0];
Gtf_cl_n = feedback(G_n*K_n,[1 0; 0 1]);
figure;
step(Gtf_cl_n);
grid on;
stepinfo(Gtf_cl_n);
%% Task 2.2 Question 2: Model uncertainty and weighting function
% Non-minimum phaase
% Uncertain parameters of the model
g_u1 = ureal('g_u1',g1,'Percentage',[-10,10]);
g_u2 = ureal('g_u2',g2,'Percentage',[-10,10]);
k_u1 = ureal('k_u1',k1,'Percentage',[-10,10]);
k_u2 = ureal('k_u2',k2,'Percentage',[-10,10]);
a_u1 = ureal('a_u1',a1,'Percentage',[-5,5]);
a_u2 = ureal('a_u2',a2,'Percentage',[-5,5]);
a_u3 = ureal('a_u3',a3,'Percentage',[-5,5]);
a_u4 = ureal('a_u4',a4,'Percentage',[-5,5]);


A_u = [-a_u1/A1*sqrt(g/(2*h10)),                          0,   a_u3/A3*sqrt(g/(2*h30)),                         0;
                              0,   -a_u2/A2*sqrt(g/(2*h20)),                         0,   a_u4/A4*sqrt(g/(2*h40));
                              0,                          0,  -a_u3/A3*sqrt(g/(2*h30)),                         0;
                              0,                          0,                         0,   -a_u4/A4*sqrt(g/(2*h40))];
B_u = [   g_u1*k_u1/A1,                  0;
                     0,       g_u2*k_u2/A2;
                     0,   (1-g_u2)*k_u2/A3;
      (1-g_u1)*k_u1/A4,                  0];
C_u = G_n.C;
D_u = G_n.D;
G_u_n = ss(A_u, B_u, C_u, D_u);
G_u_cl_n = feedback(G_u_n*K_n, [1 0; 0 1]);

%%
w = logspace(-3, 3, 100);
G_mag_n = zeros(2, 2, length(w), 1500);
for k = 1:1500
    % select the parameter in the uncetain interval via uniform
    % distribution
    p1 = (rand(1,4) * 0.2 + 0.9) .* [g1, g2, k1, k2];
    p2 = (rand(1,4) * 0.1 + 0.95) .* [a1, a2, a3, a4];
    % parateter definition
    A = [-p2(1)/A1*sqrt(g/(2*h10)),                             0,      p2(3)/A3*sqrt(g/(2*h30)),                              0;
                                 0,     -p2(2)/A2*sqrt(g/(2*h20)),                              0,      p2(4)/A4*sqrt(g/(2*h40));
                                 0,                             0,      -p2(3)/A3*sqrt(g/(2*h30)),                             0;
                                 0,                             0,                              0,      -p2(4)/A4*sqrt(g/(2*h40))];
    B = [   p1(1)*p1(3)/A1,                    0;
                         0,       p1(2)*p1(4)/A2;
                         0,   (1-p1(2))*p1(4)/A3;
         (1-p1(1))*p1(3)/A4,                   0];
    C = G_n.C;
    D = G_n.D;
    G_tmp_n = ss(A, B, C, D);
    
    [Gp_mag_n(:, : ,:, k)] = squeeze(bode(inv(G_n) * (G_tmp_n - G_n), w)); 
end


% Get Maximum Magnitude of each frequency 
Gp_mag_max_n = max(Gp_mag_n, [], 4);
FR = frd(Gp_mag_max_n, w);

% Fit weighting function,plot and display
C11.LowerBound = FR(1,1,:); C11.UpperBound = [];
C12.LowerBound = FR(1,2,:); C12.UpperBound = [];
C21.LowerBound = FR(2,1,:); C21.UpperBound = [];
C22.LowerBound = FR(2,2,:); C22.UpperBound = [];
WM11_n = fitmagfrd(FR(1, 1), 2, [], [], C11);
WM12_n = fitmagfrd(FR(1, 2), 2, [], [], C12);
WM21_n = fitmagfrd(FR(2, 1), 2, [], [], C21);
WM22_n = fitmagfrd(FR(2, 2), 2, [], [], C22);
WM_n = [WM11_n, WM12_n; WM21_n, WM22_n];
figure;
bodemag(1/WM_n ,  '-*g' ); grid on; hold on;
title('Weighting of dynamic uncertainty for non-minimum phase system');

%% Task 2.3 Question: Robust stability with dynamic uncertainty
% Non-minimum Phase

DeltaM_n = ultidyn('DeltaM_n', [2, 2], 'Bound', 1);
% Additive uncertainty
%Gp_n = G_n + WM_n * DeltaM_n;
Gp_n = G_n*(eye(2) + WM_n*DeltaM_n);
L_n = Gp_n * K_n;
T_n = feedback(L_n ,[1 0;0 1]);
[stabmarg_n,wcu_n] = robstab(T_n);

%% Task 2.4 Question: Robust stability with parametric uncertainty
% Non-minimum phase system
[stabmarg_param_n,wcu_param_n] = robstab(G_u_cl_n);
