%% Task 2: Decentralised control and uncertainty modeling(Minimum Phase System)
% Calculate Niederlinski Index(NI)
% For minimum phase system
G_tf = tf(minreal(G));
G_0 = evalfr(G_tf,0);
G_0_det = det(G_0);
NI1 = G_0_det/(G_0(1,1)*G_0(2,2)); 
NI2 = -G_0_det/(G_0(1,2)*G_0(2,1));

%% Task 2.1 Question 1: Decentralised controller based on linear and nonlinear models
% Minimum phase system
s = tf('s');
Tao1 = A1/a1*sqrt(2*h10/g);
Tao2 = A2/a2*sqrt(2*h20/g);
Tao3 = A3/a3*sqrt(2*h30/g);
Tao4 = A4/a4*sqrt(2*h40/g);
c1 = Tao1*k1/A1;
c2 = Tao2*k2/A2;

% Boundary values for K and T
K1_B = (9.2*Tao1 - 45)/(45*g1*c1);
T1_B = 4*Tao1*K1_B*g1*c1*(log(10))^2/((K1_B*g1*c1+1)^2*(pi^2 + (log(10))^2));
K2_B = (9.2*Tao2 - 45)/(45*g2*c2);
T2_B = 4*Tao2*K2_B*g2*c2*(log(10))^2/((K2_B*g2*c2+1)^2*(pi^2 + (log(10))^2));

% Choose K1,T1,K2,T2
K1 = 5;
K2 = 6;
T1 = 20;
T2 = 25;

% Calculation of closed-loop transfer function, plot and display
Gtf = [g1*c1/(1+Tao1*s)                  (1-g2)*c1/((1+Tao3*s)*(1+Tao1*s));
      (1-g1)*c2/((1+Tao4*s)*(1+Tao2*s))     g2*c2/(1+Tao2*s)               ];
K = [K1*(1+1/(T1*s)), 0; 0, K2*(1+1/(T2*s))];
Gtf_cl = feedback(G*K,[1 0; 0 1]);
figure;
step(Gtf_cl);
grid on;
stepinfo(Gtf_cl);

%% Task 2.2 Question 2: Model uncertainty and weighting function
% Minimum Phase
% Uncertain parameters of the model
g_u1 = ureal('g_u1',g1,'Percentage',[-10,10]);
g_u2 = ureal('g_u2',g2,'Percentage',[-10,10]);
k_u1 = ureal('k_u1',k1,'Percentage',[-10,10]);
k_u2 = ureal('k_u2',k2,'Percentage',[-10,10]);
a_u1 = ureal('a_u1',a1,'Percentage',[-5,5]);
a_u2 = ureal('a_u2',a2,'Percentage',[-5,5]);
a_u3 = ureal('a_u3',a3,'Percentage',[-5,5]);
a_u4 = ureal('a_u4',a4,'Percentage',[-5,5]);

% Set up state space
A_u = [-a_u1/A1*sqrt(g/(2*h10)),                          0,   a_u3/A3*sqrt(g/(2*h30)),                         0;
                              0,   -a_u2/A2*sqrt(g/(2*h20)),                         0,   a_u4/A4*sqrt(g/(2*h40));
                              0,                          0,  -a_u3/A3*sqrt(g/(2*h30)),                         0;
                              0,                          0,                         0,   -a_u4/A4*sqrt(g/(2*h40))];
B_u = [   g_u1*k_u1/A1,                  0;
                     0,       g_u2*k_u2/A2;
                     0,   (1-g_u2)*k_u2/A3;
      (1-g_u1)*k_u1/A4,                  0];
C_u = G.C;
D_u = G.D;
G_u = ss(A_u, B_u, C_u, D_u);

% Close-loop transfer function
G_u_cl = feedback(G_u*K, [1 0; 0 1]);

% Set up frequency vector
w = logspace(-3, 3, 100);
% Define frequency response buffer
G_mag = zeros(2, 2, length(w), 1500);

for k = 1:1500
    % Select the parameters in the uncetain interval through uniform distribution
    p1 = (rand(1,4) * 0.2 + 0.9) .* [g1, g2, k1, k2];
    p2 = (rand(1,4) * 0.1 + 0.95) .* [a1, a2, a3, a4];
    
    % Parateter definition,and set up state space
    A = [-p2(1)/A1*sqrt(g/(2*h10)),                             0,      p2(3)/A3*sqrt(g/(2*h30)),                              0;
                                 0,     -p2(2)/A2*sqrt(g/(2*h20)),                              0,      p2(4)/A4*sqrt(g/(2*h40));
                                 0,                             0,      -p2(3)/A3*sqrt(g/(2*h30)),                             0;
                                 0,                             0,                              0,      -p2(4)/A4*sqrt(g/(2*h40))];
    B = [   p1(1)*p1(3)/A1,                    0;
                         0,       p1(2)*p1(4)/A2;
                         0,   (1-p1(2))*p1(4)/A3;
         (1-p1(1))*p1(3)/A4,                   0];
    C = G.C;
    D = G.D;
    G_tmp = ss(A, B, C, D);
    
    % Get magnitude
    [Gp_mag(:, : ,:, k)] = squeeze(bode(inv(G) * (G_tmp - G), w)); 
end

% Get Maximum Magnitude of each frequency 
Gp_mag_max = max(Gp_mag, [], 4);
FR = frd(Gp_mag_max, w);

% Fit weighting function,plot and display
C11.LowerBound = FR(1,1,:); C11.UpperBound = [];
C12.LowerBound = FR(1,2,:); C12.UpperBound = [];
C21.LowerBound = FR(2,1,:); C21.UpperBound = [];
C22.LowerBound = FR(2,2,:); C22.UpperBound = [];
WM11 = fitmagfrd(FR(1, 1), 2, [], [], C11);
WM12 = fitmagfrd(FR(1, 2), 2, [], [], C12);
WM21 = fitmagfrd(FR(2, 1), 2, [], [], C21);
WM22 = fitmagfrd(FR(2, 2), 2, [], [], C22);
WM = [WM11, WM12; WM21, WM22];
figure;
bodemag(1/WM ,  '-*r' ); grid on; hold on;
title('Weighting of dynamic uncertainty for minimum phase system')

%% Task 2.3 Question: Robust stability with dynamic uncertainty
% Minimum Phase

DeltaM = ultidyn('DeltaM', [2, 2], 'Bound', 1);
% Additive uncertainty
Gp = G + WM * DeltaM;
%Gp = G*(eye(2) + WM*DeltaM);
L = Gp * K;
T = feedback(L ,[1 0;0 1]);
[stabmarg,wcu] = robstab(T);

%% Task 2.4 Question: Robust stability with parametric uncertainty
% Minimum Phase
[stabmarg_param,wcu_param] = robstab(G_u_cl);