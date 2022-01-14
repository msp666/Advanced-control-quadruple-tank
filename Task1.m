%% Task 1: System analysis
%% 1.2 Question 2: System representation and analysis

% A = [-a1/A1*sqrt(g/(2*h10)), 0, a3/A1*sqrt(g/(2*h30)), 0;
%     0, -a2/A2*sqrt(g/(2*h20)), 0, a4/A2*sqrt(g/(2*h40));
%     0, 0, -a3/A3*sqrt(g/(2*h30)), 0;
%     0, 0, 0, -a4/A4*sqrt(g/(2*h40))                          ];
% B = [g1*k1/A1, 0; 0, g2*k2/A2; 0, (1-g2)*k2/A3; (1-g1)*k1/A4, 0];
% C = [1 0 0 0; 0 1 0 0];
% D = [0 0;0 0];

% G = ss(A,B,C,D)

% Transfer function
G1 = tf(minreal(G));
G1_zp = zpk(minreal(G1));

% Poles and zeros
Poles = pole(G);
Zeros = tzero(G);

% Associated pole and zero directions
N = length(Poles);
pole_directions = zeros(2,2*N);
for i = 1:length(Poles)
    G_p = evalfr(G1,Poles(i)+1e-6);
    [U_p,W_p,V_p] = svd(G_p);
    pole_directions(:,2*i-1) = U_p(:,1);
    pole_directions(:,2*i) = V_p(:,1);
end


M = length(Zeros);
zero_directions = zeros(2,2*M);
for j = 1:length(Zeros)
    G_z = evalfr(G1,Zeros(j));
    [U_z,W_z,V_z] = svd(G_z);
    zero_directions(:,2*j-1) = U_z(:,1);
    zero_directions(:,2*j) = V_z(:,1);
end


%% 1.3 Question 3:Relative Gain Array (RGA)
% Calculation of G(0)'s RGA
G_0 = evalfr(G1,0);
RGA_G_0 = G_0.* pinv(G_0).';

% Compare, if RGA_G_0 equals to Matrix_lamda
lamda_m = g1*g2/(g1+g2-1);
RGA_m = [lamda_m (1-lamda_m);(1-lamda_m) lamda_m];

% Non-minimum phase condition & compare,if RGA_G_0 equals to Matrix_lamda_n
g1_n = par(2).g1;
g2_n = par(2).g2;
lamda_n = g1_n*g2_n/(g1_n + g2_n -1);
RGA_n = [lamda_n (1-lamda_n);(1-lamda_n) lamda_n];

% Dynamic RGA
lamda = 1/(1-G1(1,2)*G1(2,1)/(G1(1,1)*G1(2,2)));
RGA_dynamic = [lamda (1-lamda);(1-lamda) lamda];
figure;
bode(RGA_dynamic);
grid on;
title('dynamic RGA')

%% 1.5 Question 5: Singular values
figure;
sigma(G);





