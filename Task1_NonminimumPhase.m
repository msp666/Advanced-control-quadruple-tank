%% Task 1: System analysis (Non-minimum Phase System)
%% 1.2 Question 2: System representation and analysis
%% Non-minimum phase

% Transfer function
G_tf_n = tf(minreal(G_n));
Gtf_zp_n = zpk(minreal(G_tf_n));
% Poles and zeros
Poles_n = pole(G_n);
Zeros_n = tzero(G_n);
% Associated pole and zero directions
N = length(Poles_n);
pole_directions_n = zeros(2,2*N);
for i = 1:length(Poles_n)
    G_p_n = evalfr(G_tf_n,Poles_n(i)+1e-6);
    [U_p_n,W_p_n,V_p_n] = svd(G_p_n);
    pole_directions_n(:,2*i-1) = U_p_n(:,1);
    pole_directions_n(:,2*i) = V_p_n(:,1);
end

M = length(Zeros_n);
zero_directions_n = zeros(2,2*M);
for j = 1:length(Zeros_n)
    G_z_n = evalfr(G_tf_n,Zeros_n(j));
    [U_z_n,W_z_n,V_z_n] = svd(G_z_n);
    zero_directions_n(:,2*j-1) = U_z_n(:,1);
    zero_directions_n(:,2*j) = V_z_n(:,1);
end

%% 1.3 Question 3:Relative Gain Array (RGA)
% Non-minimum phase: Calculation of G_n(0)'s RGA  compare,if RGA_G_0 equals to Matrix_lamda_n
G_0_n = evalfr(G_tf_n,0);
RGA_G_0_n = G_0_n.* pinv(G_0_n).';

%Compare, if RGA_G_0_n equals to matrix RGA_n
g1_n = par(2).g1;
g2_n = par(2).g2;
lamda_n = g1_n*g2_n/(g1_n + g2_n -1);
RGA_n = [lamda_n (1-lamda_n);(1-lamda_n) lamda_n];

%% 1.4 Question 4: Parametric uncertainty for nonlinearities

%% 1.5 Question 5: Singular values
% Plot Singular Values of Non-minimum Phase
sigma(G_n);
title('Singular Values for non-minimum phase system');
grid on;

sigma(G_n);
title('Singular Values for non-minimum phase system');
grid on;