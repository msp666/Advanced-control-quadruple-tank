%% Task 1: System analysis (Minimum Phase System)
%% 1.2 Question 2: System representation and analysis
% Minimum phase
% Transfer function
G_tf = tf(minreal(G));
Gtf_zp = zpk(minreal(G_tf));
% Poles and zeros
Poles = pole(G);
Zeros = tzero(G);
% Associated pole and zero directions
N = length(Poles);
pole_directions = zeros(2,2*N);
for i = 1:length(Poles)
    G_p = evalfr(G_tf,Poles(i)+1e-6);
    [U_p,W_p,V_p] = svd(G_p);
    pole_directions(:,2*i-1) = U_p(:,1);
    pole_directions(:,2*i) = V_p(:,1);
end

M = length(Zeros);
zero_directions = zeros(2,2*M);
for j = 1:length(Zeros)
    G_z = evalfr(G_tf,Zeros(j));
    [U_z,W_z,V_z] = svd(G_z);
    zero_directions(:,2*j-1) = U_z(:,1);
    zero_directions(:,2*j) = V_z(:,1);
end


%% 1.3 Question 3:Relative Gain Array (RGA)
% Minimum Phase:Calculation of G(0)'s RGA
G_0 = evalfr(G_tf,0);
RGA_G_0 = G_0.* pinv(G_0).';

% Compare, if RGA_G_0 equals to matrix RGA_m
lamda_m = g1*g2/(g1+g2-1);
RGA_m = [lamda_m (1-lamda_m);(1-lamda_m) lamda_m];

% Dynamic RGA
lamda = 1/(1-G_tf(1,2)*G_tf(2,1)/(G_tf(1,1)*G_tf(2,2)));
RGA_dynamic = [lamda (1-lamda);(1-lamda) lamda];
figure;
bode(RGA_dynamic);
grid on;
title('Dynamic RGA for minimum phase system')

%% 1.4 Question 4: Parametric uncertainty for nonlinearities

%% 1.5 Question 5: Singular values
% Plot Singular Values of Minimum Phase
sigma(G);
title('Singular Values for minimum phase system');
grid on;
