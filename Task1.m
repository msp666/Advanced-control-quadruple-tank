%% Task 1
%% T1.2 TF, Poles, Zeros and directions
% transfer function
Gtf = tf(minreal(G))

% poles and zeros
Gpoles = pole(G);
Gzeros = tzero(minreal(G));

% direction of poles and zeros
poles_direction = cell(length(Gpoles), 2);
for i = 1:length(Gpoles)
    Gtf_tmp = evalfr(Gtf, Gpoles(i) + 1e-5);
    [U, S, V] = svd(Gtf_tmp);
    % Input direction
    poles_direction{i, 1} = V(:, 1);
    % Output direction
    poles_direction{i, 2} = U(:, 2);  
end

zeros_direction = cell(length(Gzeros), 2);
for i = 1: length(Gzeros)
    Gtf_tmp = evalfr(Gtf, Gzeros(i));
    [U, S, V] = svd(Gtf_tmp);
    % Input direction
    zeros_direction{i, 1} = V(:, end);
    % Output direction
    zeros_direction{i, 2} = U(:, end);
end 

%% T1.3 RGA
% RGA for s = 0
G0 = [Gtf.Numerator{1, 1}(end) / Gtf.Denominator{1, 1}(end), Gtf.Numerator{1, 2}(end) / Gtf.Denominator{1, 2}(end);
      Gtf.Numerator{2, 1}(end) / Gtf.Denominator{2, 1}(end), Gtf.Numerator{2, 2}(end) / Gtf.Denominator{2, 2}(end)];
RGA0 = G0 .* pinv(G0).';
lamda0 = g1 * g2 / (g1 + g2 - 1);

% RGA over frequency range
lamda = 1/(1-Gtf(1,2)*Gtf(2,1)/(Gtf(1,1)*Gtf(2,2))); 
RGA = zpk([lamda,   1-lamda;
       1-lamda, lamda  ]);
figure;
bode(RGA);grid on;
title('RGA magnitude')

%% T1.5
% plot the singular value over frequency range
figure
sigma(Gtf); grid on
% bode digaramm
figure 
bode(Gtf);grid on
title('Bode diagramm of plant')
figure
subplot(2,2,1); margin(Gtf(1,1)); grid on
subplot(2,2,2); margin(Gtf(1,2)); grid on
subplot(2,2,3); margin(Gtf(2,1)); grid on
subplot(2,2,4); margin(Gtf(2,2)); grid on