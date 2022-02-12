%% results of symbolic linearisation
clearvars
model = 'nlQuadtank';
% Cross sectional area of tank 
A1=30; % [cm^2]
A2=35; % [cm^2]
A3=30; % [cm^2]
A4=35; % [cm^2]

% Cross sectional area of outlet
a1=0.071; % [cm^2]
a2=0.057; % [cm^2]
a3=0.071; % [cm^2]
a4=0.057; % [cm^2]

% Graviation force
g=981; % [cm/s^2]

% Minimum phase
par(1) = struct('h10',12.1,'h20',12.6,'h30',2.5,'h40',2.49,'v10',2.99,'v20',2.97,'k1',3.33, 'k2', 3.35, 'g1',0.6,'g2',0.5);
h10 = par(1).h10;
h20 = par(1).h20;
h30 = par(1).h30;
h40 = par(1).h40;
k1 = par(1).k1;
k2 = par(1).k2;
g1 = par(1).g1;
g2 = par(1).g2;
A_mp = [ -a1/A1*sqrt( (g/(2*h10)) ), 0, a3/A1*sqrt( (g/(2*h30)) ), 0;...
         0, -a2/A2*sqrt( (g/(2*h20)) ), 0, a4/A2*sqrt( (g/(2*h40)) );...
         0, 0, -a3/A3*sqrt( (g/(2*h30)) ), 0;...
         0, 0, 0, -a4/A4*sqrt( (g/(2*h40)) ) ]
B_mp = [ (g1*k1)/A1, 0;...   
         0, (g2*k2)/A2;...
         0, ((1-g2)*k2)/A3;...
         ((1-g1)*k1)/A4, 0 ]  


% Nonminimum phase
par(2) = struct('h10',6.79,'h20',8.78,'h30',2.97,'h40',4.17,'v10',2.53,'v20',2.35,'k1',3.14, 'k2', 3.29, 'g1',0.35,'g2',0.3);
h10 = par(2).h10;
h20 = par(2).h20;
h30 = par(2).h30;
h40 = par(2).h40;
k1 = par(2).k1;
k2 = par(2).k2;
g1 = par(2).g1;
g2 = par(2).g2;
A_np = [ -a1/A1*sqrt( (g/(2*h10)) ), 0, a3/A1*sqrt( (g/(2*h30)) ), 0;...
         0, -a2/A2*sqrt( (g/(2*h20)) ), 0, a4/A2*sqrt( (g/(2*h40)) );...
         0, 0, -a3/A3*sqrt( (g/(2*h30)) ), 0;...
         0, 0, 0, -a4/A4*sqrt( (g/(2*h40)) ) ]
B_np = [ (g1*k1)/A1, 0;...   
         0, (g2*k2)/A2;...
         0, ((1-g2)*k2)/A3;...
         ((1-g1)*k1)/A4, 0 ]
