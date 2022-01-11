%% Given INCOMPLETE solution script for the quadtank coursework
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
% Nonminimum phase
par(2) = struct('h10',6.79,'h20',8.78,'h30',2.97,'h40',4.17,'v10',2.53,'v20',2.35,'k1',3.14, 'k2', 3.29, 'g1',0.35,'g2',0.3);

% Minimum phase
h10 = par(1).h10;
h20 = par(1).h20;
h30 = par(1).h30;
h40 = par(1).h40;
k1 = par(1).k1;
k2 = par(1).k2;
g1 = par(1).g1;
g2 = par(1).g2;
op_spec = operspec(model);
op_spec.Inputs(1).u = [par(1).v10; par(1).v20];  
% Allow only positive voltages
op_spec.Inputs(1).Min = zeros(2,1); 
% Define states
set(op_spec.State(1), 'x', 12.1);
set(op_spec.State(1), 'Known', true);    
[op, op_report] = findop(model, op_spec);    
G = linearize(model,getlinio(model),op(1));
U0 = op_report.Inputs(1).u;
X0 = [op_report.States(1).x; op_report.States(2).x; op_report.States(3).x; op_report.States(4).x];
Y0 = op_report.Outputs(1).y;