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

% set up flag for different task
Ad30 = 0.08;
Ad40 = 0.08;

W_act = 0.2*eye(2);
Delta_act = ultidyn('Delta_act', [2, 2], 'Bound', 1);
vals = struct;
Task4_flag = 0;
noise_setting = 0;
controller_select = 1;

% Minimum phase
par(1) = struct('h10',12.1,'h20',12.6,'h30',2.5,'h40',2.49,'v10',2.99,'v20',2.97,'k1',3.33, 'k2', 3.35, 'g1',0.6,'g2',0.5);
% Nonminimum phase
par(2) = struct('h10',6.79,'h20',8.78,'h30',2.97,'h40',4.17,'v10',2.53,'v20',2.35,'k1',3.14, 'k2', 3.29, 'g1',0.35,'g2',0.3);

% Minimum phase
mode = 1;
h10 = par(mode).h10;
h20 = par(mode).h20;
h30 = par(mode).h30;
h40 = par(mode).h40;
k1 = par(mode).k1;
k2 = par(mode).k2;
g1 = par(mode).g1;
g2 = par(mode).g2;
op_spec = operspec(model);
op_spec.Inputs(1).u = [par(mode).v10; par(mode).v20]; 
% Allow only positive voltages
op_spec.Inputs(1).Min = zeros(2,1); 
% Define states
set(op_spec.State(1), 'Known', false);
set(op_spec.State(1), 'x', h10);
set(op_spec.State(2), 'Known', false);
set(op_spec.State(2), 'x', h20);
% set(op_spec.State(2), 'Known', false); 
set(op_spec.State(3), 'x', h30);
% set(op_spec.State(3), 'Known', false); 
set(op_spec.State(4), 'x', h40);
% set(op_spec.State(4), 'Known', false); 

[op, op_report] = findop(model, op_spec);    
G = linearize(model,getlinio(model),op(1));
U0 = op_report.Inputs(1).u;
X0 = [op_report.States(1).x; op_report.States(2).x; op_report.States(3).x; op_report.States(4).x];
Y0 = op_report.Outputs(1).y;

