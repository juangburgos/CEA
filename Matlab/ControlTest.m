clear all
clc

%% Load Data

load('experiment.mat');
load('FitsParaJUAN_k.mat');
load('model.mat');
Ts = experiment.time(2)-experiment.time(1);
ruta=[12 0;12 5;7 -5;-1 4;-6 -5;-2 0];
radio=0.05;

%% Init Simulation Parameters
% Disturbances
dist             = zeros(size(dsys.c,1),max(size(time)));
% dist(1,:)        = feval(fitresultDISTURBANCES{1},time);
% dist(3,:)        = feval(fitresultDISTURBANCES{2},time);
% dist(2,:)        = feval(fitresultvelDISTURBANCES{1},time);
% dist(4,:)        = feval(fitresultvelDISTURBANCES{2},time);

% Noise
namp = 0.04;
noise      = zeros(size(dsys.c,1),max(size(time)));
% noise      = namp.*rand(size(dsys.c,1),max(size(time)))-(namp/2);


%% Unconstrained MPC Control with Optimal Robust Disturbance Model and Observer Gain

% Measure Robust LQR Computational Time
tic;

% Init Simulation Vectors
u    = zeros(size(dsys.b,2),max(size(time))); % Input
y    = zeros(size(dsys.c,1),max(size(time))); % Observer Output
xp   = zeros(size(dsys.a,1),max(size(time))); % Plant State
yp   = zeros(size(dsys.c,1),max(size(time))); % Plant Output

for t = 2:1:max(size(time))
    
    Bsch = dsys.b*[feval(fitresultROLL_k,refroll(t-1)) , 0; ...
                0, feval(fitresultPITCH_k,refpitch(t-1))];
               
    % ************** UPDATE PLANT SIMULATION ****************
    % State Equation
    xp(:,t) = dsys.a * xp(:,t-1) + Bsch * u(:,t-1); 
    % Output Equation
    yp(:,t)   = dsys.c * xp(:,t);
    % Additive Disturbance
    yp(:,t)   = yp(:,t) + dist(:,t);
    % Measurement Noise
    yp(:,t)   = yp(:,t) + noise(:,t);
    
    % ************** CONTROLLER *******************************
    % Controller goes here yp(:,t)
    [y(:,t),u(:,t)] = MatController(yp(:,t));
    
    clc
    disp(['Unconstrained MPC Progress: ',num2str(100*(t/max(size(time)))),' %']);
    
end

% Measure Robust LQR Computational Time
trobust = toc/time(end);


figure();
subplot(6,1,1),plot(experiment.time,r(1,:),'b-');
hold on;
subplot(6,1,1),plot(experiment.time,yp(1,:),'r--');
axis([experiment.time(1) experiment.time(end) 0.9*min(yp(1,:)) 1.1*max(yp(1,:))]);

subplot(6,1,2),plot(experiment.time,r(2,:),'b-');
hold on;
subplot(6,1,2),plot(experiment.time,yp(2,:),'r--');
axis([experiment.time(1) experiment.time(end) 0.9*min(yp(2,:)) 1.1*max(yp(2,:))]);

subplot(6,1,3),plot(experiment.time,r(3,:),'b-');
hold on;
subplot(6,1,3),plot(experiment.time,yp(3,:),'r--');
axis([experiment.time(1) experiment.time(end) 0.9*min(yp(3,:)) 1.1*max(yp(3,:))]);

subplot(6,1,4),plot(experiment.time,r(4,:),'b-');
hold on;
subplot(6,1,4),plot(experiment.time,yp(4,:),'r--');
axis([experiment.time(1) experiment.time(end) 0.9*min(yp(4,:)) 1.1*max(yp(4,:))]);

subplot(6,1,5),plot(experiment.time,u(1,:),'b-');
axis([experiment.time(1) experiment.time(end) -1 1]);

subplot(6,1,6),plot(experiment.time,u(2,:),'b-');
axis([experiment.time(1) experiment.time(end) -1 1]);

clc
disp(['Unconstrained MPC: ',num2str(100*trobust),' %']);