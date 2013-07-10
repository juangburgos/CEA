clear all
close all
clc

%% Load Data

load('experiment.mat');
load('FitsParaJUAN_k.mat');
load('model.mat');
load('fitresultDIST.mat');
load('fitresultvelDISTURBANCES.mat');
Ts = experiment.time(2)-experiment.time(1);
ruta=[12 0;12 5;7 -5;-1 4;-6 -5;-2 0];
% ruta=[12 0;12 5;7 -5;-1 4;-1 4;-2 0];
radio=0.05;

rutaini  = [0,0];
ruta     = [rutaini;ruta];

%% Init Simulation Parameters
% Storage Size 120s
largesize = 120/Ts;

time = 0:Ts:Ts*(largesize-1)-Ts;

% Disturbances
% dist             = zeros(size(dsys.c,1),largesize);
dist(1,:)        = feval(fitresultDISTURBANCES{1},time);
dist(3,:)        = feval(fitresultDISTURBANCES{2},time);
dist(2,:)        = feval(fitresultvelDISTURBANCES{1},time);
dist(4,:)        = feval(fitresultvelDISTURBANCES{2},time);

% Noise
% noise      = zeros(size(dsys.c,1),largesize);
namp       = 0.04;
noise      = namp.*rand(size(dsys.c,1),largesize)-(namp/2);


%% Unconstrained MPC Matlab

% Init Simulation Vectors
u    = zeros(size(dsys.b,2),largesize); % Input
y    = zeros(size(dsys.c,1),largesize); % Observer Output
xp   = zeros(size(dsys.a,1),largesize); % Plant State
yp   = zeros(size(dsys.c,1),largesize); % Plant Output

% Measure Computational Time
meast_mat   = 0; % Real Time

t = 2;
while 1
    
    Bsch = dsys.b*[feval(fitresultROLL_k,u(1,t-1)) , 0; ...
                0, feval(fitresultPITCH_k,u(2,t-1))];
               
    % ************** UPDATE PLANT SIMULATION ****************
    % State Equation
    xp(:,t) = dsys.a * xp(:,t-1) + Bsch * u(:,t-1); 
    % Output Equation
    yp(:,t)   = dsys.c * xp(:,t);
    % Additive Disturbance
    yp(:,t)   = yp(:,t) + dist(:,t);
    % Measurement Noise
    yp(:,t)   = yp(:,t) + noise(:,t);
    
    % ************** CONTROLLER ******************************
    tic;
    [y(:,t),u(:,t)] = MatController(yp(:,t),ruta(:,1),ruta(:,2),size(ruta,1),radio);
    meast_mat = meast_mat + toc;
    
    clc
    disp(['Unconstrained MPC Progress: ',num2str(100*(t/largesize)),' %']);
    
    if norm([yp(1,t),yp(3,t)]-ruta(end,:)) < (radio)
        break;
    end
    
    t = t + 1;
    
end

% Measure Computational Time
trobmat = meast_mat/(Ts*(t-1));

% Cut Storage size
u    = u(:,1:t-1); % Input
y    = y(:,1:t-1); % Observer Output
xp   = xp(:,1:t-1); % Plant State
yp   = yp(:,1:t-1); % Plant Output
time = time(:,1:t-1);

figure();
subplot(6,1,1),plot(time,y(1,:),'b-');
hold on;
subplot(6,1,1),plot(time,yp(1,:),'r--');
axis([time(1) time(end) 0.9*min(yp(1,:)) 1.1*max(yp(1,:))]);
legend('Observer Output','Plant Output');

subplot(6,1,2),plot(time,y(2,:),'b-');
hold on;
subplot(6,1,2),plot(time,yp(2,:),'r--');
axis([time(1) time(end) 0.9*min(yp(2,:)) 1.1*max(yp(2,:))]);

subplot(6,1,3),plot(time,y(3,:),'b-');
hold on;
subplot(6,1,3),plot(time,yp(3,:),'r--');
axis([time(1) time(end) 0.9*min(yp(3,:)) 1.1*max(yp(3,:))]);

subplot(6,1,4),plot(time,y(4,:),'b-');
hold on;
subplot(6,1,4),plot(time,yp(4,:),'r--');
axis([time(1) time(end) 0.9*min(yp(4,:)) 1.1*max(yp(4,:))]);

subplot(6,1,5),plot(time,u(1,:),'b-');
axis([time(1) time(end) -1 1]);

subplot(6,1,6),plot(time,u(2,:),'b-');
axis([time(1) time(end) -1 1]);

figure();
plot(ruta(:,1),ruta(:,2),'r-');
hold on;
plot(yp(1,:),yp(3,:),'b-');

clc
disp(['Unconstrained MPC Matlab: ',num2str(100*trobmat),' %']);

%% Unconstrained MPC in C

time = 0:Ts:Ts*(largesize-1)-Ts;

% Load DLL
dllname    = 'libCEA.dll';
headername = 'mydll.h';
funcon    = 'matcontroller';
if ~libisloaded( 'MYDLL' ) 
    loadlibrary( dllname, headername, 'alias', 'MYDLL' );      
end

% Init Simulation Vectors
u_c    = zeros(size(dsys.b,2),largesize); % Input
y_c    = zeros(size(dsys.c,1),largesize); % Observer Output
xp_c   = zeros(size(dsys.a,1),largesize); % Plant State
yp_c   = zeros(size(dsys.c,1),largesize); % Plant Output

% Measure Computational Time
meast_c   = 0; % Real Time

t = 2;
while 1
    
    Bsch = dsys.b*[feval(fitresultROLL_k,u_c(1,t-1)) , 0; ...
                0, feval(fitresultPITCH_k,u_c(2,t-1))];
               
    % ************** UPDATE PLANT SIMULATION ****************
    % State Equation
    xp_c(:,t) = dsys.a * xp_c(:,t-1) + Bsch * u_c(:,t-1); 
    % Output Equation
    yp_c(:,t)   = dsys.c * xp_c(:,t);
    % Additive Disturbance
    yp_c(:,t)   = yp_c(:,t) + dist(:,t);
    % Measurement Noise
    yp_c(:,t)   = yp_c(:,t) + noise(:,t);
    
    % ************** CONTROLLER ******************************
    % matcontroller(doublereal *u, doublereal *y, doublereal *ym, doublereal *wayPointX, doublereal *wayPointY, ...
    %               int numWaypoints, doublereal *radio);
    % Observer (Call DLL)
    p_u       = libpointer('doublePtr',u_c(:,t-1));
    p_ym      = libpointer('doublePtr',yp_c(:,t));
    p_y       = libpointer('doublePtr',y_c(:,t));
    p_wpx     = libpointer('doublePtr',ruta(:,1));
    p_wpy     = libpointer('doublePtr',ruta(:,2));
    p_rad     = libpointer('doublePtr',radio);
    tic;
    calllib( 'MYDLL', funcon, p_u, p_y, p_ym, p_wpx, p_wpy, size(ruta,1), p_rad);
    % Measure Computational Time
    meast_c = meast_c + toc;
    u_c(:,t)        = get(p_u,'Value');
    y_c(:,t)        = get(p_y,'Value');
    
    clc
    disp(['Unconstrained MPC Progress: ',num2str(100*(t/largesize)),' %']);
    
    if norm([yp_c(1,t),yp_c(3,t)]-ruta(end,:)) < (radio)
        break;
    end
    
    t = t + 1;
    
end

trobc = meast_c/(Ts*(t-1));

% Unload DLL
unloadlibrary('MYDLL');

% Cut Storage size
u_c    = u_c(:,1:t-1); % Input
y_c    = y_c(:,1:t-1); % Observer Output
xp_c   = xp_c(:,1:t-1); % Plant State
yp_c   = yp_c(:,1:t-1); % Plant Output
time   = time(:,1:t-1);

figure();
subplot(6,1,1),plot(time,y_c(1,:),'b-');
hold on;
subplot(6,1,1),plot(time,yp_c(1,:),'r--');
axis([time(1) time(end) 0.9*min(yp_c(1,:)) 1.1*max(yp_c(1,:))]);
legend('Observer Output','Plant Output');

subplot(6,1,2),plot(time,y_c(2,:),'b-');
hold on;
subplot(6,1,2),plot(time,yp_c(2,:),'r--');
axis([time(1) time(end) 0.9*min(yp_c(2,:)) 1.1*max(yp_c(2,:))]);

subplot(6,1,3),plot(time,y_c(3,:),'b-');
hold on;
subplot(6,1,3),plot(time,yp_c(3,:),'r--');
axis([time(1) time(end) 0.9*min(yp_c(3,:)) 1.1*max(yp_c(3,:))]);

subplot(6,1,4),plot(time,y_c(4,:),'b-');
hold on;
subplot(6,1,4),plot(time,yp_c(4,:),'r--');
axis([time(1) time(end) 0.9*min(yp_c(4,:)) 1.1*max(yp_c(4,:))]);

subplot(6,1,5),plot(time,u_c(1,:),'b-');
axis([time(1) time(end) -1 1]);

subplot(6,1,6),plot(time,u_c(2,:),'b-');
axis([time(1) time(end) -1 1]);

figure();
plot(ruta(:,1),ruta(:,2),'r-');
hold on;
plot(yp_c(1,:),yp_c(3,:),'b-');

clc
disp(['Unconstrained MPC Matlab: ',num2str(100*trobmat),' %']);
disp(['Unconstrained MPC C: ',num2str(100*trobc),' %']);