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
% ruta=[12 0;12 5;7 -5;-1 4;-6 -5;-2 0];
% ruta=[12 0;12 5;7 -5;-1 4;-1 4;-2 0]; 
load('ruta_Eval.mat');

vxmax    = 2.4;                     % max 3.5 [m/s] ~3.1
axmax    = 2.4;                       % ~1 (dont move)
jxmax    = 8;                      % ~12 (dont move)
vymax    = 0.9;                     % max 2.1 [m/s] ~1.5 (dont move)
aymax    = 0.6;                    % ~0.55 (dont move)
jymax    = 0.9;                       % ~6 (dont move)

% % vxmax    = 1;                     % max 2.1 [m/s] ~1.5 (dont move)
% % axmax    = 0.3;                    % ~0.55 (dont move)
% % jxmax    = 0.3;                       % ~6 (dont move)
% % vymax    = 2;                     % max 3.5 [m/s] ~3.1
% % aymax    = 0.6;                       % ~1 (dont move)
% % jymax    = 0.6;                      % ~12 (dont move)

% Radio
radio = 0.1; % nuevo radio

rutaini  = [0,0];
ruta     = [rutaini;ruta];

% Simulate extra time
extrat   = 200;

%% Init Simulation Parameters
% Storage Size 120s
largesize = round(220/Ts);

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
xp   = zeros(size(dsys.a,1),largesize); % Plant State
yp   = zeros(size(dsys.c,1),largesize); % Plant Output

% Measure Computational Time
meast_mat   = 0; % Real Time

t = 2;
extratime = 0;
res = 0;
while 1
    
    Bsch = dsys.b*[feval(fitresultPITCH_k,u(1,t-1)) , 0; ...
                   0, feval(fitresultROLL_k,u(2,t-1))];
               
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
    position        = [yp(1,t),yp(3,t),0];
    velocity        = [yp(2,t),yp(4,t),0];
    numAxis         = 3;
    wayPointX       = ruta(:,1);
    wayPointY       = ruta(:,2);
    numWaypoints    = size(ruta,1);
    actualWayPoint  = [yp(1,t),yp(3,t)];
    param           = [vxmax; axmax; jxmax; vymax; aymax; jymax ];
    numParam        = length(param);
    
    tic;
    u(:,t) = MatController(position,velocity,numAxis,wayPointX, ...
            wayPointY,numWaypoints,actualWayPoint,param,numParam);      
    meast_mat = meast_mat + toc;
    
    clc
    disp(['Unconstrained MPC Progress: ',num2str(100*(t/largesize)),' %']);
    
    if res+2 <= length(ruta)
        if norm([yp(1,t),yp(3,t)]-ruta(res+2,:)) < (radio)
            res = res + 1;
        end
    else
        if extratime == extrat
            break;
        end
        extratime = extratime + 1;
    end
    
    t = t + 1;
    
end

% Measure Computational Time
trobmat = meast_mat/(Ts*(t-1));

% Cut Storage size
u    = u(:,1:t-1); % Input
xp   = xp(:,1:t-1); % Plant State
yp   = yp(:,1:t-1); % Plant Output
time = time(:,1:t-1);

figure();
subplot(6,1,1),plot(time,yp(1,:),'r--');
axis([time(1) time(end) 0.9*min(yp(1,:)) 1.1*max(yp(1,:))]);
legend('Plant Output');

subplot(6,1,2),plot(time,yp(2,:),'r--');
axis([time(1) time(end) 0.9*min(yp(2,:)) 1.1*max(yp(2,:))]);

subplot(6,1,3),plot(time,yp(3,:),'r--');
axis([time(1) time(end) 0.9*min(yp(3,:)) 1.1*max(yp(3,:))]);

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
funcon    = 'Control';
if ~libisloaded( 'MYDLL' ) 
    loadlibrary( dllname, headername, 'alias', 'MYDLL' );      
end

% Init Simulation Vectors
u_c    = zeros(size(dsys.b,2),largesize); % Input
xp_c   = zeros(size(dsys.a,1),largesize); % Plant State
yp_c   = zeros(size(dsys.c,1),largesize); % Plant Output

% Measure Computational Time
meast_c   = 0; % Real Time

t = 2;
extratime = 0;
res = 0;
while 1
    
    Bsch = dsys.b*[feval(fitresultPITCH_k,u_c(1,t-1)) , 0; ...
                0, feval(fitresultROLL_k,u_c(2,t-1))];
               
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
    position        = [yp_c(1,t),yp_c(3,t),0];
    velocity        = [yp_c(2,t),yp_c(4,t),0];
    numAxis         = 3;
    wayPointX       = ruta(:,1);
    wayPointY       = ruta(:,2);
    numWaypoints    = size(ruta,1);
    actualWayPoint  = [yp_c(1,t),yp_c(3,t)];
    param           = [vymax; aymax; jymax; vxmax; axmax; jxmax];
    numParam        = length(param);
    % void __cdecl Control (double *position, double *velocity, 
    %              double *action, int numAxis, double *wayPointX, 
    %              double *wayPointY, int numWaypoints, double *actualWayPoint, 
    %              double *param, int numParam, double *seconds)
    % Controller (Call DLL)
    p_position       = libpointer('doublePtr',position);
    p_velocity       = libpointer('doublePtr',velocity);
    p_action         = libpointer('doublePtr',zeros(2,1));
    p_wayPointX      = libpointer('doublePtr',wayPointX);
    p_wayPointY      = libpointer('doublePtr',wayPointY);
    p_actualWayPoint = libpointer('doublePtr',actualWayPoint);
    p_param          = libpointer('doublePtr',param);
    tic;
    calllib( 'MYDLL', funcon, p_position, p_velocity, p_action, numAxis, ...
                              p_wayPointX, p_wayPointY, numWaypoints, ...
                              p_actualWayPoint, p_param, numParam);
    % Measure Computational Time
    meast_c  = meast_c + toc;
    u_c(:,t) = get(p_action,'Value');
    
    clc
    disp(['Unconstrained MPC Progress: ',num2str(100*(t/largesize)),' %']);
    
    if res+2 <= length(ruta)
        if norm([yp_c(1,t),yp_c(3,t)]-ruta(res+2,:)) < (radio)
            res = res + 1;
        end
    else       
        if extratime == extrat
            break;
        end
        extratime = extratime + 1;
    end
    
    t = t + 1;
    
end

trobc = meast_c/(Ts*(t-1));

% Unload DLL
unloadlibrary('MYDLL');

% Cut Storage size
u_c    = u_c(:,1:t-1); % Input
xp_c   = xp_c(:,1:t-1); % Plant State
yp_c   = yp_c(:,1:t-1); % Plant Output
time   = time(:,1:t-1);

figure();
subplot(6,1,1),plot(time,yp_c(1,:),'r--');
axis([time(1) time(end) 0.9*min(yp_c(1,:)) 1.1*max(yp_c(1,:))]);

subplot(6,1,2),plot(time,yp_c(2,:),'r--');
axis([time(1) time(end) 0.9*min(yp_c(2,:)) 1.1*max(yp_c(2,:))]);

subplot(6,1,3),plot(time,yp_c(3,:),'r--');
axis([time(1) time(end) 0.9*min(yp_c(3,:)) 1.1*max(yp_c(3,:))]);

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