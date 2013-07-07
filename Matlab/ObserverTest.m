clear all
clc

%% Load Data

load('experiment.mat');
load('FitsParaJUAN_k.mat');
load('model.mat');
load('distmod.mat');

ym = [experiment.x'; ...
      experiment.vx'; ...
      experiment.y'; ...
      experiment.vy'];
refroll  = experiment.refroll;
refpitch = experiment.refpitch;
time = experiment.time;
Ts = time(2)-time(1);


%% Robust Observer MATLAB

u  = zeros(size(dsys.b,2),max(size(time))); % Input
y  = zeros(size(dsys.c,1),max(size(time)));   % Output

% % Define a time to open the observer loop in order to assess
% % the quality of the prediction
% open = 1;
% if open
%     opentheloop = floor(max(size(time))/2); % Open at Half of Experiment
% else
%     opentheloop = max(size(time));          % Never Open
% end

tic
for t = 2:1:max(size(time))
    
    % Original Input
    u(:,t-1)   = [refroll(t-1); ...
                  refpitch(t-1)];
    % Observer
    y(:,t)     = MatObserver(u(:,t-1),ym(:,t));
    
    clc
    disp(['Progress: ',num2str(100*(t/max(size(time)))),' %']);
    
end
mat_t = toc/(Ts*(t-1));

% tline = [experiment.time(opentheloop),experiment.time(opentheloop)];

figure();
subplot(4,1,1),plot(experiment.time,ym(1,:),'b-');
hold on;
subplot(4,1,1),plot(experiment.time,y(1,:),'r--');
% hold on;
% subplot(4,1,1),plot(tline,[min(ym(1,:)),max(ym(1,:))],'g-','Linewidth',2);

subplot(4,1,2),plot(experiment.time,ym(2,:),'b-');
hold on;
subplot(4,1,2),plot(experiment.time,y(2,:),'r--');
% hold on;
% subplot(4,1,2),plot(tline,[min(ym(2,:)),max(ym(2,:))],'g-','Linewidth',2);

subplot(4,1,3),plot(experiment.time,ym(3,:),'b-');
hold on;
subplot(4,1,3),plot(experiment.time,y(3,:),'r--');
% hold on;
% subplot(4,1,3),plot(tline,[min(ym(3,:)),max(ym(3,:))],'g-','Linewidth',2);

subplot(4,1,4),plot(experiment.time,ym(4,:),'b-');
hold on;
subplot(4,1,4),plot(experiment.time,y(4,:),'r--');
% hold on;
% subplot(4,1,4),plot(tline,[min(ym(4,:)),max(ym(4,:))],'g-','Linewidth',2);

clc

%% Robust Observer C

% Load DLL
dllname    = 'libCEA.dll';
headername = 'mydll.h';
funcobs    = 'matobserver';
if ~libisloaded( 'MYDLL' ) 
    loadlibrary( dllname, headername, 'alias', 'MYDLL' );      
end

u_c  = zeros(size(dsys.b,2),max(size(time))); % Input
y_c  = zeros(size(dsys.c,1),max(size(time)));   % Output

% % Define a time to open the observer loop in order to assess
% % the quality of the prediction
% open = 1;
% if open
%     opentheloop = floor(max(size(time))/2); % Open at Half of Experiment
% else
%     opentheloop = max(size(time));          % Never Open
% end

tic
for t = 2:1:max(size(time))
    
    % Original Input
    u_c(:,t-1)   = [refroll(t-1); ...
                  refpitch(t-1)];
    % Observer (Call DLL)
    p_u       = libpointer('doublePtr',u_c(:,t-1));
    p_ym      = libpointer('doublePtr',ym(:,t));
    p_y       = libpointer('doublePtr',y_c(:,t));
    calllib( 'MYDLL', funcobs, p_u, p_ym, p_y );
    u_c(:,t-1)        = get(p_u,'Value');
    y_c(:,t)          = get(p_y,'Value');
    
    clc
    disp(['Progress: ',num2str(100*(t/max(size(time)))),' %']);
    
end
cprg_t = toc/(Ts*(t-1));

% Unload DLL
unloadlibrary('MYDLL');

% tline = [experiment.time(opentheloop),experiment.time(opentheloop)];

figure();
subplot(4,1,1),plot(experiment.time,ym(1,:),'b-');
hold on;
subplot(4,1,1),plot(experiment.time,y_c(1,:),'r--');
% hold on;
% subplot(4,1,1),plot(tline,[min(ym(1,:)),max(ym(1,:))],'g-','Linewidth',2);

subplot(4,1,2),plot(experiment.time,ym(2,:),'b-');
hold on;
subplot(4,1,2),plot(experiment.time,y_c(2,:),'r--');
% hold on;
% subplot(4,1,2),plot(tline,[min(ym(2,:)),max(ym(2,:))],'g-','Linewidth',2);

subplot(4,1,3),plot(experiment.time,ym(3,:),'b-');
hold on;
subplot(4,1,3),plot(experiment.time,y_c(3,:),'r--');
% hold on;
% subplot(4,1,3),plot(tline,[min(ym(3,:)),max(ym(3,:))],'g-','Linewidth',2);

subplot(4,1,4),plot(experiment.time,ym(4,:),'b-');
hold on;
subplot(4,1,4),plot(experiment.time,y_c(4,:),'r--');
% hold on;
% subplot(4,1,4),plot(tline,[min(ym(4,:)),max(ym(4,:))],'g-','Linewidth',2);

clc

disp(['Matlab Program: ',num2str(mat_t)]);
disp(['C Program: ',num2str(cprg_t)]);