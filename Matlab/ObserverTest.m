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


%% Robust Observer

u  = zeros(size(dsys.b,2),max(size(time))); % Input
x  = zeros(size(dsys.a,1),max(size(time))); % Nominal Model States
xd = zeros(size(Jd.a,1),max(size(time))); % Disturbance Model States
y  = zeros(size(dsys.c,1),max(size(time)));   % Output

% x(:,1) = dsys.c\ym(:,1); % Initial Conditions
% y(:,1) = ym(:,1);        % Initial Output

% % Define a time to open the observer loop in order to assess
% % the quality of the prediction
% open = 1;
% if open
%     opentheloop = floor(max(size(time))/2); % Open at Half of Experiment
% else
%     opentheloop = max(size(time));          % Never Open
% end

for t = 2:1:max(size(time))
    
    % Original Input
    u(:,t-1)   = [refroll(t-1); ...
                  refpitch(t-1)];
    % Observer
    y(:,t)     = MatObserver(u(:,t-1),ym(:,t));
    
    
%     % 1) PREDICT
%     % Scheduled B matrix
%     Bsch = dsys.b*[feval(fitresultROLL_k,refroll(t-1)) , 0; ...
%                 0, feval(fitresultPITCH_k,refpitch(t-1))];     
%     % Nominal Model State Equation
%     x(:,t)  = dsys.a * x(:,t-1) + Jd.c * xd(:,t-1) + Bsch * u(:,t-1);
%     % Disturbance Model State Equation
%     xd(:,t) = Jd.a * xd(:,t-1);
%     % Output Equation
%     y(:,t)   = dsys.c * x(:,t);
%     
%     % Closed/Open Loop Prediction
%     if t <= opentheloop
%         % 2) CORRECT
%         % Update Nominal Model States
%         x(:,t)    = x(:,t)  + Jd.d*(ym(:,t)-y(:,t));
%         % Update Nominal Model States
%         xd(:,t)   = xd(:,t) + Jd.b*(ym(:,t)-y(:,t));
%         % Update Output
%         y(:,t)   = dsys.c * x(:,t);
%     end
    
    clc
    disp(['Progress: ',num2str(100*(t/max(size(time)))),' %']);
    
end

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