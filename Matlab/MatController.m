function [y,u] = MatController(ym,ruta,radio)
% void __cdecl Control (
%
% double *position, Vector de 3 elementos que contiene la posición actual 
%                   del ARDrone: -position[0]: posición actual en X
%                                -position[1]: posición actual en Y
%                                ­position[2]: posición actual en Z
%
% double *velocity, Vector de 3 elementos que contiene la velocidad actual
%                   del ARDrone: -velocity[0]: velocidad actual en X
%                                -velocity[1]: velocidad actual en Y
%                                -velocity[2]: velocidad actual en Z
%
% double *action,   Vector de 2 elementos que contiene la acción de control
%                   del ARDrone: -action[0]: RefPitch
%                                -action[1]: RefRoll
%
% int numAxis,      Número de ejes, es una constante de variable igual a 3
%
% double *wayPointX,Vector de dimensión igual a numWaypoints que contiene
%                   la componente X de cada uno de los waypoints que forman
%                   la trayectoria (en orden).
%
% double *wayPointY,Vector de dimensión igual a numWaypoints que contiene 
%                   la componente Y de cada uno de los waypoints que
%                   forman la trayectoria (en orden).
%
% int numWaypoints, Número de WayPoints de la trayectoria. 
%
% double *actualWayPoint, Variable para comunicarle al ARDrone hacia que 
%                         Waypoint nos estamos intentando aproximar.
%
% double *param, Vector de dimensión igual a numParam que se encuentra
%                a disposición del usuario para poder transmitir parámetros
%                ente la aplicación TrackDrone Lite® y la dll en tiempo real.
%
% int numParam) Número de parámetros del vector param.

% Controller Status Variable (Init = 0, Run = 1)
persistent cStatus;
if isempty(cStatus)
    cStatus = 0;
end
% Old Control Input
persistent ukm1;
% Internal System Model State
persistent x;
% Internal System Model A Matrix
persistent sysA;
% Internal System Model B Matrix
persistent sysB;
% Internal System Model C Matrix
persistent sysC;
% Internal Disturbance Model State
persistent xd;
% Internal Disturbance Model A Matrix
persistent disA;
% Internal Disturbance Model B Matrix
persistent disB;
% Internal Disturbance Model C Matrix
persistent disC;
% Internal Disturbance Model D Matrix
persistent disD;
% Roll Scheduling Parameters
persistent sRoll;
% Roll Scheduling Parameters
persistent sPitch;
% Internal System Model C Position Matrix
persistent sysCpos;
% Prediction Matrix Phi
persistent pPhi;
% Prediction Matrix Jota
persistent pJota;
% Prediction Matrix Theta
persistent pTheta;
% Prediction Matrix Gamma
persistent pGamma;
% Prediction Matrix Omega
persistent pOmega;
% Prediction Matrix Psi
persistent pPsi;
% Internal Time
persistent iTime;
% Reference State
persistent res;
% Time when the last Waypoint was achieved
persistent lastT;
% References
persistent tr_x;
persistent tr_y;
persistent pr_x;
persistent pr_y;

% DEFINES
% Sampling Time
persistent Ts;
if isempty(Ts)
    Ts = 0.06;
end
% Prediction Horizon
persistent N;
if isempty(N)
    N = 50;
end

% If Observer Initialization
if cStatus == 0
    % Initlialize Parameters and States
    load('QU_Controller_Parameters.mat');
    sysA = Ap;
    sysB = Bp;
    sysC = Cp;
    disA = Aj;
    disB = Bj;
    disC = Cj;
    disD = Dj;
    sRoll   = rollparams;
    sPitch  = pitchparams;
    sysCpos = Cpos;
    pPhi    = phi;
    pJota   = jota;
    pTheta  = theta;
    pGamma  = gama;
    pOmega  = omega;
    pPsi    = pssi;
    iTime   = 0;
    res     = 0; % This index shound be 0 in C
    lastT   = 0;
    % Initial Conditions
    x     = sysC\ym;
    xd    = zeros(size(Aj,1),1);
    ukm1  = zeros(size(Bp,2),1);
    % <----------- START CONTROL ----------------------------------------
    % SUPERVISORY CONTROL
    % Reference Design
    dx       = ruta(res+2,1)-ruta(res+1,1); % In C it should be +1 and +0
    dy       = ruta(res+2,2)-ruta(res+1,2); % In C it should be +1 and +0
    vymax    = 1.5;                     % max 2.1 [m/s] ~1.5 (dont move)
    aymax    = 0.55;                    % ~0.55 (dont move)
    jymax    = 6;                       % ~6 (dont move)
    vxmax    = 3.1;                     % max 3.5 [m/s] ~3.1
    axmax    = 1;                       % ~1 (dont move)
    jxmax    = 12;                      % ~12 (dont move)
    if abs(dy) > abs(0.8*(vymax/vxmax)*dx)
        [t,py] = thirdord(dy,vymax,aymax,jymax,Ts);
        if dy < 0
            py = -py;
        end
        px = (dx/dy)*py;
    else
        [t,px] = thirdord(dx,vxmax,axmax,jxmax,Ts);
        if dx < 0
            px = -px;
        end
        py = (dy/dx)*px;
    end
    tr_x = lastT + t;
    tr_y = lastT + t;
    pr_x = ruta(res+1,1) + px;
    pr_y = ruta(res+1,2) + py;
% % %     % Check weather next point has been reached
% % %     if norm([ym(1),ym(3)]-ruta(res+2,:)) < (radio)
% % %         if res+1 < size(ruta,1)
% % %             % Update Reference State and Last Time
% % %             res      = res+1;
% % %             lastT    = iTime;
% % %         end
% % %     end
    % OBSERVER
    % 1) PREDICT
    unew1 = sinfit(ukm1(1),sRoll);
    unew2 = sinfit(ukm1(2),sPitch);
    % Scheduled B matrix
    Bsch     = sysB*[unew1 , 0; ...
                     0, unew2];    
    % Nominal Model State Equation
    x       = sysA * x + disC * xd + Bsch * ukm1;
    % Disturbance Model State Equation
    xd      = disA * xd;
    % Output Equation
    y       = sysC * x;
    % 2) CORRECT
    % Update Nominal Model States
    x       = x + disD*(ym-y);
    % Update Nominal Model States
    xd      = xd + disB*(ym-y);
    % Update Output
    y       = sysC * x;
    % STATE FEEDBACK
    % 1) Compute Reference Vector
    rtmp = zeros(size(sysCpos,1),N);
    mytime = zeros(1,N);
    for i = 1:1:N
        mytime(i) = iTime+(i-1)*Ts;
    end
    rtmp(1,:) = interp1(tr_x, pr_x, mytime(1:N),'linear','extrap'); % X Position
    rtmp(2,:) = interp1(tr_y, pr_y, mytime(1:N),'linear','extrap'); % Y Position
    ref = zeros(N*size(sysCpos,1),1);
    for i = 1:1:N
        ref((i-1)*size(sysCpos,1)+1:i*size(sysCpos,1),1) = [rtmp(1,i);rtmp(2,i)];
    end
    % 2) Compute Unconstrained MPC Input
    filas    = size(Bsch,1);
    columnas = size(Bsch,2);
    be = zeros(N*filas,N*columnas);
    for k1=1:N
        be((k1-1)*filas+1:(k1)*filas, (k1-1)*columnas+1:(k1)*columnas ) = Bsch;
    end
    ntheta = pTheta*Bsch;  % Warning, do not overwrite pTheta
    ngamma = pGamma*be;    % Warning, do not overwrite pGamma
    G   = 2*(pPsi+ngamma'*pOmega*ngamma);
    F   = 2*ngamma'*pOmega*(pPhi*x+pJota*xd+ntheta*ukm1-ref);
    sol = -G\F;
    du  = sol(1:size(Bsch,2),1);
    ukm1   = ukm1 + du;
    % 3) Limit the Inputs
    ukm1(1) = min(ukm1(1),1);
    ukm1(1) = max(ukm1(1),-1);
    ukm1(2) = min(ukm1(2),1);
    ukm1(2) = max(ukm1(2),-1);
    % 4) Write Input
    u       = ukm1;
    % Next internal time
    iTime   = iTime + Ts;
    % <----------- END CONTROL ----------------------------------------
    % Initialization Finished
    cStatus = 1;
% If Observer Running
elseif cStatus == 1
    % <----------- START CONTROL ----------------------------------------
    % SUPERVISORY CONTROL
    % Check weather next point has been reached
    if norm([ym(1),ym(3)]-ruta(res+2,:)) < (radio)
        if res+2 < size(ruta,1)
            % Update Reference State and Last Time
            res      = res+1;
            lastT    = iTime;
            % Reference Design
            dx       = ruta(res+2,1)-ruta(res+1,1); % In C it should be +1 and +0
            dy       = ruta(res+2,2)-ruta(res+1,2); % In C it should be +1 and +0
            vymax    = 1.5;                     % max 2.1 [m/s] ~1.5 (dont move)
            aymax    = 0.55;                    % ~0.55 (dont move)
            jymax    = 6;                       % ~6 (dont move)
            vxmax    = 3.1;                     % max 3.5 [m/s] ~3.1
            axmax    = 1;                       % ~1 (dont move)
            jxmax    = 12;                      % ~12 (dont move)
            if abs(dy) > abs(0.8*(vymax/vxmax)*dx)
                [t,py] = thirdord(dy,vymax,aymax,jymax,Ts);
                if dy < 0
                    py = -py;
                end
                px = (dx/dy)*py;
            else
                [t,px] = thirdord(dx,vxmax,axmax,jxmax,Ts);
                if dx < 0
                    px = -px;
                end
                py = (dy/dx)*px;
            end
            tr_x = lastT + t;
            tr_y = lastT + t;
            pr_x = ruta(res+1,1) + px;
            pr_y = ruta(res+1,2) + py;
        end
    end
    % OBSERVER
    % 1) PREDICT
    unew1 = sinfit(ukm1(1),sRoll);
    unew2 = sinfit(ukm1(2),sPitch);
    % Scheduled B matrix
    Bsch     = sysB*[unew1 , 0; ...
                     0, unew2];    
    % Nominal Model State Equation
    x       = sysA * x + disC * xd + Bsch * ukm1;
    % Disturbance Model State Equation
    xd      = disA * xd;
    % Output Equation
    y       = sysC * x;
    % 2) CORRECT
    % Update Nominal Model States
    x       = x + disD*(ym-y);
    % Update Nominal Model States
    xd      = xd + disB*(ym-y);
    % Update Output
    y       = sysC * x;
    % STATE FEEDBACK
    % 1) Compute Reference Vector
    rtmp = zeros(size(sysCpos,1),N);
    mytime = zeros(1,N);
    for i = 1:1:N
        mytime(i) = iTime+(i-1)*Ts;
    end
    rtmp(1,:) = interp1(tr_x, pr_x, mytime(1:N),'linear','extrap'); % X Position
    rtmp(2,:) = interp1(tr_y, pr_y, mytime(1:N),'linear','extrap'); % Y Position
    ref = zeros(N*size(sysCpos,1),1);
    for i = 1:1:N
        ref((i-1)*size(sysCpos,1)+1:i*size(sysCpos,1),1) = [rtmp(1,i);rtmp(2,i)];
    end
    % 2) Compute Unconstrained MPC Input
    filas    = size(Bsch,1);
    columnas = size(Bsch,2);
    be = zeros(N*filas,N*columnas);
    for k1=1:N
        be((k1-1)*filas+1:(k1)*filas, (k1-1)*columnas+1:(k1)*columnas ) = Bsch;
    end
    ntheta = pTheta*Bsch;  % Warning, do not overwrite pTheta
    ngamma = pGamma*be;    % Warning, do not overwrite pGamma
    test1 = ngamma'*pOmega*ngamma;
    test2 = ntheta(1:2*N,1:2)*ukm1(1:2,1)-ref(1:2*N,1);
    G   = 2*(pPsi+test1(1:2*N,1:2*N));
    F   = 2*ngamma'*pOmega*(pPhi*x+pJota*xd+test2);
    sol = -G\F;
    du  = sol(1:size(Bsch,2),1);
    ukm1   = ukm1 + du;
    % 3) Limit the Inputs
    ukm1(1) = min(ukm1(1),1);
    ukm1(1) = max(ukm1(1),-1);
    ukm1(2) = min(ukm1(2),1);
    ukm1(2) = max(ukm1(2),-1);
    % 4) Write Input
    u       = ukm1;
    % Next internal time
    iTime   = iTime + Ts;
    % <----------- END CONTROL ----------------------------------------
end



%% Sinfit Function for Scheduling
function out = sinfit(x,coeff)
out=0;
for i=1:3:max(size(coeff))
    out=out+coeff(i)*sin(coeff(i+1)*x+coeff(i+2));
end

%% Third Order Setpoint Designer Function
function [tx,xp]=thirdord(p,v,a,j,Ts)

% PART 1
p=abs(p);
v=abs(v);
a=abs(a);
j=abs(j);

tol = eps;   % tolerance required for continuous time calculations
jd = j;      % required for discrete time calculations

% Calculation t1
t1 = (p/(2*j))^(1/3) ; % largest t1 with bound on jerk
if Ts>0  
   t1 = ceil(t1/Ts)*Ts; 
   jd = 1/2*p/(t1^3); 
end
% velocity test
if v < jd*t1^2         % v bound violated ?
   t1 = (v/j)^(1/2) ;  % t1 with bound on velocity not violated
   if Ts>0  
       t1 = ceil(t1/Ts)*Ts; 
       jd = v/(t1^2); 
   end
end
% acceleration test
if a < jd*t1     % a bound violated ?
   t1 = a/j ;    % t1 with bound on acceleration not violated
   if Ts>0  
       t1 = ceil(t1/Ts)*Ts; 
       jd = a/t1; 
   end
end
j = jd;  % as t1 is now fixed, jd is the new bound on jerk

% Calculation t2
t2 = (t1^2/4+p/j/t1)^(1/2) - 3/2*t1 ;   % largest t2 with bound on acceleration
if Ts>0  
   t2 = ceil(t2/Ts)*Ts; 
   jd = p/( 2*t1^3 + 3*t1^2*t2 + t1*t2^2 ); 
end
if abs(t2)<tol, t2=0; end % for continuous time case
% velocity test
if v < (jd*t1^2 + jd*t1*t2)   % v bound violated ?
   t2 = v/(j*t1) - t1 ;       % t2 with bound on velocity not violated
   if Ts>0  
       t2 = ceil(t2/Ts)*Ts; 
       jd = v/( t1^2 + t1*t2 ); 
   end
end
if abs(t2)<tol, t2=0; end % for continuous time case
j = jd;  % as t2 is now fixed, jd is the new bound on jerk

% Calculation t3
t3 = (p - 2*j*t1^3 - 3*j*t1^2*t2 - j*t1*t2^2)/v ; % t3 with bound on velocity
if Ts>0  
   t3 = ceil(t3/Ts)*Ts; 
   jd = p/( 2*t1^3 + 3*t1^2*t2 + t1*t2^2 + t1^2*t3 + t1*t2*t3 ); 
end
if abs(t3)<tol, t3=0; end % for continuous time case

% All time intervals are now calculated
t=[ t1 t2 t3 ] ;

% PART 2

tt=   [0 1 1 2 2 3 3 4 ]*t(1) ...
    + [0 0 1 1 1 1 2 2 ]*t(2) ...
    + [0 0 0 0 1 1 1 1 ]*t(3) ;

ttest=[tt 1.5*tt(8)];
len = round(1.2*tt(8)/Ts + 1); % length of profiles
xj = zeros(len,1);
xa = xj;
xv = xj;
xp = xj;
xj(1) = jd;
tx=0:Ts:1.2*tt(8)+Ts/2;
tx = tx(1:len);
for time=Ts:Ts:(1.2*tt(8)+Ts/2)
  i = find( (time + Ts/2) <= ttest ); 
  i = i(1)-1;
  k = round(time/Ts);
  if i==1 || i==7 
      xj(k+1) =  jd;
  elseif i==3 || i==5 
      xj(k+1) = -jd;
  else
      xj(k+1) =  0;
  end
  xa(k+1) = xa(k) + xj(k)*Ts;
  xv(k+1) = xv(k) + xa(k)*Ts;
  xp(k+1) = xp(k) + xv(k)*Ts;
end