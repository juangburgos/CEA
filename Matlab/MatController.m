function [y,action] = MatController(position,velocity,numAxis,wayPointX, ...
            wayPointY,numWaypoints,actualWayPoint,param,numParam)
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
%                a %disposición del usuario para poder transmitir parámetros
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
persistent tr;
persistent pr_x;
persistent pr_y;

% DEFINES
% Sampling Time
persistent Ts;
if isempty(Ts)
    Ts = 0.06;
end
% Radio
persistent radio;
if isempty(radio)
    radio = 0.05;
end
% Prediction Horizon
persistent N;
if isempty(N)
    N = 50;
end

% Load Params
vymax    = param(1);                     % max 2.1 [m/s] ~1.5 (dont move)
aymax    = param(2);                    % ~0.55 (dont move)
jymax    = param(3);                       % ~6 (dont move)
vxmax    = param(4);                     % max 3.5 [m/s] ~3.1
axmax    = param(5);                       % ~1 (dont move)
jxmax    = param(6);                      % ~12 (dont move)

% Read Meas
ym(1,1) = position(1);
ym(2,1) = velocity(1);
ym(3,1) = position(2);
ym(4,1) = velocity(2);

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
    x     = mypinv(sysC)*ym;
    xd    = zeros(size(Aj,1),1);
    ukm1  = zeros(size(Bp,2),1);
    % <----------- START CONTROL ----------------------------------------
    % SUPERVISORY CONTROL
    % Reference Design
    dx       = wayPointX(res+2)-wayPointX(res+1);% In C it should be +1 and +0
    dy       = wayPointY(res+2)-wayPointY(res+1);% In C it should be +1 and +0
    if abs(dy) > abs(0.8*(vymax/vxmax)*dx)
        [tr,pr_y] = thirdord(dy,vymax,aymax,jymax,Ts);
        if dy < 0
            pr_y = -pr_y;
        end
        pr_x = (dx/dy)*pr_y;
    elseif dx == 0 && dy == 0
            tr  = [0,Ts];
            pr_x = [0,0];
            pr_y = pr_x;
    else
        [tr,pr_x] = thirdord(dx,vxmax,axmax,jxmax,Ts);
        if dx < 0
            pr_x = -pr_x;
        end
        pr_y = (dy/dx)*pr_x;
    end
    tr   = lastT + tr;
    pr_x = wayPointX(res+1) + pr_x;
    pr_y = wayPointY(res+1) + pr_y;
    
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
    mytime = zeros(N,1);
    for i = 1:1:N
        mytime(i) = iTime+(i-1)*Ts;
    end
    rtmpx = interpola(tr, pr_x, mytime); % X Position
    rtmpy = interpola(tr, pr_y, mytime); % Y Position
    %disp(['Inperpolated: ',mat2str(rtmpx)]);
    ref = zeros(N*size(sysCpos,1),1);
    for i = 1:1:N
        ref((i-1)*size(sysCpos,1)+1:i*size(sysCpos,1),1) = [rtmpx(i);rtmpy(i)];
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
    sol = myls(G,F);
    du  = sol(1:size(Bsch,2),1);
    ukm1   = ukm1 - du;
    % 3) Limit the Inputs
    ukm1(1) = min(ukm1(1),1);
    ukm1(1) = max(ukm1(1),-1);
    ukm1(2) = min(ukm1(2),1);
    ukm1(2) = max(ukm1(2),-1);
    % 4) Write Input
    action       = ukm1;
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
    if sqrt( (ym(1)-wayPointX(res+2))^2 + (ym(3)-wayPointY(res+2))^2 ) < radio
        if res+2 < numWaypoints
            % Update Reference State and Last Time
            res      = res+1;
            lastT    = iTime;
            % Reference Design
            dx       = wayPointX(res+2)-wayPointX(res+1);
            dy       = wayPointY(res+2)-wayPointY(res+1);
            if abs(dy) > abs(0.8*(vymax/vxmax)*dx)
                [tr,pr_y] = thirdord(dy,vymax,aymax,jymax,Ts);
                if dy < 0
                    pr_y = -pr_y;
                end
                pr_x = (dx/dy)*pr_y;
            elseif dx == 0 && dy == 0
                    tr  = [0,Ts];
                    pr_x = [0,0];
                    pr_y = pr_x;
            else
                [tr,pr_x] = thirdord(dx,vxmax,axmax,jxmax,Ts);
                if dx < 0
                    pr_x = -pr_x;
                end
                pr_y = (dy/dx)*pr_x;
            end
            tr   = lastT + tr;
            pr_x = wayPointX(res+1) + pr_x;
            pr_y = wayPointY(res+1) + pr_y;

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
    % Update Disturbance Model States
    xd      = xd + disB*(ym-y);
    % Update Output
    y       = sysC * x;
    % STATE FEEDBACK
    % 1) Compute Reference Vector
    mytime = zeros(N,1);
    for i = 1:1:N
        mytime(i) = iTime+(i-1)*Ts;
    end
    rtmpx = interpola(tr, pr_x, mytime); % X Position
    rtmpy = interpola(tr, pr_y, mytime); % Y Position
    ref = zeros(N*size(sysCpos,1),1);
    for i = 1:1:N
        ref((i-1)*size(sysCpos,1)+1:i*size(sysCpos,1),1) = [rtmpx(i);rtmpy(i)];
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
    sol = myls(G,F);
    du  = sol(1:size(Bsch,2),1);
    ukm1   = ukm1 - du;
    % 3) Limit the Inputs
    ukm1(1) = min(ukm1(1),1);
    ukm1(1) = max(ukm1(1),-1);
    ukm1(2) = min(ukm1(2),1);
    ukm1(2) = max(ukm1(2),-1);
    % 4) Write Input
    action       = ukm1;
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
% keyboard;
try
    % PART 1
    p=abs(p);
    v=abs(v);
    a=abs(a);
    j=abs(j);

    %disp('--------- t1 --------');
    % Calculation t1
    t1 = (p/(2*j))^(1/3) ; % largest t1 with bound on jerk
    %disp(num2str(t1));
    t1 = ceil(t1/Ts)*Ts; 
    %disp(num2str(t1));
    jd = 1/2*p/(t1^3); 
    % velocity test
    if v < jd*t1^2         % v bound violated ?
       t1 = (v/j)^(1/2) ;  % t1 with bound on velocity not violated
       %disp(num2str(t1));
       t1 = ceil(t1/Ts)*Ts; 
       %disp(num2str(t1));
       jd = v/(t1^2); 
    end
    % acceleration test
    if a < jd*t1     % a bound violated ?
       t1 = a/j ;    % t1 with bound on acceleration not violated
       %disp(num2str(t1));
       t1 = ceil(t1/Ts)*Ts; 
       %disp(num2str(t1));
       jd = a/t1; 
    end
    j = jd;  % as t1 is now fixed, jd is the new bound on jerk
    %disp(num2str(t1));
    %disp('--------- t2 --------');

    % Calculation t2
    t2 = ((t1^2)/4+p/j/t1)^(1/2) - (3/2)*t1 ;   % largest t2 with bound on acceleration
    %disp(num2str(t2));
    t2 = ceil(t2/Ts)*Ts;
    %disp(num2str(t2));
    jd = p/( 2*t1^3 + 3*t1^2*t2 + t1*t2^2 ); 

    % velocity test
    if v < (jd*t1^2 + jd*t1*t2)   % v bound violated ?
       t2 = v/(j*t1) - t1 ;       % t2 with bound on velocity not violated
       %disp(num2str(t2));
       t2 = ceil(t2/Ts)*Ts;
       %disp(num2str(t2));
       jd = v/( t1^2 + t1*t2 ); 
    end
    j = jd;  % as t2 is now fixed, jd is the new bound on jerk
    %disp('--------- t3 --------');
    % Calculation t3
    t3 = (p - 2*j*t1^3 - 3*j*t1^2*t2 - j*t1*t2^2)/v ; % t3 with bound on velocity
    %disp(num2str(t3));
    t3 = ceil(t3/Ts)*Ts; 
    %disp(num2str(t3));
    jd = p/( 2*t1^3 + 3*t1^2*t2 + t1*t2^2 + t1^2*t3 + t1*t2*t3 ); 

    % All time intervals are now calculated
    t=[ t1 t2 t3 ] ;
    %disp('t1,t2,t3 and jd ----');
    %disp([num2str(t1),' ',num2str(t2),' ',num2str(t3)]);
    %disp(num2str(jd));

    % PART 2
    
    tt = t*[0 1 1 2 2 3 3 4 ; ...
            0 0 1 1 1 1 2 2 ; ...
            0 0 0 0 1 1 1 1 ];

    ttest=[tt 1.5*tt(8)];
    %disp(['tt_end: ',num2str(tt(8))]);
    len = round(1.2*tt(8)/Ts + 1); % length of profiles (only malloc!)
    %disp(['len: ',num2str(len)]);
    
    xj = zeros(len,1);
    xa = xj;
    xv = xj;
    xp = xj;
    xj(1) = jd;
    tx=0:Ts:Ts*(len-1);

    for k=1:(len-1)  
      i = 0;
      for j=1:length(ttest)
          if ( Ts*(k + 1/2) <= ttest(j) )
              i = j - 1;
              break;
          end
      end     
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
    %disp(['xp: ',mat2str(xp)]);
    %disp('--------- END --------');

catch err
    keyboard;
end

%% Interpola Function for Interpolating
function my_pr = interpola(tr_x, pr_x, my_tr)
try
    my_pr = zeros(size(my_tr));
    tr_mem = 1;
    for i = 1:length(my_tr)
        my_flag = 0;
        for j = tr_mem:length(tr_x)     
            if my_tr(i) < tr_x(j)
                m = (pr_x(j)-pr_x(j-1))/(tr_x(j)-tr_x(j-1));
                b = pr_x(j) - m*tr_x(j);
                my_pr(i) = m*my_tr(i) + b;
                my_flag = 1;
                tr_mem = j;
                break;
            elseif my_tr(i) == tr_x(j)
                my_pr(i) = pr_x(j);
                my_flag = 1;
                tr_mem = j;
                break;
            else
                continue;
            end
            
        end       
        if my_flag == 1
            continue;
        else
            my_pr(i) = pr_x(end);
        end        
    end   
catch err
    keyboard;
end

%% Least Squares Optimization
function x = myls(R,F)
% Solves the equation R*x = F
% using the QR decomposition method
% R must be square of size (n,n) while
% F and x must be tall matrices of size (n,1)

try
% % QR Decomposition

% Initialize variables
n = size(R,1);
c = zeros(n,1);
d = zeros(n,1);

% Perform Householder QR decomposition
for k = 1:(n-1)
    scale = 0.0;
    for i = k:n
        scale = max(scale, abs(R(i,k)));
    end
    if (scale == 0.0)
        c(k) = 0.0;
        d(k) = 0.0;
    else
        for i = k:n
            R(i,k) = R(i,k) / scale;
        end
        sum = 0.0;
        for i = k:n
            sum = sum + R(i,k) * R(i,k);
        end
        sigma = sqrt(sum) * sign(R(k,k));
        R(k,k) = R(k,k) + sigma;
        c(k) = sigma * R(k,k);
        d(k) = -1.0 * scale * sigma;
        for j = (k+1):n
            sum = 0.0;
            for i = k:n
                sum = sum + R(i,k) * R(i,j);
            end
            tau = sum / c(k);
            for i = k:n
                R(i,j) = R(i,j) - tau * R(i,k);
            end
        end
    end
end
d(n) = R(n,n);

% Construct Q and erase temporary data in R
% TODO: Could inegrate with loop in next step
Q = eye(n);
for j = (n-1):-1:1
    Q(j:n,:) = Q(j:n,:) - ((1.0 / c(j)) * R(j:n,j)) * (R(j:n,j)' * Q(j:n,:));
    R((j+1):n,j) = 0;
    R(j,j) = d(j);
end

% % Compute Least Squares Solution

b = Q'*F;

x=zeros(n,1);
x(n)=b(n)/R(n,n);
for i = n-1:-1:1
    term = 0;
    for j = n:-1:i+1
        term = term + R(i,j)*x(j);
    end
    x(i) = (1/R(i,i))*(b(i)-term);
end

catch err
    keyboard;
end