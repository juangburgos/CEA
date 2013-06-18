function y = MatObserver(u,ym)
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

% If Observer Initialization
if cStatus == 0
    % Initlialize Parameters and States
    load('QU_Controller_Parameters.mat');
    x    = zeros(size(Ap,1),1);
    sysA = Ap;
    sysB = Bp;
    sysC = Cp;
    xd   = zeros(size(Aj),1);
    disA = Aj;
    disB = Bj;
    disC = Cj;
    disD = Dj;
    sRoll   = rollparams;
    sPitch  = pitchparams;
    % Initial Conditions
    x     = sysC\ym;
    % 1) PREDICT
    unew1 = sinfit(u(1),sRoll);
    unew2 = sinfit(u(2),sPitch);
    % Scheduled B matrix
    Bsch     = sysB*[unew1 , 0; ...
                     0, unew2];    
    % Nominal Model State Equation
    x       = sysA * x + disC * xd + Bsch * u;
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
    % Initialization Finished
    cStatus = 1;
% If Observer Running
elseif cStatus == 1
    % 1) PREDICT
    unew1 = sinfit(u(1),sRoll);
    unew2 = sinfit(u(2),sPitch);
    % Scheduled B matrix
    Bsch     = sysB*[unew1 , 0; ...
                     0, unew2];    
    % Nominal Model State Equation
    x       = sysA * x + disC * xd + Bsch * u;
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
end

% Sinfit Function for Scheduling
function out = sinfit(x,coeff)
out=0;
for i=1:3:max(size(coeff))
    out=out+coeff(i)*sin(coeff(i+1)*x+coeff(i+2));
end