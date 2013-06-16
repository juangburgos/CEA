function [u,xest,xdist,yest,rs,lastout] = fcn(ymeas,ukm1,xpkm1,xdkm1,Ap,Bp,Cp, ...
                               Aj,Bj,Cj,Dj,rollparams,pitchparams, ...
                               Cpos,phi,jota,theta,gamma,omega,psi, ...
                               ruta, time, rskm1,lastin,radio)
                           
%% General Parameters

% Sampling Time
Ts = 0.06;
% Prediction Horizon
N = 50;
                           
%% Reference Designer

% % Third Order Setpoint
% Init Parameters
rutaini  = [0,0];
ruta     = [rutaini;ruta];
% Algorithm
dx       = ruta(rskm1+2,1)-ruta(rskm1+1,1);
dy       = ruta(rskm1+2,2)-ruta(rskm1+1,2);
vymax    = 1.5;                     % max 2.1 [m/s] ~1.5 (dont move)
aymax    = 0.55;                    % ~0.55 (dont move)
jymax    = 6;                       % ~6 (dont move)
vxmax    = 3.1;                     % max 3.5 [m/s] ~3.1
axmax    = 1;                       % ~1 (dont move)
jxmax    = 12;                      % ~12 (dont move)
% if dy ~= 0
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
    
tr_x = lastin + t;
tr_y = lastin + t;
pr_x = ruta(rskm1+1,1) + px;
pr_y = ruta(rskm1+1,2) + py;
% Update Reference State and Last Time
if norm([ymeas(1),ymeas(3)]-ruta(rskm1+2,:)) < (radio/1)
    if rskm1+1 < size(ruta,1)
        rs      = rskm1+1;
        lastout = time;
    else
        rs      = rskm1;
        lastout = lastin;
    end
else
    rs      = rskm1;
    lastout = lastin;
end

%% Gain Scheduling

unew1 = sinfit(ukm1(1),rollparams);
unew2 = sinfit(ukm1(2),pitchparams);
% Scheduled B matrix
Bsch     = Bp*[unew1 , 0; ...
               0, unew2];
           
%% Observer

% 1) PREDICT   
% Nominal Model State Equation
xest  = Ap * xpkm1 + Cj * xdkm1 + Bsch * ukm1;
% Disturbance Model State Equation
xdist = Aj * xdkm1;
% Output Equation
yest  = Cp * xest; 
% 2) CORRECT
% Update Nominal Model States
xest    = xest  + Dj*(ymeas-yest);
% Update Nominal Model States
xdist   = xdist + Bj*(ymeas-yest);
% Update Output
yest    = Cp * xest;

%% State Feedback

% 1) Compute Reference Vector
rtmp = zeros(size(Cpos,1),N);
mytime = zeros(1,N);
for i = 1:1:N
    mytime(i) = time+(i-1)*Ts;
end
rtmp(1,:) = interp1(tr_x, pr_x, mytime(1:N),'linear','extrap'); % X Position
rtmp(2,:) = interp1(tr_y, pr_y, mytime(1:N),'linear','extrap'); % Y Position
ref = zeros(N*size(Cpos,1),1);
for i = 1:1:N
    ref((i-1)*size(Cpos,1)+1:i*size(Cpos,1),1) = [rtmp(1,i);rtmp(2,i)];
end
% 2) Compute Unconstrained MPC Input
filas    = size(Bsch,1);
columnas = size(Bsch,2);
be = zeros(N*filas,N*columnas);
for k1=1:N
    be((k1-1)*filas+1:(k1)*filas, (k1-1)*columnas+1:(k1)*columnas ) = Bsch;
end
theta = theta*Bsch;
gamma = gamma*be;
test1 = gamma'*omega*gamma;
test2 = theta(1:2*N,1:2)*ukm1(1:2,1)-ref(1:2*N,1);
G   = 2*(psi+test1(1:2*N,1:2*N));
F   = 2*gamma'*omega*(phi*xest+jota*xdist+test2);
sol = -G\F;
du  = sol(1:size(Bsch,2),1);
u   = ukm1 + du;
% 3) Limit the Inputs
u(1) = min(u(1),1);
u(1) = max(u(1),-1);
u(2) = min(u(2),1);
u(2) = max(u(2),-1);


%% Third Order Setpoint Generator

function [tx,xp]=thirdord(p,v,a,j,Ts)
% function [tx,xp,xv,xa,xj]=mycode(p,v,a,j,Ts)

%% Part 1
% keyboard

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

%% Part 2

tt=   [0 1 1 2 2 3 3 4 ]*t(1) ...
    + [0 0 1 1 1 1 2 2 ]*t(2) ...
    + [0 0 0 0 1 1 1 1 ]*t(3) ;

ttest=[tt 1.5*tt(8)];
len = round(1.2*tt(8)/Ts + 1); % length of profiles
% xj = zeros(len,1);
xj = zeros(12*60,1);
xj = xj(1:len);
xa = xj;
xv = xj;
xp = xj;
xj(1) = jd;
% tx=0:Ts:1.2*tt(8)+Ts/2;
tx=0:Ts:12*60/Ts;
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


function out = sinfit(x,coeff)
    out=0;
    for i=1:3:max(size(coeff))
        out=out+coeff(i)*sin(coeff(i+1)*x+coeff(i+2));
    end
