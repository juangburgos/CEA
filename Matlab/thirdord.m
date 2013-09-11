%% Third Order Setpoint Designer Function
function [tx,xp,xv,xa,xj]=thirdord(p,v,a,j,Ts)
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