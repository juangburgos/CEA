vxmax    = 3;                     % max 3.5 [m/s] ~3.1
axmax    = 1.5;                       % ~1 (dont move)
jxmax    = 1;                      % ~12 (dont move)
vymax    = 1.5;                     % max 2.1 [m/s] ~1.5 (dont move)
aymax    = 0.75;                    % ~0.55 (dont move)
jymax    = 0.5;                       % ~6 (dont move)

xchange = 0.5;
[time,xptraj,xvtraj,xatraj,xjtraj] = thirdord(xchange,vxmax,axmax,jxmax,0.06);

figure();
subplot(4,1,1),plot(time,xptraj,'b--','LineWidth',2);
subplot(4,1,2),plot(time,xvtraj,'b--','LineWidth',2);
subplot(4,1,3),plot(time,xatraj,'b--','LineWidth',2);
subplot(4,1,4),plot(time,xjtraj,'b--','LineWidth',2);

ychange = 0.5;
[time,yptraj,yvtraj,yatraj,yjtraj] = thirdord(ychange,vymax,aymax,jymax,0.06);

figure();
subplot(4,1,1),plot(time,yptraj,'b--','LineWidth',2);
subplot(4,1,2),plot(time,yvtraj,'b--','LineWidth',2);
subplot(4,1,3),plot(time,yatraj,'b--','LineWidth',2);
subplot(4,1,4),plot(time,yjtraj,'b--','LineWidth',2);