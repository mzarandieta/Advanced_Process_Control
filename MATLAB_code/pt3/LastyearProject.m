%Building HVAC System r for room, s for solid, SSfor steady state
clear all
Vr = 60;  %m^3
Vs = 10;  %m^3
Ar = 200; %m^2
As = 5;   %m^2
U = 5.69; %W/m^2 K
Ro = 40;  %mole/m^3
Cp = 29;  %J/mole K
RoCps = 2000; %J/m^3 K
Tc = 20;  %degree C
Cfre = 0;   %ppm fre=fresh
SSTr = 26; %degree C
SSTs = 26; %degree C
SSCr = 400; %ppm
SSFrcy = 0.453; %m^3/s rcy=recylce
SSFfre = 0.0375; %m^3/s fre=fresh
Tout = 29; %degree C  out=outside
SDTout = 3; %degree C standard deviation of Toutside 
Sc = 0.25; %ppm/s 
SDSc = 0.1; %ppm/s standard deviation of Sc
T = 24*3600;    %h, Simulation time
dt = 0.5*3600;  %h, sample time
xp0 = [34;30;390]-[26;26;400]; %initial value of physical process
x0 = [26;26;400;29;0.25]-[26;26;400;29;0.25]; %initial value of Compound System

Ap = zeros(3,3);
Ap(1,1) = -(SSFrcy+SSFfre)/Vr-(U*Ar/(Vr*Ro*Cp))-(U*As/(Vr*Ro*Cp));
Ap(2,1) = (U*As)/(Vs*RoCps);
Ap(1,2) = (U*As)/(Vr*Ro*Cp);
Ap(2,2) = -Ap(2,1);
Ap(3,3) = -(SSFfre/Vr);

Bp = zeros(3,2);
Bp(1,1) = (Tc-SSTr)/Vr;
Bp(1,2) = Bp(1,1);
Bp(3,2) = (Cfre-SSCr)/Vr;

Gp = zeros(3,2);
Gp(1,1) = (U*Ar)/(Vr*Ro*Cp);
Gp(3,2) = 1;

Dxp = zeros(4,3);
Dxp(1,1) = 1;
Dxp(2,3) = 1;

Dup = zeros(4,2);
Dup(3,1) = 1;
Dup(4,2) = 2;

Dwp = zeros(4,2);

C=[0 1 1];

tspan=[0 24*3600]; s0=xp0;
    
[t,s] = ode45('HVACEqs', tspan, s0); % Simulate the NonLinear Model

figure(1);
plot(t,s(:,1),'k-','LineWidth',4);
legend('Troom'); 
title('Non-linear Model','FontSize',14,'Fontname','Times New Roman');
ylabel('State Variables','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Seconds)','FontSize',14,'Fontname','Times New Roman');
figure(2);
plot(t,s(:,2),'k-','LineWidth',4);
legend('Tsolid'); 
title('Non-linear Model','FontSize',14,'Fontname','Times New Roman');
ylabel('State Variables','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Seconds)','FontSize',14,'Fontname','Times New Roman');
figure(3);
plot(t,s(:,3),'k-','LineWidth',4);
legend('Croom'); 
title('Non-linear Model','FontSize',14,'Fontname','Times New Roman');
ylabel('State Variables','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Seconds)','FontSize',14,'Fontname','Times New Roman');

% sample and hold
nxp=3; nup=2; Ndt=200; ddt=dt/Ndt; sum=zeros(nxp);
for jjj=1:Ndt
    sum=sum+expm(Ap*jjj*ddt);
end
Adp=expm(Ap*dt); Bdp=sum*Bp*ddt; Gdp=sum*Gp*ddt;

%{
%YALMIP
yalmip('clear');C1=[];C2=[];
%Define Variables
Z0=sdpvar(nxp);Z1=sdpvar(nup,nxp);
C1=[0<=Z0]; %Z0>0
C2=[[Z0-eye(nxp) (Adp*Z0-Bdp*Z1)
    (Adp*Z0-Bdp*Z1)' Z0]>=0];
Constraints=[C1,C2];Objective=[];
options=sdpsettings('verbose',0,'solver','mosek'); sol=optimize(Constraints,Objective,options);
if sol.problem==0 %problem is feasible
    Z0feas=value(Z0);Z1feas=value(Z1);LLfeas=Z1feas/Z0feas,eig(Adp-Bdp*LLfeas)
elseif sol.problem==1 %problem is infeasible
    display('Infesible Problem');
else 
    display('Hmm,Something went wrong!');
    sol.info;
end
 %}

% simulate the system using continuous-time model
NNN=T/dt; ttt=zeros(1,NNN); xxxCT=zeros(nxp,NNN); 
xxxCT(:,1)=xp0;
for kk=1:NNN-1
    ttt(kk+1)=dt*kk;
    xxxCT(:,kk)=expm(Ap*dt*(kk-1))*xp0;
end
xxxCT(:,kk+1)=expm(Ap*dt*(kk))*xp0;

figure(4);
plot(ttt,xxxCT(1,:),'k*-',ttt,xxxCT(2,:),'k--',ttt,xxxCT(3,:),'ko-');
legend('Troom','Tsolid','Croom'); 
title('Continuous-time Model','FontSize',14,'Fontname','Times New Roman');
ylabel('State Variables','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Seconds)','FontSize',14,'Fontname','Times New Roman');

%simulate the system using the discrete-time model
xxxDT=zeros(nxp,NNN); xxxDT(:,1)=xp0; 
for kk=1:NNN-1
    ttt(kk+1)=dt*kk;
    xxxDT(:,kk+1)=Adp*xxxDT(:,kk)
end

figure(5);
plot(ttt/3600,xxxDT(1,:),'k*-',ttt/3600,xxxDT(2,:),'k--',ttt/3600,xxxDT(3,:),'ko-');
legend('Troom','Tsolid','Croom'); 
title('Discrete-time Model','FontSize',14,'Fontname','Times New Roman');
ylabel('State Variables','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');

%define shaping filter
taoTout = 24*3600;
taoSc = 2*3600;
SwTout = 2*taoTout*SDTout^2;
SwSc = 2*taoSc*SDSc^2;
Af=[-1/taoTout 0;0 -1/taoSc]; Gf=[1/taoTout 0;0 1/taoSc]; 
Dxf = eye(2); Dwf = zeros(2); 
Swf=[SwTout 0;0 SwSc];

A = [Ap Gp*Dxf;0 0 0 -1/taoTout 0;0 0 0 0 -1/taoSc]; B = [Bp;0 0;0 0]; G=[Gp*Dwf;Gf];
Dx=[Dxp Dwp*Dxf]; Du=Dup; Dw=Dwp*Dwf;

% sample and hold
nx=5; Ndt=200; ddt=dt/Ndt; sum=zeros(nx); Sigw=Swf./dt;
for jjj=1:Ndt
    sum=sum+expm(A*jjj*ddt);
end
Ad=expm(A*dt); Bd=sum*B*ddt; Gd=sum*G*ddt;


%simulate the compound system for open loop
randn('state',2^6-1); xxxCS=zeros(nx,NNN); xxxCS(:,1)=x0; 
zzzOL=zeros(4,NNN);
for kk=1:NNN-1
    ttt(kk+1)=dt*kk;
    wk=sqrt(Sigw)*randn(2,1); 
    xxxCS(:,kk+1)=Ad*xxxCS(:,kk)+Gd*wk; 
    zzzOL(:,kk)=Dx*xxxCS(:,kk)+Dw*wk; 
end
zzzOL(:,NNN)=zzzOL(:,NNN-1);

figure(6);
plot(ttt/3600,xxxCS(1,:)+26,'k*-');
legend('Troom'); 
title('Compound system for Discrete-time Model,Open Loop','FontSize',14,'Fontname','Times New Roman');
ylabel('State Variables','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');
figure(7);
plot(ttt/3600,xxxCS(2,:)+26,'k*-');
legend('Tsolid'); 
title('Compound system for Discrete-time Model,Open Loop','FontSize',14,'Fontname','Times New Roman');
ylabel('State Variables','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');
figure(8);
plot(ttt/3600,xxxCS(3,:)+400,'ko-');
legend('Croom'); 
title('Compound system for Discrete-time Model,Open Loop','FontSize',14,'Fontname','Times New Roman');
ylabel('State Variables','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');
figure(9);
plot(ttt/3600,zzzOL(1,:)+26,'k*-');
legend('Troom'); 
title('Compound system for Discrete-time Model,Open-loop','FontSize',14,'Fontname','Times New Roman');
ylabel('Performance Output','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');
figure(10);
plot(ttt/3600,zzzOL(2,:)+400,'ko-');
legend('Croom'); 
title('Compound system for Discrete-time Model,Open-loop','FontSize',14,'Fontname','Times New Roman');
ylabel('Performance Output','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');
figure(11);
plot(ttt/3600,zzzOL(3,:)+0.453,'k*-');
legend('Frecycle'); 
title('Compound system for Discrete-time Model,Open-loop','FontSize',14,'Fontname','Times New Roman');
ylabel('Performance Output','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');
figure(12);
plot(ttt/3600,zzzOL(4,:)+0.0375,'ko-');
legend('Ffresh'); 
title('Compound system for Discrete-time Model,Open-loop','FontSize',14,'Fontname','Times New Roman');
ylabel('Performance Output','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');

%LQOC gain in discrete-time
%determin Pd and Ld using dare 
Qd=eye(nx).*dt; 
Rd=[1 0;0 1].*dt;

[Pd,Lamd,Ld]=dare(Ad,Bd,Qd,Rd);

% Output covariance in discrete-time with FSI (page 284)
Sigx=dlyap(Ad-Bd*Ld,Gd*Sigw*Gd');
SigzFSIDT=(Dx-Du*Ld)*Sigx*(Dx-Du*Ld)';

%closed-loop system
%u=[SSFrcy;SSFfre]; 
uuuCL=zeros(2,NNN); zzzCL=zeros(4,NNN); xxxCL=zeros(nx,NNN); xxxCL(:,1)=x0;
for kk=1:NNN-1
    ttt(kk+1)=dt*kk;
    wk=sqrt(Sigw)*randn(2,1); 
    xxxCL(:,kk+1)=(Ad-Bd*Ld)*xxxCL(:,kk)+Gd*wk; 
    zzzCL(:,kk)=(Dx-Du*Ld)*xxxCL(:,kk)+Dw*wk;  
    uuuCL(:,kk)=-Ld*xxxCL(:,kk);
    %vk=sqrt(Sigv)*randn(3,1);
    %yyy(:,kk)=C*xxx(:,kk)+vk
end
zzzCL(:,NNN)=zzzCL(:,NNN-1);
%design Optimal Filter
%Sigeplus=dare(Ad',C',Gd*Sigw*Gd',Sigv);
%Sige=inv(inv(Sigeplus)+C'*inv(Sigv)*C);

%implement Optimal Filter
%xxxhat=zeros(nx,NNN);xxxhatplus=zeros(nx,NNN);
%xxxhatplus(:,1)=zeros(nx,1);
%for kk=1:NNN-1
%xxxhat(:,kk)=xxxhatplus(:,kk)+K*(yyy(:,kk)-C*xxxhatplus(:,kk));
%xxxhatplus(:,kk+1)=Ad*xxxhat(:,kk);
%end
%xxxhat(:,NNN)=xxxhatplus(:,NNN)+K*(yyy(:,NNN)-C*xxxhatplus(:,NNN));

% Define EDOR
% solve for steady state
Sigz12=[SigzFSIDT(1,1) SigzFSIDT(1,3);SigzFSIDT(3,1) SigzFSIDT(3,3)];
qbarDT=[24;1];
alpha=2; Sigz12inv=inv(Sigz12);
%determine range of independent variable
xindmax=sqrt(alpha^2*Sigz12inv(2,2)/det(Sigz12inv));
xindmin=-xindmax;
N=200;  % how to define N?
xxxind=zeros(1,N);yyyup=zeros(1,N);yyylo=zeros(1,N);
step=(xindmax-xindmin)/(N-1);
% calculate upper and lower curve values
for iii=1:N
    xind=xindmin+step*(iii-1);
    xxxind(iii)=xind;
    bbb=Sigz12inv(1,2)*xind/Sigz12inv(2,2);
    ccc=(Sigz12inv(1,1)*xind*xind-alpha^2)/Sigz12inv(2,2);
    yyyup(iii)=-bbb+sqrt(bbb^2-ccc);
    yyylo(iii)=-bbb-sqrt(bbb^2-ccc);
end
xxxind=xxxind+qbarDT(1);yyyup=yyyup+qbarDT(2);yyylo=yyylo+qbarDT(2);


%Simulate Process
NN=(1000*3600)/dt;
tt=zeros(1,NN); xx=zeros(nx,NN); zz=zeros(4,NN); xx0=zeros(5,1); xx(:,1)=xx0;
for kk=1:NN-1
    tt(kk+1)=dt*kk
    wk=sqrt(Sigw)*randn(2,1); 
    xx(:,kk+1)=(Ad-Bd*Ld)*xx(:,kk)+Gd*wk;
    zz(:,kk+1)=(Dx-Du*Ld)*xx(:,kk)+Dw*wk;
end
zz(:,NN)=zz(:,NN-1)

figure(13);
plot(ttt/3600,xxxCL(1,:)+26,'k*-');
legend('Troom'); 
title('Compound system for Discrete-time Model,closed-Loop','FontSize',14,'Fontname','Times New Roman');
ylabel('State Variables','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');
figure(14);
plot(ttt/3600,xxxCL(2,:)+26,'k*-');
legend('Tsolid'); 
title('Compound system for Discrete-time Model,closed-Loop','FontSize',14,'Fontname','Times New Roman');
ylabel('State Variables','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');
figure(15);
plot(ttt/3600,xxxCL(3,:)+400,'ko-');
legend('Croom'); 
title('Compound system for Discrete-time Model,closed-Loop','FontSize',14,'Fontname','Times New Roman');
ylabel('State Variables','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');
%movegui(figure(6),'northwest')%
figure(16);
plot(ttt/3600,zzzCL(1,:)+26,'k*-',ttt/3600,26+2*(SigzFSIDT(1,1)^0.5),'-k.',ttt/3600,26-2*(SigzFSIDT(1,1)^0.5),'-k.');
legend('Troom','E[z_k]\pm2*\Sigma_z^{1/2}'); 
title('closed-loop','FontSize',14,'Fontname','Times New Roman');
ylabel('Performance Output','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');
figure(17);
plot(ttt/3600,zzzCL(2,:)+400,'ko-',ttt/3600,400+2*(SigzFSIDT(2,2)^0.5),'-k.',ttt/3600,400-2*(SigzFSIDT(2,2)^0.5),'-k.');
legend('Croom','E[z_k]\pm2*\Sigma_z^{1/2}'); 
title('closed-loop','FontSize',14,'Fontname','Times New Roman');
ylabel('Performance Output','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');
figure(18);
plot(ttt/3600,zzzCL(3,:)+0.453,'k*-',ttt/3600,0.453+2*(SigzFSIDT(3,3)^0.5),'-k.',ttt/3600,0.453-2*(SigzFSIDT(3,3)^0.5),'-k.');
legend('Frecycle','E[z_k]\pm2*\Sigma_z^{1/2}'); 
title('closed-loop','FontSize',14,'Fontname','Times New Roman');
ylabel('Performance Output','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');
figure(19);
plot(ttt/3600,zzzCL(4,:)+0.0375,'ko-',ttt/3600,0.0375+2*(SigzFSIDT(4,4)^0.5),'-k.',ttt/3600,0.0375-2*(SigzFSIDT(4,4)^0.5),'-k.');
legend('Ffresh','E[z_k]\pm2*\Sigma_z^{1/2}'); 
title('closed-loop','FontSize',14,'Fontname','Times New Roman');
ylabel('Performance Output','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');
figure(20);
plot(ttt/3600,uuuCL(1,:),'k*-');
legend('Frecycle','E[z_k]\pm2*\Sigma_z^{1/2}'); 
title('closed-loop','FontSize',14,'Fontname','Times New Roman');
ylabel('Manipulated Variables','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');
figure(21);
plot(ttt/3600,uuuCL(2,:),'ko-');
legend('Ffresh'); 
title('closed-loop','FontSize',14,'Fontname','Times New Roman');
ylabel('Manipulated Variables','FontSize',14,'Fontname','Times New Roman');
xlabel('Time(Hours)','FontSize',14,'Fontname','Times New Roman');
figure(22)
plot(qbarDT(1),qbarDT(2),'pk',xxxind,yyyup,'k--',xxxind,yyylo,'k--', 'LineWidth',4)
hold on
plot(zz(1,:)+24,zz(3,:)+1,'k.','MarkerSize',4)
ylabel('Croom','FontSize',14,'Fontname','Times New Roman');
xlabel('Troom','FontSize',14,'Fontname','Times New Roman');


