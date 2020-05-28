clc; close all; clear all 

%%%%%%constants
Q_ss = 2.845*10^6;       % heat (SSOP) [kcal/h]
C_A0_ss = 1;             % inlet conc. (SSOP) [kmole/m^3]
dt=0.1;                  % sample time [h]
s_ss = [362.2 0.0154 0.9846 449.8 385.6];
q_ss = [s_ss(4); Q_ss];  % assume no limitations on performance output q

pbar = 1;                % Assume the disturbance,CA0, is characterized as being colored noise with a mean of 1 kmol e/m
st_div = 0.2;            % standard diviation [kmole/m^3]
tau = 0.25;              % correlation time [h]
Sw = 2*tau*st_div^2;     % spectral density

Sv_1 = 40;               % zero mean white noise with spectral density 1
Sv_2 = 10;               % zero mean white noise with spectral density 2

%matrices
A = [-9.3750,0,0,0,6.875;
    0,-162.538579499525,0,-0.0122619049667084,0;
    0,160.038579499525,-2.5,0.0122619049667084,0;
    2.5,-32007.7158999051,0,-4.95238099334167,0;
    6.875,0,0,2.5,-9.375];

B = [0;0;0;0.00025;0];

G = [0;2.5;0;0;0];

%z=Dx*x+ Du*u+ Dw*w or  q=[T_3 ;Q]=Dx*[s]+ Du*[m]+ Dw*[0]  // where// z=q-qss / w=p-p22 disturbance inputs [Ca0]/ u=m-mss manipulated variables [Q]/ 
Dx = [0 0 0 1 0; 0 0 0 0 0];
Du = [0; 1];
Dw = [0; 0];

%shaping filter (used to model the disturbance)
st_div = 0.2;                   % standard diviation [kmole/m^3]
tau = 0.25;                     % correlation time [h]
Sw = 2*tau*st_div^2;            % spectral density

Swf = [Sw];                     % spectral density with filter
Af = [-1/tau]; 
Gf = [1/tau]; 
Dxf = eye(1);                   % Dx with filter
Dwf = zeros(1);                 % Dw with filter

%Continuous time 
An = [A G*Dxf;0 0 0 0 0 Af];    % new matrix A eq. 5.28
Bn = [B; 0];                    % new matrix B
Gn =[G*Dwf; Gf];                % new matrix G
Dx_new=[Dx Dw*Dxf];             % new matrix Dx
Du_new=Du; 
Dwn=Dw*Dwf;

% sample and hold --> (can be done either before or after addinf shaping filter) Eq 5.136  
nx=6; 
Ndt=200; 
ddt=dt/Ndt; 
sum=zeros(nx); 
Sigw=Swf./dt;
for jjj=1:Ndt
    sum=sum+expm(An*jjj*ddt);
end
Ad_new = expm(An*dt);                 % new matrix Ad
Bd_new = sum*Bn*ddt;                  % new matrix Bd
Gd_new = sum*Gn*ddt;                  % new matrix Gd

L = (10^4)*[0.0047 -1.6407 0 0.0084 0.0036 -0.7817];

%%% Part (i)discrete-time Kalman filter and simulate the closed-loop (10 hours with initial condition at steady-state)
display('i')
dt = 0.1;
Sv = diag([40 10]);
Sigv = Sv/dt;

% outout equation: y=Cx+v (in our case it is T1, T3 + white nosie) 
C = [1 0 0 0 0 0; 0 0 0 1 0 0];

% Simulate with Stochastic Disturbance
NNN=round(10/dt); 
ttt=zeros(1,NNN); 
xxx=zeros(nx,NNN);
xxx(:,1)=zeros(nx,1);
yyy=zeros(2,NNN); 
randn('state',2^6-1);

for ii=1:NNN-1
ttt(ii+1)=dt*ii;
ww=randn*sqrt(Sigw); 
xxx(:,ii+1)=Ad_new*xxx(:,ii)+Gd_new*ww;
vv=sqrt(Sigv)*randn(2,1); 
yyy(:,ii)=C*xxx(:,ii)+vv;
end

%%%% Design Optimal Filter
Sige_plus = dare(Ad_new',C',Gd_new*Sigw*Gd_new',Sigv)
Sige = inv(inv(Sige_plus)+C'*inv(Sigv)*C)
K = Sige*C'*inv(Sigv)

% Implement Optimal Filter
xxx_hat=zeros(nx,NNN); 
xxx_hat_plus=zeros(nx,NNN);
xxx_hat_plus(:,1)=zeros(nx,1);

for ii=1:NNN-1
xxx_hat(:,ii)=xxx_hat_plus(:,ii)+K*(yyy(:,ii)-C*xxx_hat_plus(:,ii));
xxx_hat_plus(:,ii+1)=Ad_new*xxx_hat(:,ii);
end
xxx_hat(:,NNN)=xxx_hat_plus(:,NNN)+K*(yyy(:,NNN)-C*xxx_hat_plus(:,NNN));

%Plot Estimates
figure(1)
subplot(1,2,1)
 plot(ttt,xxx(1,:),'r-',ttt,xxx_hat(1,:),'b--','linewidth',1)
legend('Actual T_1',' Estimate of T_1')
title('T_1 Comparaison')
xlabel('Time [h]','FontSize',10,'FontName','Times New Roman');
subplot(1,2,2)
plot(ttt,xxx(4,:),'r-',ttt,xxx_hat(4,:),'b--','linewidth',1)
legend('Actual T_3',' Estimate of T_3')
title('T_3 Comparaison')
xlabel('Time [h]','FontSize',10,'FontName','Times New Roman');

%%%Part (ii)one-step prediction form of Kalman filter (eq 6.76 or 6.82=ss)
display('ii')
Kplus = Ad_new*Sige_plus*C'*inv(C*Sige_plus*C'+Sigv) %% Sige_plus from Design Optimal Filter -->I used eq 6.76

% Implement Optimal Filter
xxx_hat=zeros(nx,NNN); 
xxx_hat_plus=zeros(nx,NNN);
xxx_hat_plus(:,1)=zeros(nx,1);

for ii=1:NNN-1
xxx_hat(:,ii)=xxx_hat_plus(:,ii)+Kplus*(yyy(:,ii)-C*xxx_hat_plus(:,ii));
xxx_hat_plus(:,ii+1)=Ad_new*xxx_hat(:,ii);
end
xxx_hat(:,NNN)=xxx_hat_plus(:,NNN)+Kplus*(yyy(:,NNN)-C*xxx_hat_plus(:,NNN));

%Plot Estimates
figure(2)
subplot(1,2,1) 
plot(ttt,xxx(1,:),'r-',ttt,xxx_hat(1,:),'b--','linewidth',1) %is this T1?
legend('Actual T_1',' Estimate of T_1')
title('One step Predictor T_1 Comparaison')
xlabel('Time [h]','FontSize',10,'FontName','Times New Roman');
subplot(1,2,2)
plot(ttt,xxx(4,:),'r-',ttt,xxx_hat(4,:),'b--','linewidth',1) %is this T3?
legend('Actual T_3',' Estimate of T_3')
title('One step Predictor T_3 Comparaison')
xlabel('Time [h]','FontSize',10,'FontName','Times New Roman');

%%%Part (iii)Make plots of the estimation errors and two standard deviation envelopes of error variances for all state variables for both the estimator and the one-step predicator.
%%%third plot on pg 214 textbook (e_k and +/- 2*SigmaE^1/2
%error_puls = x_k - x_hat_k_plus


%make two figures, one for estimator, other for step predictor 
error=zeros(nx,NNN);
meanerror=zeros(nx,NNN);
sdpos1=zeros(nx,NNN);
sdneg1=zeros(nx,NNN);
error_plus=zeros(nx,NNN);
meanerror_plus=zeros(nx,NNN);
sdpos2=zeros(nx,NNN);
sdneg2=zeros(nx,NNN);
ii=0;
for ii=1:NNN
error(:,ii)=xxx(:,ii)-xxx_hat(:,ii);
end 

% Stehart chart with upper limit and lower limit 
jjj=0;
    for jjj=1:nx
    meanerror(jjj,:)=mean(error(jjj,:));
    end

ii=0;
for ii=1:NNN
for jjj=1:nx
sdpos1(jjj,ii)=meanerror(jjj,ii)+2*Sige(jjj,jjj)^(0.5);
sdneg1(jjj,ii)=meanerror(jjj,ii)-2*Sige(jjj,jjj)^(0.5);
end
end 

titles=["T_1";"C_A";"C_B";"T_3";"T_4";"L"];

for j=1:nx
    figure()
    plot(ttt,error(j,:),'k-o' ,ttt,sdpos1(j,:),'r--',ttt,sdneg1(j,:),'r--','linewidth',1)
    legend('error','\pm 2 \Sigma^{1/2}')
    title(titles(j))
    xlabel('Time [h]','FontSize',10,'FontName','Times New Roman');
    
end 

for ii=1:NNN
error_plus(:,ii)=xxx(:,ii)-xxx_hat_plus(:,ii);
end 

% Stehart chart with upper limit and lower limit 
jjj=0;
    for jjj=1:nx
    meanerror_plus(jjj,:)=mean(error_plus(jjj,:));
    end

ii=0;
for ii=1:NNN
for jjj=1:nx
sdpos2(jjj,ii)=meanerror_plus(jjj,ii)+2*Sige_plus(jjj,jjj)^(0.5);
sdneg2(jjj,ii)=meanerror_plus(jjj,ii)-2*Sige_plus(jjj,jjj)^(0.5);
end
end 


titles=["T_1+";"C_A+";"C_B+";"T_3+";"T_4+";"L+"];

for j=1:nx
    figure()
    plot(ttt,error_plus(j,:),'k-o',ttt,sdpos2(j,:),'r--',ttt,sdneg2(j,:),'r--','linewidth',1)
    legend('error','\pm 2 \Sigma^{1/2}')
    title(titles(j))
    xlabel('Time [h]','FontSize',10,'FontName','Times New Roman');
    
end 
