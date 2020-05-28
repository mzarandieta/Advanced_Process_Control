function [SigxOL,SigzOL, SigxCL, SigzCL, An, Bn, Gn, Dx_new, Dun, Dwn,L]= cont_time_EDOR()

%%%%%%constants
Q_ss = 2.845*10^6;       % heat (SSOP) [kcal/h]
C_A0_ss = 1;             % inlet conc. (SSOP) [kmole/m^3]
dt=0.1;
s_ss = [362.2 0.0154 0.9846 449.8 385.6];
q_ss = [s_ss(4); Q_ss];  % assume no limitations on performance output q
%m = [Q];  
%p = [C_A0];

pbar = 1;                % Assume the disturbance,CA0, is characterized as being colored noise with a mean of 1 kmol e/m
st_div = 0.2;            % standard diviation [kmole/m^3]
tau = 0.25;              % correlation time [h]
Sw = 2*tau*st_div^2;     % 

%matrices
A = [-9.3750,0,0,0,6.875;
    0,-162.538579499525,0,-0.0122619049667084,0;
    0,160.038579499525,-2.5,0.0122619049667084,0;
    2.5,-32007.7158999051,0,-4.95238099334167,0;
    6.875,0,0,2.5,-9.375];

B = [0;0;0;0.00025;0];

G = [0;2.5;0;0;0];


%z=Dx*x+ Du*u+ Dw*w or  q=[T_3 ;Q]=Dx*[s]+ Du*[m]+ Dw*[0]  // where// z=q-qss / w=p-p22 disturbance inputs [Ca0]/ u=m-mss manipulated variables [Q]/ 
Dx=[0 0 0 1 0; 0 0 0 0 0];
Du=[0; 1];
Dw = [0; 0];

%define shaping filter
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
Dun=Du; 
Dwn=Dw*Dwf;

% sample and hold --> for discrete time only
nx=6; Ndt=200; ddt=dt/Ndt; sum=zeros(nx); Sigw=Swf./dt;
for jjj=1:Ndt
    sum=sum+expm(An*jjj*ddt);
end
Ad_new = expm(An*dt);                 % new matrix Ad
Bd_new = sum*Bn*ddt;                  % new matrix Bd
Gd_new = sum*Gn*ddt;                  % new matrix Gd


%%%(i) linearizing around ss (open-loop) + ss covariance matrix (x) %u = 0;
SigxOL = lyap(An,Gn*Sw*Gn');
SigzOL = (Dx_new - Dun)*SigxOL*(Dx_new - Dun)'  ; 
I = eye(6);                  %I think use Ad_new here instead of An? Yes
Q = 1000*I;                              %I = A*inv(A)
R = 0.001;
x0 = [0 0.5 -0.5 -30 -10 0]';    

%%%(ii) repeat part (i) with closed-loop --> u = -Lx(t)
    P=dare(An,Bn,Q,R); 
    L=inv(R+Bn'*P*Bn)*Bn'*P*An;

SigxCL = lyap(An-Bn*L, Gn*Sw*Gn');

SigzCL = (Dx_new - Dun*L)*SigxCL*(Dx_new - Dun*L)' ;              

% %%%(iii) plot closed-loop cont.-time EDOR ellipse T3 and Q 
% % Define EDOR
% zbar=qbarCC; SigZ=[SigzCC(1,1) SigzCC(1,3); SigzCC(3,1) SigzCC(3,3)];
% alpha = 2; SigZinv=inv(SigZ);
% % Determine range of independent variable
% xind_max=sqrt(alpha^2*SigZinv(2,2)/det(SigZinv)); xind_min=-xind_max;
% N=200; xxxind=zeros(1,N); yyyup=zeros(1,N); yyylo=zeros(1,N);
% step = (xind_max-xind_min)/(N-1);
% % Calculate upper and lower curve values
% for iii=1:N
% xind=xind_min+step*(iii-1);
% xxxind(iii)=xind;
% bbb=SigZinv(1,2)*xind/SigZinv(2,2);
% ccc=(SigZinv(1,1)*xind*xind-alpha^2)/SigZinv(2,2);
% yyyup(iii)=-bbb+sqrt(bbb^2-ccc);
% yyylo(iii)=-bbb-sqrt(bbb^2-ccc);
% end
% % Plot upper and lower curves
% xxxind=xxxind+zbar(1); yyyup=yyyup+zbar(2); yyylo=yyylo+zbar(2);
% figure(100)
% hold off
% plot(zbar(1),zbar(2),'*k', xxxind,yyyup,'k--',xxxind,yyylo,'k--','LineWidth',4)
% xlabel('Output 1','FontSize',12)
% ylabel('Output 3','FontSize',12)