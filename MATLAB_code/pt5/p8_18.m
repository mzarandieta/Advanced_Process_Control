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
Dw_new=Dw*Dwf;

% sample and hold 
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

%%%Part(i_8.18) I feel like we did that already 
Ad_new
Bd_new
Gd_new
Sigw


%%%Part(ii_8.18) compare the trajectories of the actual state variables -->using hw7 p8.6   +  example 8.3 book
%%%obtained from the deterministic and stochastic LQOC over 10h
I = eye(6);            %Ad_new*inv(Ad_new);                  %I think use Ad_new here instead of An? Yes
Q = 1000*I;                              %I = A*inv(A)
R = 0.001;
x0 = [0 0.5 -0.5 -30 -10 0]';             

[P, Lam, L]=dare(Ad_new,Bd_new,Q,R) 
Sigx=dlyap(Ad_new-Bd_new*L, Gd_new*Sigw*Gd_new')
Sigz=(Dx_new-Du_new*L)*Sigx*(Dx_new-Du_new*L)'

randn('state', 2^6-1); 
T=10; 
NN=T/dt;
nz=6;                                               %added this
wbar=0;                                             %added this
% tt=zeros(1,NN); 
% xx=zeros(nz,NN); 
% zz=zeros(nz,NN); 
% xx(:,1)=x0;                                         %x0 has to be modified from 5x1
% for kk=1:NN-1 
%     tt(kk+1)=dt*kk; 
%     wk=sqrt(Sigw)*randn+wbar; 
%     xx(:,kk+1)=(Ad_new-Bd_new*L)*xx(:,kk)+Gd_new*wk; 
% %     zz(:,kk)=(Dx_new-Du_new*L)*xx(:,kk); 
% end
%%%Part(iii) controler from (ii) over 10h 
s_ss = [362.22; 0.0154; 0.9846; 449.79; 385.57; 0];  % s_ss = [T_1 C_A C_B T_3 T_4]
% % Implement controller L (hw7 p.7.9)
X=zeros(nz,NN);
S=zeros(nz,NN);
X(:,1)=x0;
t=zeros(1,NN);
%u=zeros(1, 2);
dt=0.1;
for k=1:NN-1 
    t(k+1)=k*dt; 
         % for i=1:nx
%     U(:,k)=-L(k,:)*x(:,k); 
         %x(:,k+1)=(abs(Ad)-abs(Bd)*Liid)*x(:,k); 
    X(:,k+1)=(Ad_new-Bd_new*L)*X(:,k); 
end
S=X+s_ss;
titles=["T_1";"C_A";"C_B";"T_3";"T_4";"L"];
%%% Plots 
%temperature plots 
figure()
plot(t, S(1,:), 'r', t, S(4,:), 'b', t, S(5, :), 'g')
legend('T{1}','T{3}', 'T{4}') 
title('Temperatures') 
xlabel('Time [hours]','FontSize',10,'FontName','Times New Roman');
%PLOT CONCENTRATIONS 
figure()
plot(t, S(2,:), 'r', t, S(3,:), 'b')
legend('C{A}','C{B}') 
title('Concentrations') 
xlabel('Time [hours]','FontSize',10,'FontName','Times New Roman');
figure()
plot(t, S(6,:), 'r')
legend('No idea what this is') 
title('No idea') 
xlabel('Time [hours]','FontSize',10,'FontName','Times New Roman');


%%%Part(iii_8.18) plot T3 versus Q over 400h; scatter plot and EDOR of stochastic LQOC
Sigx=dlyap(Ad_new-Bd_new*L, Gd_new*Sigw*Gd_new')
Sigz=(Dx_new-Du_new*L)*Sigx*(Dx_new-Du_new*L)'




