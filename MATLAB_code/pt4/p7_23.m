clc; close all; clear all 

%%%%%%constants (no disturbance)
Q_ss = 2.845*10^6;       % heat (SSOP) [kcal/h]
s_ss = [362.2 0.0154 0.9846 449.8 385.6];
q_ss = [s_ss(4); Q_ss];  % assume no limitations on performance output q

%matrices
A = [-9.3750,0,0,0,6.875;
    0,-162.538579499525,0,-0.0122619049667084,0;
    0,160.038579499525,-2.5,0.0122619049667084,0;
    2.5,-32007.7158999051,0,-4.95238099334167,0;
    6.875,0,0,2.5,-9.375];

B = [0;0;0;0.00025;0];

G = [0;2.5;0;0;0];

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

%%%Part(i) Continous-time LQOC (pg254 -->ex 7.12) -->care
I = A*inv(A);
Q = 1000*I;                      %I = A*inv(A)
R = 0.001;
dt=0.1;
display('Part(i)') 
[P_i,Lam_i,L_i]=care(A,B,Q,R)


%%%Part(ii) discrete-time LQOC; use sample and hold method for system -->dare
%%%Report the values of Ad , Bd, Gd
Ad=expm(A*dt)
sum=zeros(5); 
Ndt=2000; 
ddt=dt/Ndt;
for ii=1:Ndt
    sum=sum+expm(A*ii*ddt); 
end
Bd=sum*B*ddt
Gd=sum*G*ddt

display('Part(ii)') 
[Pii,Lamii,Lii]=dare(Ad,Bd,Q,R)

display('Part(ii) with Q_d and R_d') 
[Piid,Lamiid,Liid]=dare(Ad,Bd,Q*dt,R*dt)


%%%Part(iii) controler from (ii) over 10h 
x0 = [0; 0.5; -0.5; -30; -10];
s_ss = [362.22; 0.0154; 0.9846; 449.79; 385.57];  % s_ss = [T_1 C_A C_B T_3 T_4]
% % Implement controller L (hw7 p.7.9)
nx=5;
x0=[0; 0.5; -0.5; -30; -10]; 
NN=100; 
X=zeros(nx,NN);
S=zeros(nx,NN);
X(:,1)=x0;
t=zeros(1,NN);
%u=zeros(1, 2);
dt=0.1;
for k=1:NN-1 
    t(k+1)=k*dt; 
         % for i=1:nx
%     U(:,k)=-L(k,:)*x(:,k); 
         %x(:,k+1)=(abs(Ad)-abs(Bd)*Liid)*x(:,k); 
    X(:,k+1)=(Ad-Bd*Liid)*X(:,k); 
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

% %%% Plots 
% %temperature plots 
% figure()
% plot(t, x(1,:), 'r', t, x(4,:), 'b', t, x(5, :), 'g')
% legend('T{1}','T{2}', 'T{3}') 
% title('Temperatures') 
% xlabel('Time [hours]','FontSize',10,'FontName','Times New Roman');
% %PLOT CONCENTRATIONS 
% figure()
% plot(t, x(2,:), 'r', t, x(3,:), 'b')
% % legend('C{A}','C{B}') 
% % title('Concentrations') 
% % xlabel('Time [hours]','FontSize',10,'FontName','Times New Roman');
% 
% %%%Comparing part (i) & (ii) Continuous tims vs discrete (dt=0.1) 
% %%% solution to deterministic LQOC. Plotting eigenvaleus 
%  figure()
%  x=[1, 2, 3, 4, 5, 6];
%  plot( x, Lami(:,1), 'r-o',x,  Lamii(:,1),'b-*',x, Lamiid(:,1),'g-o')
% legend('Continuous time LQOC', 'Discrete time LQOC', 'Discrete time LQOC with Qd and RD', 'Location','southwest')
% title('Comparaison between Eigenvalues')
% xlabel('\lambda','FontSize',10,'FontName','Times New Roman');
% 
%  figure()
%  x=[1, 2, 3, 4, 5, 6];
%  plot( x, Li(:,1), 'r-o',x,  Lii(:,1),'b-*',x, Liid(:,1),'g-o')
% legend('Continuous time LQOC', 'Discrete time LQOC', 'Discrete time LQOC with Qd and RD', 'Location','southwest')
% title('Comparaison between L')
% xlabel('\lambda','FontSize',10,'FontName','Times New Roman');

