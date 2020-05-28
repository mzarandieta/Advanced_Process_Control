function [SigxdOL,SigzdOL, SigxdCL, SigzdCL, Ad_new, Bd_new, Gd_new, Dx_new, Du_new, Dwn,L]= disc_time_EDOR()
%%%%%%constants
Q_ss = 2.845*10^6;       % heat (SSOP) [kcal/h]
C_A0_ss = 1;             % inlet conc. (SSOP) [kmole/m^3]
dt=0.1;
s_ss = [362.2 0.0154 0.9846 449.8 385.6];
q_ss = [s_ss(4); Q_ss];  % assume no limitations on performance output q
%m = [Q];  
%p = [C_A0];

%matrices
A = [-9.37500000000000,0,0,0,6.87500000000000;
    0,-162.538579499525,0,-0.0122619049667084,0;
    0,160.038579499525,-2.50000000000000,0.0122619049667084,0;
    2.50000000000000,-32007.7158999051,0,-4.95238099334167,0;
    6.87500000000000,0,0,2.50000000000000,-9.37500000000000];

B = [0;0;0;0.000250000000000000;0];

G = [0;2.50000000000000;0;0;0];

Dx = [0 0 0 1 0 ; 0 0 0 0 0 ];

Du=[0;1];

Dw = [0; 0];

%define shaping filter
st_div = 0.2;                   % standard diviation [kmole/m^3]
tau = 0.25;                     % correlation time [h]
Sw = 2*tau*st_div^2;            % spectral density

Swf = [Sw];                     % spectral density with filter
Af = [-1/tau];                  % A matrix with filter
Gf=[1/tau];                     % G matrix with filter
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

%%%(iv)open-loop (comapred to cont-time = super close yay)
SigxdOL=dlyap(Ad_new-Bd_new,Gd_new*Sigw*Gd_new');
SigzdOL=(Dx_new-Du_new)*SigxdOL*(Dx_new-Du_new)';

I = eye(6);                  %I think use Ad_new here instead of An? Yes
Q = 1000*I;                              %I = A*inv(A)
R = 0.001;
x0 = [0 0.5 -0.5 -30 -10 0]';    

%%%(ii) repeat part (i) with closed-loop --> u = -Lx(t)
    P=dare(Ad_new,Bd_new,Q,R); 
    L=inv(R+Bd_new'*P*Bd_new)*Bd_new'*P*Ad_new;
% SigxdCL=dlyap(Ad_new-Bd_new*L,Gd_new*Sigw*Gd_new');
% SigzdCL=(Dx_new-Du_new*L)*SigxdCL*(Dx_new-Du_new*L)';
SigxdCL=dlyap(Ad_new-Bd_new*L, Gd_new*Sigw*Gd_new')
SigzdCL=(Dx_new-Du_new*L)*SigxdCL*(Dx_new-Du_new*L)'

%%%(v) plot closed-loop discrete-time EDOR ellipse T3 and Q + use 'hold on' to compare cont. with disc. times plots



%%%(vi) simulate closed-loop discrete time; add appropriate outputs to EDOR + verify scatter plot corresponds to discrete-time ellipse

% %I took this from "help_past"; I fixed xx0 from (5,1) to (6,1) and NN bc we have in h and they had dt in sec
% NN=1000/dt;
% tt=zeros(1,NN); 
% xx=zeros(nx,NN); 
% zz=zeros(2,NN); 
% xx0=zeros(6,1); 
% xx(:,1)=xx0;
% for kk=1:NN-1
%     tt(kk+1)=dt*kk
%     wk=sqrt(Sigw)*randn(1,6); 
%     xx(:,kk+1)=(Ad_new-Bd_new*Ld)*xx(:,kk)+Gd_new*wk;
%     zz(:,kk+1)=(Dx_new-Du_new*Ld)*xx(:,kk)+Dwn*wk;
% end
% zz(:,NN)=zz(:,NN-1);
% 
% %not sure how helpful it, there are also "plots" this is at the very end
% of the help code "%Simulate Process"



