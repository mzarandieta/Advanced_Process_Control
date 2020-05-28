%%Main function 
clc
close all 
clear all 

%%% A. CONTINUOUS TIME 

% Linearize system around SS values And calulate OPEN LOOP cov Matrix 
    [SigxOL,SigzOL, SigxCL, SigzCL, An, Bn, Gn, Dx_new, Dun, Dwn,L]= cont_time_EDOR();
    disp("CONTINUOUS TIME PROCESS:")
    disp("Open Loop with Filter New Compound Matrixes:")
    disp("An*x+Bn*u+Gn*w")
    An
    Bn
    Gn
disp("z=Dx_new*x+ Dun*u+ Dwn*w")
    Dx_new
    Dun
    Dwn

disp("Open Loop Covariance Matrixes:")
    SigxOL
    SigzOL

% Linearize system around SS values And calulate CLOSED LOOP cov Matrix 

disp("With controller u(t)=-L(x(t) and L=")
    L
disp("The values of the Covariance Matrixes in Closed Loop are:")
    SigxCL
    SigzCL

%Plot CLOSE LOOP continuous time EDOR 
disp("Standard deviation parameter alpha=2")

%Solve for steady state
pbar=1;

%%SIGNS 

sbarCC=inv(An)*Gn*pbar;
qbarCC=(Dx_new - Dun*L)*sbarCC;

%Define EDOR
zbar=qbarCC; 
sigZ=SigzCL;
alpha = 2;
Ezinv=inv(sigZ);

%Determie range of independent variable 

xind_max=sqrt(alpha^2*Ezinv(2,2)/det(Ezinv));
xind_min=-xind_max;
N=200; 
xind2=zeros(1,N);
yup=zeros(1,N);
ydown=zeros(1,N);
step = (xind_max-xind_min)/(N-1);

%calculate upper and lower curve values

for i=1:N
xind=xind_min+step*(i-1);
xind2(i)=xind;
b=Ezinv(1,2)*xind/Ezinv(2,2);
c=(Ezinv(1,1)*xind*xind-alpha^2)/Ezinv(2,2);
yup(i)=-b+sqrt(b^2-c);
ydown(i)=-b-sqrt(b^2-c);
end

%Plot upper and lower curves 

xind2=xind2+zbar(1); yup=abs(yup+zbar(2)); ydown=abs(ydown+zbar(2));
figure(1)
subplot(2,1,1);
hold off
plot(abs(zbar(1)),abs(zbar(2)),'*', xind2,yup,'r',xind2,ydown,'r','LineWidth',2)
xlabel('T_3')
ylabel('Q')
title('Continuous Time Ellipse')

%%% B. DISCRETE TIME 

% Linearize system around SS values And calulate OPEN LOOP cov Matrix 
    [SigxdOL,SigzdOL, SigxdCL, SigzdCL, Ad_new, Bd_new, Gd_new, Dx_new, Du_new, Dwn]= disc_time_EDOR;
    disp("Discrete TIME PROCESS:")
    disp("Open Loop with Filter New Compound Matrixes:")
    disp("Ad*x+Bd*u+Gd*w")
    Ad_new
    Bd_new
    Gd_new
disp("z=Dx_new*x+ Dun*u+ Dwn*w")
    Dx_new
    Dun
    Dwn

disp("Open Loop Covariance Matrixes:")
    SigxdOL
    SigzdOL

% Linearize system around SS values And calulate CLOSED LOOP cov Matrix 

disp("With controller u(t)=-L(x(t) and L=")
    L
disp("The values of the Covariance Matrixes in Closed Loop are:")
    SigxdCL
    SigzdCL

%Plot CLOSE LOOP continuous time EDOR 
disp("Standard deviation parameter alpha=2")

%Solve for steady state

sbarDD=-inv(eye(size(Ad_new))-Ad_new)*Gd_new*pbar; 
avdd=(Dx_new - Dun*L)*sbarDD;


%Define EDOR 
zbar=avdd; 
dt=0.1;
Sw=0.02;

% Change this too 

sigZd=SigzdCL;
Ezinvd=inv(sigZd);
%Determine ranges of independent variable 

xind_maxd=sqrt(alpha^2*Ezinvd(2,2)/det(Ezinvd)); xind_mind=-xind_maxd;
N=200; xind2d=zeros(1,N); yupd=zeros(1,N); yd=zeros(1,N);
step = (xind_maxd-xind_mind)/(N-1);

%Calculate upper and lower curve values 
for ii=1:N
xindd=xind_mind+step*(ii-1);
xind2d(ii)=xindd;
b=Ezinvd(1,2)*xindd/Ezinvd(2,2);
c=(Ezinvd(1,1)*xindd*xindd-alpha^2)/Ezinvd(2,2);
yupd(ii)=-b+sqrt(b^2-c);
yd(ii)=-b-sqrt(b^2-c);
end

%Plot upper and lower curves 
xind2d=xind2d+zbar(1); 
yupd=abs(yupd+zbar(2)); 
yd=abs(yd+zbar(2));
figure(1)
subplot(2,1,2);
hold on 
plot(abs(zbar(1)),abs(zbar(2)),'pk', xind2d,yupd,'b',xind2d,yd,'b','LineWidth',2)
xlabel('T_3')
ylabel('Q')
title('Discrete Time Ellipse')

% Simulate process
randn('state',2^6-1);
Ld = (10^4)*[0.0047 -1.6407 0 0.0084 0.0036 -0.7817];
NN=round(250/dt);
Sig_p = Sw/dt;
tt=zeros(1,NN); 
ss=zeros(6,NN); 
qq=zeros(2,NN); 
ss0=sbarDD; 
ss(:,1)=ss0;

for kk=1:NN-1
tt(kk+1)=dt*kk;
pk=sqrt(Sig_p)*randn+pbar;
ss(:,kk+1)=(Ad_new)*ss(:,kk)+Gd_new*pk;
%ss(:,kk+1)=(Ad_new-Bd_new*Ld)*ss(:,kk)+Gd_new*pk;
qq(:,kk)=(Dx_new - Dun*Ld)*ss(:,kk);
%qq(:,kk)=(Dx_new)*ss(:,kk);
end
qq(:,NN)=qq(:,NN-1);

figure()
hold off
plot( xind2,yup,'r',abs(zbar(1)),abs(zbar(2)),'*',xind2,ydown,'r','LineWidth',4)
% legend({'Continuous time'},'Location','southwest')
xlabel('T_3')
ylabel('Q')
hold on 
plot(abs(zbar(1)),abs(zbar(2)),'pk', xind2d,yupd,'b',xind2d,yd,'b','LineWidth',2)
hold on
plot(abs(qq(1,:)), abs(qq(2,:)),'k.');
title('Continuous Time and Discrete Time Ellipses')
% legend({'~','Continuous time','Discrete time', 'Simulated discrete time process'},'Location','southwest')

