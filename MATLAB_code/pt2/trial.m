

Dx=eye(size(A));
pbar=1; %Sw=0.1;

%not sure what that means --> Assume the disturbance,CA0, is characterized as being colored noise with a mean of 1 kmol e/m
st_div = 0.2;           % standard diviation [kmole/m^3]
tau = 0.25;             % correlation time [h]
Sw = 2*tau*st_div^2;    %not 100% sure about this eq but other project used it
dt=0.1;

%Solve for steady state
sbarCC=-inv(A)*G*pbar;
qbarCC=Dx*sbarCC
SigxOL = lyap(A, G*Sw*G')

SigzCC=Dx*SigxOL*Dx'

%Define EDOR
zbar=qbarCC; 
%%%Change for our values T3 (4th value) and Q (??)
sigZ=[SigzCC(1,1) SigzCC(1,3); SigzCC(3,1) SigzCC(3,3)];
alpha = 2;
Ezinv=inv(sigZ);
%Determie range of independent variable 

xind_max=sqrt(alpha^2*Ezinv(2,2)/det(Ezinv));
xind_min=-xind_max;
N=200; 
xind2=zeros(1,N);
yup=zeros(1,N);
y=zeros(1,N);
step = (xind_max-xind_min)/(N-1);

%calculate upper and lower curve values

for i=1:N
xind=xind_min+step*(i-1);
xind2(i)=xind;
b=Ezinv(1,2)*xind/Ezinv(2,2);
c=(Ezinv(1,1)*xind*xind-alpha^2)/Ezinv(2,2);
yup(i)=-b+sqrt(b^2-c);
y(i)=-b-sqrt(b^2-c);
end

%Plot upper and lower curves 

xind2=xind2+zbar(1); yup=yup+zbar(2); y=y+zbar(2);
figure()
hold off
plot(zbar(1),zbar(2),'*', xind2,yup,':',xind2,y,'--','LineWidth',4)
xlabel('Output 1')
ylabel('Output 3')

%Convert to Discrete time 
dt=0.1;
Ndt=200; ddt=dt/Ndt; sum=zeros(size(A));
for h=1:Ndt sum=sum+expm(A*h*ddt); end
Ad=expm(A*dt)
Gd=sum*G*ddt
Sig_p = Sw/dt;

%Solve for steady state
sbarDD=inv(eye(size(Ad))-Ad)*Gd*pbar; 
avdd=Dx*sbarDD

%Solve Covariance Matrix with dylap OPEN LOOP
Exdd=dlyap(Ad,Gd*Sig_p*Gd'); 
Ezdd=Dx*Exdd*Dx'

%Define EDOR 
zbar=avdd; 
% Change this too 
sigZ=[Ezdd(1,1) Ezdd(1,3); Ezdd(3,1) Ezdd(3,3)];
alpha = 2; Ezinv=inv(sigZ);
%Determine ranges of independent variable 
xind_max=sqrt(alpha^2*Ezinv(2,2)/det(Ezinv)); xind_min=-xind_max;
N=200; xind2=zeros(1,N); yup=zeros(1,N); y=zeros(1,N);
step = (xind_max-xind_min)/(N-1);

%Calculate upper and lowe curve values 
for ii=1:N
xind=xind_min+step*(ii-1);
xind2(ii)=xind;
b=Ezinv(1,2)*xind/Ezinv(2,2);
c=(Ezinv(1,1)*xind*xind-alpha^2)/Ezinv(2,2);
yup(ii)=-b+sqrt(b^2-c);
y(ii)=-b-sqrt(b^2-c);
end

%Plot upper and lower curves 
xind2=xind2+zbar(1); yup=yup+zbar(2); y=y+zbar(2);

figure()
hold on
plot(zbar(1),zbar(2),'pk', xind2,yup,':',xind2,y,'--','LineWidth',4)

%Simulate process
randn('state',2^6-1); NN=round(5000/dt);

tt=zeros(1,NN); ss=zeros(5,NN); qq=zeros(5,NN); ss0=sbarDD; ss(:,1)=ss0;
for kk=1:NN-1
tt(kk+1)=dt*kk;
pk=sqrt(Sig_p)*randn+pbar;
ss(:,kk+1)=Ad*ss(:,kk)+Gd*pk;
qq(:,kk)=Dx*ss(:,kk);
end

qq(:,NN)=qq(:,NN-1);
figure(1)
hold on
plot(qq(1,:),qq(3,:),'k.','MarkerSize',4);
title('dt=0.1')
