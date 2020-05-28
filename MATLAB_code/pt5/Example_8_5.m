clear all 
close all 
clc
% Mass spring damper model
    Ad=[0.941 0.192; -0.576 0.903]; 
    Bd=[.020; 0.192]; 
    Gd=Bd;
    Dx=[1 0;0 1;0 0]; 
    Du=[0;0;1]; 
    Q=[1 0;0 1]; 
    R=0.005; 
    C=[1 0];
    dt=0.2; 
    NN=200/dt; 
    Sigw=2.5; 
    Sigv=0.0025*eye(1);
% Find LQOC Gain
    P=dare(Ad,Bd,Q,R); 
    L=inv(R+Bd'*P*Bd)*Bd'*P*Ad;

% Find Optimal Estimation Gains
    Sige_plus=dare(Ad',C',Gd*Sigw*Gd',Sigv)
    K_plus=Ad*Sige_plus*C'*inv(C*Sige_plus*C'+Sigv);
    Sige=inv(inv(Sige_plus)+C'*inv(Sigv)*C); K=Sige*C'*inv(Sigv);
% Generate White Noise Sequences
    randn('state',2^6-1); 
    ww=sqrt(Sigw)*randn(NN,1); 
    vv=sqrt(Sigv)*randn(1,NN);
% Simulate with FSI
    tt=zeros(1,NN+1); xx=zeros(2,NN+1); uu=zeros(1,NN); xx(:,1)=[1;0];
   for kk=1:NN
        tt(kk)=(kk-1)*dt;
        uu(:,kk)=-L*xx(:,kk); xx(:,kk+1)=Ad*xx(:,kk)+Bd*uu(:,kk)+Gd*ww(kk);
    end
    tt(NN+1)=(NN)*dt; uu(:,NN+1)=-L*xx(:,NN+1);
    figure(1); plot(xx(2,:),xx(1,:),'k.');
    figure(2); plot(uu(1,:),xx(1,:),'k.')
% Simulate with PSI and without computational delay
    xx=zeros(2,NN+1); xhat=zeros(2,NN+1); xhat_plus=zeros(2,NN+1);
    uu=zeros(1,NN); xx(:,1)=[1;0]; xhat_plus(:,1)=[0;0];
    for kk=1:NN
        yykk=C*xx(:,kk)+vv(:,kk);
        xhat(:,kk)=xhat_plus(:,kk)+K*(yykk-C*xhat_plus(:,kk));
        uu(:,kk)=-L*xhat(:,kk); xx(:,kk+1)=Ad*xx(:,kk)+Bd*uu(:,kk)+Gd*ww(kk);
        xhat_plus(:,kk+1)=Ad*xhat(:,kk)+Bd*uu(:,kk);
    end


    tt(NN+1)=(NN)*dt; xhat(:,NN+1)=xhat_plus(:,NN+1);
    uu(:,NN+1)=-L*xhat(:,NN+1);\
    
    figure(11); 
    plot(xx(2,:),xx(1,:),'k.');
    figure(22); 
    plot(uu(1,:),xx(1,:),'k.')
    
% Simulate with PSI and with computational delay
    xx=zeros(2,NN+1);
    xhat=zeros(2,NN+1);
    xhat_plus=zeros(2,NN+1);
    uu=zeros(1,NN);
    xx(:,1)=[1;0];
    xhat_plus(:,1)=[0;0];
for kk=1:NN
    yy=C*xx(:,kk)+vv(:,kk);
    xhat(:,kk)=xhat_plus(:,kk)+K*(yy-C*xhat_plus(:,kk));
    uu(:,kk)=-L*xhat_plus(:,kk);
    xx(:,kk+1)=Ad*xx(:,kk)+Bd*uu(:,kk)+Gd*ww(kk);
    xhat_plus(:,kk+1)=Ad*xhat(:,kk)+Bd*uu(:,kk);
end
    xhat(:,NN+1)=xhat_plus(:,NN+1); u
    u(:,NN+1)=-L*xhat(:,NN+1);
    figure(111); 
    plot(xx(2,:),xx(1,:),'k.');
    figure(222); 
    plot(uu(1,:),xx(1,:),'k.')
    
% Find Output Covariance: FSI
    Sigx=dlyap(Ad-Bd*L,Gd*Sigw*Gd')
    Sigz_FSI=(Dx-Du*L)*Sigx*(Dx-Du*L)'
% Find Output Covariance: PSI without delay
    Sigxhat=dlyap(Ad-Bd*L,Sige_plus-Sige)
    Sigz_PSInD=(Dx-Du*L)*Sigxhat*(Dx-Du*L)'+Dx*Sige*Dx'
% Find Output Covariance: PSI with delay
    Sigxhat_plus=dlyap(Ad-Bd*L,Ad*(Sige_plus-Sige)*Ad')
    Sigz_PSIwD=(Dx-Du*L)*Sigxhat_plus*(Dx-Du*L)'+Dx*Sige_plus*Dx'