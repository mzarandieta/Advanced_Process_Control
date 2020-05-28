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

    C = [0 0 0 0 0 1; 1 0 0 0 0 0];   % differnet C now bc measurments are T0 and T1 ??


    %z=Dx*x+ Du*u+ Dw*w or  q=[T_3 ;Q]=Dx*[s]+ Du*[m]+ Dw*[0]  // where// z=q-qss / w=p-p22 disturbance inputs [Ca0]/ u=m-mss manipulated variables [Q]/
    Dx = [0 0 0 1 0; 0 0 0 0 0];
    %Dx = [0 0 0 0 1; 1 0 0 0 0];
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

    Sv = diag([40 10]);
    Sigv = Sv/dt;

    %%%Part(i_8.19) plot T3 versus Q over 200h; scatter plot and EDOR of stochastic PSI LQOC -->hw 7p8.11; book ex8.6
    % T0 = inlet temperature, which is the state of the shaping filter
    I = An*inv(An);                  %I think use Ad_new here instead of An? Yes
    Q = 1000*I;                      %I = A*inv(A)
    Qd = Q*dt;
    R = 0.001;
    Rd = R*dt;
    x0 = [0 0.5 -0.5 -30 -10 0]';

    % LQOC gain in continuous-time
    [P,Lam_CT,L]=care(An,Bn,Q,R);

    % Output covariance in continuous-time with FSI
    Sigx_CT=lyap(An-Bn*L,Gn*Sw*Gn');
    Sigz_FSI_CT=(Dx_new-Du_new*L)*Sigx_CT*(Dx_new-Du_new*L)';

    % Output covariance in continuous-time with PSI
    [Sige_CT,Lam_CT,K_CT]=care(An',C',Gn*Sw*Gn',Sv);
    K_CT=K_CT';
    Sigxhat_CT=lyap(An-Bn*L,Sige_CT*C'*inv(Sv)*C*Sige_CT);
    Sigz_PSI_CT=(Dx_new-Du_new*L)*Sigxhat_CT*(Dx_new-Du_new*L)'+Dx_new*Sige_CT*Dx_new';

    % LQOC gain in disctrete-time
    [Pd, Lamd, Ld]=dare(Ad_new,Bd_new,Qd,Rd);

    % Output covariance in discrete-time with FSI
    Sigx=dlyap(Ad_new-Bd_new*Ld, Gd_new*Sigw*Gd_new');
    Sigz=(Dx_new-Du_new*Ld)*Sigx*(Dx_new-Du_new*Ld)';

    % Optimal estmator gains
    [Sige_plus,Lam,K_plus]=dare(Ad_new',C',Gd_new*Sigw*Gd_new',Sigv);
    K_plus=K_plus';

    Sige=inv(inv(Sige_plus)+C'*inv(Sigv)*C);
    K=Sige*C'*inv(Sigv);

    % Output covariance in discrete-time with PSI and no delay
    Sigxhat=dlyap(Ad_new-Bd_new*Ld,Sige_plus-Sige);
    Sigz_PSI_noD=(Dx_new-Du_new*Ld)*Sigxhat*(Dx_new-Du_new*Ld)'+Dx_new*Sige*Dx_new'

    %%%Part(ii_8.19) plot T3 versus Q over 200h; scatter plots and EDORs of stochastic PSI LQOC with computational delay -->hw 7p8.11

    % Output covariance in discrete-time with PSI and with delay
    Sigxhat_plus=dlyap(Ad_new-Bd_new*Ld,Ad_new*(Sige_plus-Sige)*Ad_new');
    Sigz_PSI_withD=(Dx_new-Du_new*Ld)*Sigxhat_plus*(Dx_new-Du_new*Ld)'+Dx_new*Sige_plus*Dx_new'

    %Plot CLOSE LOOP continuous time EDOR
    disp("Standard deviation parameter alpha=2")
    %Solve for steady state

    I = eye(6);  %An*inv(An);                  %I think use Ad_new here instead of An? Yes
    Q = 1000*I;                              %I = A*inv(A)
    R = 0.001;
    x0 = [0 0.5 -0.5 -30 -10 0]';

    %%%(ii) repeat part (i) with closed-loop --> u = -Lx(t)
    P=dare(Ad_new,Bd_new,Q,R);
    L=inv(R+Bd_new'*P*Bd_new)*Bd_new'*P*Ad_new;

    sbarDD=-inv(eye(size(Ad_new))-Ad_new)*Gd_new*pbar;
    avdd=(Dx_new - Du_new*L)*sbarDD;

    %Define EDOR
    zbar=avdd;
    dt=0.1;
    Sw=0.02;
    
    %% NO DELAY PLOT
    % Change this too

    sigZd=Sigz_PSI_noD;
    Ezinvd=inv(sigZd);
    
    %Determine ranges of independent variable
    alpha=2;
    xind_maxd=sqrt(alpha^2*Ezinvd(2,2)/det(Ezinvd));
    xind_mind=-xind_maxd;
    N=200;
    xind2d=zeros(1,N);
    yupd=zeros(1,N);
    yd=zeros(1,N);
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
    %subplot(2,1,2);
    figure()
    hold on
    plot(abs(zbar(1)),abs(zbar(2)),'pk', xind2d,yupd,'b',xind2d,yd,'b','LineWidth',2)
    xlabel('T_3')
    ylabel('Q')
    title('Without Delay')

    %% WITH DELAY PLOT
    sigZd=Sigz_PSI_withD;
    Ezinvd=inv(sigZd);
    zbarWD=avdd;
    %Determine ranges of independent variable
    alpha=2;
    xind_maxdWD=sqrt(alpha^2*Ezinvd(2,2)/det(Ezinvd));
    xind_mindWD=-xind_maxdWD;
    N=200;
    xind2dWD=zeros(1,N);
    yupdWD=zeros(1,N);
    ydWD=zeros(1,N);
    step = (xind_maxdWD-xind_mindWD)/(N-1);

    %Calculate upper and lower curve values
    for ii=1:N
        xinddWD=xind_mindWD+step*(ii-1);
        xind2dWD(ii)=xinddWD;
        bWD=Ezinvd(1,2)*xinddWD/Ezinvd(2,2);
        cWD=(Ezinvd(1,1)*xinddWD*xinddWD-alpha^2)/Ezinvd(2,2);
        yupdWD(ii)=-bWD+sqrt(bWD^2-cWD);
        ydWD(ii)=-bWD-sqrt(bWD^2-cWD);
    end

    %Plot upper and lower curves
    xind2dWD=xind2dWD+zbarWD(1);
    yupdWD=abs(yupdWD+zbarWD(2));
    ydWD=abs(ydWD+zbarWD(2));
%      ydWD=(ydWD+zbarWD(2));
    %subplot(2,1,2);

    hold on
    plot(abs(zbarWD(1)),abs(zbarWD(2)),'pk', xind2dWD,yupdWD,'r',xind2dWD,ydWD,'r','LineWidth',2)
    xlabel('T_3')
    ylabel('Q')
    title('PSI LQOC')
    legend('With Delay', 'Without Delay')

    %% PROCESS SIMULATION

    % Simulate process
    randn('state',2^6-1);
    % Ld = (10^4)*[0.0047 -1.6407 0 0.0084 0.0036 -0.7817];
    NN=round(400/dt);
    Sig_p = Sw/dt;
    tt=zeros(1,NN);
    ss=zeros(6,NN);
    qq=zeros(2,NN);
    ss0=sbarDD;
    ss(:,1)=ss0;

    for kk=1:NN-1
    tt(kk+1)=dt*kk;
    pk=sqrt(Sig_p)*randn+pbar;
    %ss(:,kk+1)=(Ad_new)*ss(:,kk)+Gd_new*pk;
    ss(:,kk+1)=(Ad_new-Bd_new*L)*ss(:,kk)+Gd_new*pk;
    %qq(:,kk)=(Dx_new - Du_new*L)*ss(:,kk);
%     qq(:,kk)= (C-Sv*L)*ss(:,kk);
   qq(:,kk)=(Dx_new-Du_new*L)*ss(:,kk);
 end
    %
    qq(:,NN)=qq(:,NN-1);    
    Q_ss = 2.845*10^6;       % heat (SSOP) [kcal/h]
    s_ss = [362.2 0.0154 0.9846 449.8 385.6];
    qq0=[362.22;Q_ss]*ones(1,NN);
    qqf=qq+qq0;
    %% PLOT FOR COMPARING WITH AND WITHOUT DELAYS
    figure()
    hold off
    % plot( xind2,yup,'r',abs(zbar(1)),abs(zbar(2)),'*',xind2,ydown,'r','LineWidth',4)
    % legend({'Continuous time'},'Location','southwest')
    xlabel('T_3')
    xlim([220 420])
    ylabel('Q')
    ylim([3*10^4 6*10^4])
    hold on
    plot(abs(zbarWD(1)),abs(zbarWD(2)),'pk', xind2dWD,yupdWD,'r',xind2dWD,ydWD,'r','LineWidth',2)
    title('PSI LQOC with Delay')
    hold on
   % plot(abs(zbar(1)),abs(zbar(2)),'pk', xind2d,yupd,'b',xind2d,yd,'b','LineWidth',2)
figure()
    hold on
    plot(qqf(1,:), qqf(2,:),'k.');
    title('Continuous Time and Discrete Time Ellipses')
    legend('With Delay', 'Without Delay', 'process simulation');
    % legend({'~','Continuous time','Discrete time', 'Simulated discrete time process'},'Location','southwest')
    %% SIMULATIN PROCESS

    % Generate White Noise Sequences
    randn('state',2^6-1); ww=sqrt(Sigw)*randn(NN,1); vv=sqrt(Sigv)*randn(2,NN);

    %% Simulate with FSI
    tt=zeros(1,NN+1); xx=zeros(6,NN+1); uu=zeros(1,NN); xx(:,1)=x0;
    for kk=1:NN
        tt(kk)=(kk-1)*dt;
        uu(:,kk)=-L*xx(:,kk);
        xx(:,kk+1)=Ad_new*xx(:,kk)+Bd_new*uu(:,kk)+Gd_new*ww(kk);
        qq(:,kk)=(Dx_new-Du_new*L)*xx(:,kk);
    end
    tt(NN+1)=(NN)*dt; 
    uu(:,NN+1)=-L*xx(:,NN+1);
    qq(:,NN)=qq(:,NN-1);
    %figure(31); plot(xx(1,:),xx(6,:),'k.');
    %figure(32); plot(abs(uu(1,:)),abs(xx(1,:)),'k.')
%figure(32); plot(qq(1,:),qq(2,:),'k.')
    %% Simulate with PSI and without computational delay
    xx=zeros(6,NN+1); xhat=zeros(6,NN+1); xhat_plus=x0;
    %zeros(6,NN+1);
    uu=zeros(1,NN); xx(:,1)=x0; xhat_plus(:,1)=zeros(6,1);
    qq=zeros(2,NN);
    for kk=1:NN
        yykk=C*xx(:,kk)+vv(:,kk);
        xhat(:,kk)=xhat_plus(:,kk)+K*(yykk-C*xhat_plus(:,kk));
        uu(:,kk)=-L*xhat(:,kk); 
        xx(:,kk+1)=Ad_new*xx(:,kk)+Bd_new*uu(:,kk)+Gd_new*ww(kk);
        xhat_plus(:,kk+1)=Ad_new*xhat(:,kk)+Bd_new*uu(:,kk);
       % qq(:,kk)=(Dx_new-Du_new*L)*xhat_plus(:,kk);
        qq(:,kk)=(Dx_new-Du_new*L)*xx(:,kk);
    end
    tt(NN+1)=(NN)*dt; 
    xhat(:,NN+1)=xhat_plus(:,NN+1);
    uu(:,NN+1)=-L*xhat(:,NN+1);
    qq(:,NN)=qq(:,NN-1);
    %figure(11); plot(xx(1,:),xx(6,:),'k.');
    %figure(22); plot(uu(1,:),xx(1,:),'k.')
    %figure(22); plot(qq(1,:),qq(2,:),'k.')
    
    %% Simulate with PSI and with computational delay
    xx=zeros(6,NN+1); xhat=zeros(6,NN+1); xhat_plus=x0;
    %zeros(6,NN+1);
    uu=zeros(1,NN); xx(:,1)=x0; xhat_plus(:,1)=zeros(6,1);
    qq=zeros(2,NN);
    for kk=1:NN
        yy=C*xx(:,kk)+vv(:,kk);
        xhat(:,kk)=xhat_plus(:,kk)+K*(yy-C*xhat_plus(:,kk));
        uu(:,kk)=-L*xhat_plus(:,kk);
        xx(:,kk+1)=Ad_new*xx(:,kk)+Bd_new*uu(:,kk)+Gd_new*ww(kk);
        xhat_plus(:,kk+1)=Ad_new*xhat(:,kk)+Bd_new*uu(:,kk);
       % qq(:,kk)=(Dx_new-Du_new*L)*xhat_plus(:,kk);
         qq(:,kk)=(Dx_new-Du_new*L)*xx(:,kk);
    end
    xhat(:,NN+1)=xhat_plus(:,NN+1); 
    uu(:,NN+1)=-L*xhat(:,NN+1);

   % figure(111); plot(xx(1,:),xx(6,:),'k.');
    figure(223); plot(qqf(1,:),qqf(2,:),'k.')
    %figure(222); plot(uu(1,:),xx(1,:),'k.')



