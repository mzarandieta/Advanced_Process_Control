function [X]= rxns_discrete(t,x)
%%%SSOP constants
Q_SS = 2.845e6;            % manipulated measurment --> Q_ss heat (SSOP) [kcal/h]
C_A0SS = 1;                % disturbance --> C_A0_ss inlet conc. (SSOP) [kmole/m^3]

Q = 2.845e6;               % manipulated measurment --> Q_ss heat (SSOP) [kcal/h]
C_A0 = 1;                  % disturbance --> C_A0_ss inlet conc. (SSOP) [kmole/m^3]

if ( t<0 )
    Q = 2.845e6;           % manipulated measurment --> Q_ss heat (SSOP) [kcal/h]
    C_A0 = 1;              % disturbance --> C_A0_ss inlet conc. (SSOP) [kmole/m^3]
elseif ( t>0 && t<10)
    Q = 2.7e6;             % manipulated measurment --> Q_ss heat (SSOP) [kcal/h]
    C_A0 = 1.2;            % disturbance --> C_A0_ss inlet conc. (SSOP) [kmole/m^3]
elseif ( t > 10 && t < 20)
    Q = 3.1e6;             % manipulated measurment --> Q_ss heat (SSOP) [kcal/h]
    C_A0 = 1.2;            % disturbance --> C_A0_ss inlet conc. (SSOP) [kmole/m^3]
end

Ad =[0.491325149210661,-7.90655159101652,0,0.0438761361021651,0.291511365549888;-1.07755691972936e-05,0.0117974188021300,0,-5.90284905365730e-05,-3.02893791379800e-06;1.07755691972713e-05,0.767003364269276,0.778800783071405,5.90284905365708e-05,3.02893791387334e-06;0.148330942772703,-154.084309104058,0,0.771030816653235,0.0438761361021508;0.307466324132495,-28.1278772314054,0,0.148330942772699,0.491325149210660];
Bd = [4.32082172648495e-07;-1.55694729737992e-09;1.55694729737991e-09;2.18354097417844e-05;2.18525181527072e-06];
Gd = [-0.729216424576047;0.0182031047226589;0.202982287542898;-40.6415861980369;-3.87057984837156];

nx=5;
N=20/dt; td=zeros(1,N); X=zeros(nx, N);
x0=zeros(nx,1); X(:,1)=x0;
for k=1:N-1
td(k+1)=dt*k; uk=[2.7e6-2.845e6]; wk=[0.2];
if (td(k)>10) uk=[3.1e6-2.845e6]; wk=[0.2];end
X(:,k+1)=Ad*X(:,k)+Bd*uk+Gd*wk;
end

X; 

end
