function [dxdt]= rxns_ode_lin(t,x)
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
 
A = [-9.37500000000000,0,0,0,6.87500000000000;
    0,-162.538579499525,0,-0.0122619049667084,0;
    0,160.038579499525,-2.50000000000000,0.0122619049667084,0;
    2.50000000000000,-32007.7158999051,0,-4.95238099334167,0;
    6.87500000000000,0,0,2.50000000000000,-9.37500000000000];

B = [0;0;0;0.000250000000000000;0];

G = [0;2.50000000000000;0;0;0];

dxdt = A*x + B*(Q-Q_SS) + G*(C_A0-C_A0SS);
end
