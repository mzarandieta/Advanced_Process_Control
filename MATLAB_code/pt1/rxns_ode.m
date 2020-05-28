function [dsdt,qout] = rxns_ode(t,s)
%%%constants
V_1 = 4;               %volume [m^3]
V_3 = V_1;
V_4 = V_1;

v_0 = 10;              % volumetric flow ate [m^3/h]
k_0 = 1.5e3;        % reaction rate constant[1/h]
E = 2e3;            % activation energy for rxn [kcal/kmol]
R = 1.987;             % mass-specific gas constant [kcal/kmole/K]
rho = 1000;            % density [kg/m^3]
C_p = 1;                % heat capacity with const. P [kcal/kg/K]
T_0 = 298;             % initial temp. [K]
neg_del_H = -2e5;   % negative delta enthalpy - realeased heat [kcal/kmole]
U = 550;               % internal energy [kcal/hr/K/m^2]
a = 50;                % area [m^2]

%%%SSOP constants
Q = 2.845e6;   % manipulated measurment --> Q_ss heat (SSOP) [kcal/h]
C_A0 = 1;        % disturbance --> C_A0_ss inlet conc. (SSOP) [kmole/m^3]

if ( t<0 )
    Q = 2.845e6;   % manipulated measurment --> Q_ss heat (SSOP) [kcal/h]
    C_A0 = 1;        % disturbance --> C_A0_ss inlet conc. (SSOP) [kmole/m^3]
elseif ( t>0 && t<10)
    Q = 2.7e6;   % manipulated measurment --> Q_ss heat (SSOP) [kcal/h]
    C_A0 = 1.2;        % disturbance --> C_A0_ss inlet conc. (SSOP) [kmole/m^3]
elseif ( t > 10 && t < 20)
    Q = 3.1e6;   % manipulated measurment --> Q_ss heat (SSOP) [kcal/h]
    C_A0 = 1.2;        % disturbance --> C_A0_ss inlet conc. (SSOP) [kmole/m^3]
end

ds1dt = ( (v_0*rho*C_p*(T_0-s(1))) + (U*a*(s(5)-s(1))) ) / (rho*C_p*V_1);
ds2dt = ( (v_0*(C_A0-s(2))) - (k_0*exp(-((E)/(R*s(4))))*s(2)*V_3) ) / (V_3);
ds3dt = ( -(v_0*s(3)) + (k_0*exp(-((E)/(R*s(4))))*s(2)*V_3  )) / (V_3);
ds4dt = ( (v_0*rho*C_p*((s(1)+Q/(v_0*rho*C_p))-s(4))) + (neg_del_H*k_0*exp(-((E)/(R*s(4)))))*s(2)*V_3 ) / (rho*C_p*V_3);
ds5dt = ( (v_0*rho*C_p*(s(4)-s(5))) - (U*a*(s(5)-s(1))) ) / (rho*C_p*V_4);

dsdt=[ds1dt;ds2dt;ds3dt;ds4dt;ds5dt]; 

end