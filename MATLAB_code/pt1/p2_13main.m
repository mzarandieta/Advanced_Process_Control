clear; close all force; clc; 

%%%constants
V_1 = 4;               %volume [m^3]
V_2 = V_1;
V_3 = V_1;
V_4 = V_1;

v_0 = 10;              % volumetric flow ate [m^3/h]
k_0 = 1.5*10^3;        % reaction rate constant[1/h]
E = 2*10^3;            % activation energy for rxn [kcal/kmol]
R = 1.987;             % mass-specific gas constant [kcal/kmole/K]
rho = 1000;            % density [kg/m^3]
C_p = 1;                % heat capacity with const. P [kcal/kg/K]
T_0 = 298;             % initial temp. [K]
neg_del_H = -2*10^5;   % negative delta enthalpy - realeased heat [kcal/kmole]
U = 550;               % internal energy [kcal/hr/K/m^2]
a = 50;                % area [m^2]

dt = 0.1;              % sample time [h]

%%%SSOP constants
Q_ss = 2.845*10^6;                 %heat (SSOP) [kcal/h]
C_A0_ss = 1;                       % inlet conc. (SSOP) [kmole/m^3]

%%from part ii --> myfun.m
T_2_ss = 646.72;                                  % [K]
s_ss = [362.22; 0.0154; 0.9846; 449.79; 385.57];  % s_ss = [T1 c_a c_b T3 T4]
q_ss = [449.79; 2.845*10^6];                      % qout = [s(4); Q_ss]
T_3= s_ss(4);
C_A= s_ss(2);



A =[ -(U*a + C_p*rho*v_0)/(C_p*V_1*rho),                                         0,        0,                                                                                  0,                (U*a)/(C_p*V_1*rho);
                                   0,      -(v_0 + V_3*k_0*exp(-E/(R*T_3)))/V_3,        0,                                             -(C_A*E*k_0*exp(-E/(R*T_3)))/(R*T_3^2),                                  0;
                                   0,                       k_0*exp(-E/(R*T_3)), -v_0/V_3,                                              (C_A*E*k_0*exp(-E/(R*T_3)))/(R*T_3^2),                                  0;
                             v_0/V_3, (k_0*neg_del_H*exp(-E/(R*T_3)))/(C_p*rho),        0, -(C_p*rho*v_0 - (C_A*E*V_3*k_0*neg_del_H*exp(-E/(R*T_3)))/(R*T_3^2))/(C_p*V_3*rho),                                  0;
                 (U*a)/(C_p*V_4*rho),                                         0,        0,                                                                            v_0/V_4, -(U*a + C_p*rho*v_0)/(C_p*V_4*rho)];

B = [0; 0; 0; 1/(rho*C_p*V_3); 0];

G = [0; v_0/V_3; 0; 0; 0];

%%%compare the linearized model with the original nonlinear model
%Call function simpleode45()

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
% 
% %%%compare the discrete-time simulation with the linear continuous-time simulation (from part iv) 
%       %I GIVE UP%