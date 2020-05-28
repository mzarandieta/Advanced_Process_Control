%Linearization of model around SSOP, Report A, B and G 
%dx/dt= A*X+ B*u+ G*w

function [A,B,G] = lin1()

syms v_0 rho C_p T_0 T_1 T_3 T_4 U a V_1 V_3 V_4 C_A0 C_A C_B k_0 E R Q neg_del_H

f1 = ( (v_0*rho*C_p*(T_0-T_1)) + (U*a*(T_4-T_1)) ) / (rho*C_p*V_1);
f2 = ( (v_0*(C_A0-C_A)) - (k_0*exp(-((E)/(R*T_3)))*C_A*V_3) ) / (V_3);
f3 = ( -(v_0*C_B) + (k_0*exp(-((E)/(R*T_3)))*C_A*V_3  )) / (V_3);
f4 = ( (v_0*rho*C_p*((T_1+Q/(v_0*rho*C_p))-T_3)) + (neg_del_H*k_0*exp(-((E)/(R*T_3))))*C_A*V_3 ) / (rho*C_p*V_3);
f5 = ( (v_0*rho*C_p*(T_3-T_4)) - (U*a*(T_4-T_1)) ) / (rho*C_p*V_4);

A = jacobian ([f1; f2; f3; f4; f5], [T_1 C_A C_B T_3 T_4])

B = jacobian ([f1; f2; f3; f4; f5], [Q])

G = jacobian ([f1; f2; f3; f4; f5], [C_A0])

%%%constants
V_1 = 4;            %volume [m^3]
V_3 = V_1;
V_4 = V_1;

v_0 = 10;           % volumetric flow ate [m^3/h]
k_0 = 1.5e3;        % reaction rate constant[1/h]
E = 2e3;            % activation energy for rxn [kcal/kmol]
R = 1.987;          % mass-specific gas constant [kcal/kmole/K]
rho = 1000;         % density [kg/m^3]
C_p = 1;            % heat capacity with const. P [kcal/kg/K]
T_0 = 298;          % initial temp. [K]
neg_del_H = -2e5;   % negative delta enthalpy - realeased heat [kcal/kmole]
U = 550;            % internal energy [kcal/hr/K/m^2]
a = 50;             % area [m^2]

s_ss = [362.22; 0.0154; 0.9846; 449.79; 385.57];  % s_ss = [T_1 C_A C_B T_3 T_4]
C_A = s_ss(2);
T_3 = s_ss(4);


A = eval(A);

B = eval(B);

G = eval(G);

end