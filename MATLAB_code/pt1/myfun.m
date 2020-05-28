function F = myfun(x)

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
Cp = 1;                % heat capacity with const. P [kcal/kg/K]
T_0 = 298;             % initial temp. [K]
neg_del_H = -2*10^5;   % negative delta enthalpy - realeased heat [kcal/kmole]
U = 550;               % internal energy [kcal/hr/K/m^2]
a = 50;                % area [m^2]

dt = 0.1;              % sample time [h]

%%%SSOP constants
Q_ss = 2.845*10^6;   % heat (SSOP) [kcal/h]
C_A0_ss = 1;         % inlet conc. (SSOP) [kmole/m^3]

F = [v_0*rho*Cp*(T_0-x(1))+U*a*(x(5)-x(1)); 
    (v_0*(C_A0_ss-x(2)))-(k_0*exp(-E/(R*x(4)))*V_3*x(2)); 
    (k_0*exp(-E/(R*x(4)))*x(2)*V_3)-(v_0*x(3)); 
    (neg_del_H*k_0*exp(-E/(R*x(4)))*x(2)*V_3)+(v_0*rho*Cp*(x(6)-x(4))); 
    (v_0*rho*Cp*(x(4)-x(5)))-(U*a*(x(5)-x(1))); 
    (x(1)+(Q_ss/(v_0*rho*Cp)))-x(6)];
end
