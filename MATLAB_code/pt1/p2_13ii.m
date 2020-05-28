clear; close all force; clc; 

%%%constants
V_1 = 4;               %volume [m^3]
V_2 = V_1;
V_3 = V_1;
V_4 = V_1;

v_0 = 10;              % volumetric flow ate [m^3/h]
k_0 = 1.5*10^3;        % reaction rate constant[1/h]
E = 2*10^3;            % enthalpy [kcal/kmol]
R = 1.987;             % mass-specific gas constant [kcal/kmole/K]
rho = 1000;            % density [kg/m^3]
Cp = 1;                % heat capacity with const. P [kcal/kg/K]
T_0 = 298;             % initial temp. [K]
neg_del_H = -2*10^5;   % negative delta enthalpy - realeased heat [kcal/kmole]
U = 550;               % [kcal/hr/K/m^2]
a = 50;                % area [m^2]

dt = 0.1;              % sample time [h]

%%%SSOP constants
Q_ss = 2.845*10^6;   %heat (SSOP) [kcal/h]
C_A0_ss = 1;         % inlet conc. (SSOP) [kmole/m^3]

%%%Steady-state Calculations --> myfun.m
x0 = [350; 0.5; 0.5; 350; 350; 350];

[x,fval] = fsolve(@myfun,x0) % x(1)=T1 x(2)=Ca x(3)=Cb x(4)=T3 x(5)=T4 x(6)=T2

%%%Linearize the nonlinear model around the SSOP; report A, B, and G --> lin1.m
A = [-9.3750         0         0         0    6.8750;
         0   -2.5000         0         0         0;
         0         0   -2.5000         0         0;
         0         0         0   -2.5000         0;
    6.8750         0         0    2.5000   -9.3750];

B =     [0;
         0;
         0;
    0.00025;
         0];

G = [0;
     2.5;
     0;
     0;
     0];
 

