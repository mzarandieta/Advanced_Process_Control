% Continuous-time Model
nx=3; nu=2;AA=[1  -3  0;  0  -2  0; 2  -2  -1];BB=[1  1;1  0;0  1];
% BB=[1  0;1  0;0  1];
%Convert to Discrete-time
dt=0.1; Ndt=200;  ddt=dt/Ndt;  sum=zeros(3);
for jjj=1:Ndt  
    sum=sum+expm(AA*jjj*ddt); 
end; 
AAd=expm(AA*dt), 
BBd=sum*BB*ddt, 
%part(i)
display('part(i)'); 
eig(AAd)
%part(ii)
display('part(ii)'); 
Lc=[BBd  AAd*BBd    AAd*AAd*BBd];
size(Lc), Q = orth(Lc) 
%part(iii)
yalmip('clear'); C1=[]; C2=[];
% Define Variables
Z0=sdpvar(nx); Z1=sdpvar(nu,nx);C1=[0<=Z0];
% Z0 > 0
% [ Z0-I             (Ad*Z0-Bd*Z1)]% [                   ] >0
% [(Ad*Z0-Bd*Z1)'              Z0]
C2=[ ( Z0-eye(nx)        (AAd*Z0-BBd*Z1)(AAd*Z0-BBd*Z1)'       Z0  )>=0];
Constraints=[C1,C2]; 
Objective=[];
% Define Constraints and Objective
options = sdpsettings('verbose',0,'solver','mosek');
sol = optimize(Constraints,Objective,options);
% Determine Feasibility
if sol.problem == 0   
    % problem is feasible
    Z0feas=value(Z0); 
    Z1feas=value(Z1); 
    LLfeas=Z1feas/Z0feas,eig(AAd-BBd*LLfeas)
elseif sol.problem == 1 
    % problem is infeasible
    display('Infeasible Problem');
else
    display('Hmm, something went wrong!'); 
    sol.info; pause; 
end
