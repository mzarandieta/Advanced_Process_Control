close all 

tspan = [0 20];

%NON LINEAR 
%Non Linear Initial value = s_ss 

s_ss = [362.22; 0.0154; 0.9846; 449.79; 385.57];
[t1,y1] = ode45('rxns_ode', tspan, s_ss);


%LINEAR 
[t2,y2] = ode45('rxns_ode_lin', tspan, zeros(1,5));
titles=["T_1";"C_A";"C_B";"T_3";"T_4"];

%Plot all variables 

for n=1:5
    y2 (:,n) = y2(:,n) +s_ss(n);
%     figure()
  subplot(2,3,n);
    title(titles(n))
    hold on 
    plot(t1,y1(:,n),'r')
    plot(t2,y2(:,n),'b')
%     legend('nonlinear model','linearized model')
    xlabel('time [h]')
end
legend('nonlinear model','linearized model')
legend('Location','southeast')
legend('boxoff')

