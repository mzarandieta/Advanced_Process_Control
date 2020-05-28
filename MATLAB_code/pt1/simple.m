function [dydt] = simple(t,y)
dydt1=t^2;
dydt2=3*t;
dydt=[dydt1 ;dydt2]
end