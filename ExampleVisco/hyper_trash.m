function [] = hyper_trash()



x0 = [2,-1,0,0]';
fun = @(t,x) sourceTerm(t,x);

[t,y]= ode45(fun, [0, 2*pi], x0);
figure(88); clf;
plot(t,y)

function [y] = sourceTerm(t, x)
t1_2=x(3);
t2_1=x(4);
t1_1= x(1);
t2_2 = x(2);
y(1,1)=  t1_2 + t2_1;
y(2,1) = - t1_2 - t2_1;
y(3,1) =    t2_2 - t1_1;
y(4,1)=  t2_2 - t1_1;