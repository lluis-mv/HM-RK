clear all

syms y real
syms t real


u =  0.1*y^2*(y-1)^2*t^2
pw = 0.01*y*(y-1)*t^2;



ff = diff(diff(diff(u,y),y) - diff(pw, y), t)

ff2 = +0*diff( diff(u,y),t) - diff(diff(pw, y), y)