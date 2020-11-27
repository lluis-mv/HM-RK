clear all

syms y real


solutionpw = y^2*(y^1-1);


forceWP = diff( diff(solutionpw, y),y)


solutionU = 0.01*y^2*(y-1);

f = diff( diff(solutionU, y),y) - diff(solutionpw, y)



syms t real


u =  0.01*y^2*(y-1)^2*t

pw = y*(y-1)*t;


ff = diff(diff(diff(u,y),y)-diff(pw, y), t)

ff2 = diff( diff(u,y),t) - diff(diff(pw, y), y)