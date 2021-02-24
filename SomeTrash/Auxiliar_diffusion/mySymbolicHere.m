clear all

syms y real


solutionpw = y^2*(y^1-1);


forceWP = diff( diff(solutionpw, y),y)


solutionU = 0.01*y^2*(y-1);

f = diff( diff(solutionU, y),y) - diff(solutionpw, y)



syms t real


u =  0.01*y^2*(y-1)*t;
pw = 0.1*y^2*(y-1);
pw = 0;

f = diff(diff(diff(u,y),y)-diff(pw, y), t)

f2 = diff( diff(u,y),t) - diff(pw, y)