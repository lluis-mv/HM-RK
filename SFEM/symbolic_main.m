



syms x real
syms y real


N = [1-x-y, x, y]

subs( subs(N, x, 0), y, 0)
 
subs( subs(N, x, 1), y, 0)
  
subs( subs(N, x, 0), y, 1)

res = int(N, y, 0, 1-x)
int(res, x, 0, 1)






hold off

plot([0, 1, 0, 0], [0, 0, 1, 0], 'k', 'linewidth', 3)
axis equal
hold on
xx = [0, 1];
yy = 0.5*(1-xx);
plot(xx, yy, 'r')

xx = [0, 0.5];
yy = 1-2*xx;
plot(xx, yy, 'r')


xx = [0, 0.5];
yy = xx;
plot(xx, yy, 'r')





% Integration over Patch1
res = int( N, y, 0,  0.5*(1-x));
patch1 = int(res, x, 0, 1/3);

res = int(N, y, 0, 1-2*x);
patch1 = patch1 + int(res, x, 1/3, 1/2)


% Integration over patch2
res = int(N, y, 0, x);
patch2 = int(res, x, 1/3, 1/2);

res = int(N, y, 0, 1-x);
patch2 = patch2 + int(res, x, 0.5,1);

res = int(N, y, 0, 1-2*x);
patch2 = patch2 - int(res, x, 1/3, 0.5)

% integration over pathc3
res = int(N, y, 0, 1-x);
patch3 = int(res, x, 0,  0.5);

res = int(N, y, 0, 0.5*(1-x));
patch3 = patch3 - int(res, x, 0, 1/3);

res = int(N, y, 0, x);
patch3 = patch3 - int(res, x, 1/3, 1/2)

