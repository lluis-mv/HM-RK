syms x real


y = (x-0)*(x-1);
diff(diff(y, x),x)

syms h real


N = [h-x, x]/h

int( diff(N, x)'*diff(N,x), x, 0, h)

int(N, x, 0, h)



syms h1 real
syms h2 real



N1 = x/h1;
N2 = (h2-x)/h2


int(N1, x, 0, h1)
int(N2, x, 0, h2)


int( diff(N1, x), x, 0, h1)
int( diff(N2, x), x, 0, h2)