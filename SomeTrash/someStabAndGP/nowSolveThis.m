function [] = nowSolveThis
figure(1)
hold off

R1 = [];
R2 = [];

x0 = 1/3*ones(18,1);
% x0([1:3,7:18]) = 0;
x0 = x0 + 0.1*rand(18,1)
x0 = 5*(rand(18,1)-0.5);
equations = 1:36;
x0 = [ 2/5, -1/5, -1/5, 3/5, -1/5, 3/5 -3/5, 3/5, 0, 0, 4/5, -4/5 -3/5, 0, 3/5, -4/5, 4/5, 0]';
 
x0 = x0+0.03*rand(18,1)
x = x0;
iter = 0;
for i = 1:100
    
    res = residual(x);
    
    thisResidual = res(equations);
    
    normResidual = norm(res);
    normThisResidual = norm(thisResidual);
    
    if ( normThisResidual < 1E-12 && iter > 1)
        break;
    end
    
    R1 = [R1, normResidual];
    R2 = [R2, normThisResidual];
    figure(1)
    semilogy(0:iter, R1, 'b*-.')
    hold on
    semilogy(0:iter, R2, 'r*-.')
    
    J = jacobian(x);
    
    J = J(equations,:);
    dX = -J\thisResidual;
    
    alfa = linspace(-1,3);
    for j = 1:length(alfa)
        xx(j) = norm(residual(x + alfa(j)*dX));
    end
    pol = polyfit( alfa, xx, 2);
    
    
    al = linspace(-1,2,100);
    rr = polyval(pol, al);
    
    figure(900)
    plot(alfa, xx, '*-.')
    hold on
    plot(al, rr)
%     hold off;
     
    a = pol(1); b = pol(2); c = pol(3);
    plot(al, a*al.^2 + b*al+c, 'k')
    hold off
    
    m = -b/(2*a)
    
    [r, ind] = min(xx);
    m = alfa(ind);
    m = 1;
    x = x+m*dX;
    iter = iter+1;
end

al = linspace(0,1,50)
be = linspace(0,1,50)
[al, be] = meshgrid(al, be);
T = delaunay(al,be);
%trimesh(T,x,y,z)
for i = 1:6
    a = x(1+i);
    b = x(6+i);
    c = x(12+i);
    na = a+b*al+c*be;
    
    [rs] = find( be > 1-al -0.001);
    na(rs) = nan;
    
    figure(232+i)
    trimesh(T, al, be, na);
end
    


return
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
x5 = x(5);
x6 = x(6);
x7 = x(7);
x8 = x(8);
x9 = x(9);
x10 = x(10);
x11 = x(11);
x12 = x(12);
x13 = x(13);
x14 = x(14);
x15 = x(15);
x16 = x(16);
x17 = x(17);
x18 = x(18);