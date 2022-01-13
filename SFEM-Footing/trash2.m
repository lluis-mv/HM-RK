function [] = trash2();

close all


X = [0,0;
    0.5, 0;
    1, 0;
    0, 0.5;
    0.5, 0.5;
    0,1];


plot(X(:,1), X(:,2), 'k*')
axis equal
hold on

plot([0, 1], [0.25, 0.25], 'k-.')
plot(  [0, 1], [0.75, 0.75], 'k-.')
plot( [0.25, 0.25], [0, 1], 'k-.')
plot( [0.75, 0.75], [0, 1], 'k-.')
for i = 1:1000
    x = rand(); y = rand();
    [number, color] = GetNumber(X, x,y);
    if ( number > 0)
        plot(x, y, [color, '*']);
    end
end



function [number, color] = GetNumber(X, x, y)
number = -1;
color = '';
if ( y > 1-x)
    return;
end

norma = zeros(6,1);
d = X-ones(6,1)*[x,y];
for i = 1:6
    norma(i) = norm(d(i,:));
end
number = find(norma == min(norma));

if ( number == 1)
    color = 'r';
elseif (number == 2)
    color  = 'g';
elseif (number == 3)
    color  = 'b';
elseif (number == 4)
    color  = 'c';
elseif (number == 5)
    color  = 'm';
elseif (number == 6)
    color  = 'y';
end


