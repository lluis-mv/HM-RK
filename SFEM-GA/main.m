% xi = linspace(-6,2,300);
% yi = linspace(-4,4,300);
% [X,Y] = meshgrid(xi,yi);
% Z = ps_example([X(:),Y(:)]);
% Z = reshape(Z,size(X));
% surf(X,Y,Z,'MeshStyle','none')
% colormap 'jet'
% view(-26,43)
% xlabel('x(1)')
% ylabel('x(2)')
% title('ps\_example(x)')
% 
% rng default % For reproducibility
% options = optimoptions('ga','PlotFcn', @gaplotbestf, 'Display','iter','PopulationSize',400);
% [x,fval,exitflag,output,population,scores] =  ga(@ps_example,2,options);


% AlphaVal = linspace(0,4,20);
% 
% i = 1;
% for alpha = AlphaVal
%     Norm(i) = this_function(alpha);
%     figure(232);
%     semilogy(AlphaVal(1:i), Norm(1:i), '*-.')
%     i = i+1;
% end
% 
% return;

figure(878); clf;
fun = @(x) this_function(x)
% [x,fval,exitflag,output,population,scores] =  ga(fun,1,[],[],[],[],[-0.1],[0.1])

X = fminsearch(fun,6)





