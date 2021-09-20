
syms k positive
syms Mbulk positive
syms AlphaStab positive
syms he positive
syms dt positive

syms x real

k = 1;
Mbulk = 1;
he = 0.5;
% dt = 1E0;
% AlphaStab = 0;


N = [(he-x)/he, (x)/he]';
gradN = diff(N, x);

K = int( gradN*gradN'*Mbulk, x, 0, he);
Q = -int( gradN*N', x, 0, he);
He = -int( k*gradN*gradN', x, 0, he);
Mse = AlphaStab*he/12*[1, -1; -1, 1];


% Construct elemental solution
nNodes = 3;

A1 = sym(zeros(2*nNodes, 2*nNodes));
A2 = A1;
for i = 1:nNodes-1
    ii = 2*(i-1);
    A1( ii+[1,3], ii+[1,3] ) = A1( ii+[1,3], ii+[1,3] ) + K;
    A1( ii+[1,3], ii+[2,4] ) = A1( ii+[1,3], ii+[2,4] ) + Q;
    A1( ii+[2,4], ii+[1,3] ) = A1( ii+[2,4], ii+[1,3] ) + Q';
    
    A1( ii+[2,4], ii+[2,4] ) = A1( ii+[2,4], ii+[2,4] ) + Mse;
    
    A2( ii+[2,4], ii+[2,4] ) = A2( ii+[2,4], ii+[2,4] ) + He;
end



Xn = zeros(2*nNodes,1);
f = Xn; 
f(1) = 1;


AA = A1+dt*A2;
f = A1*Xn+f;
f(2*nNodes-1) = 0;
AA(2,:) = 0; AA(2,2) = 1;
AA(2*nNodes-1,:) = 0; AA(2*nNodes-1,2*nNodes-1) = 1;
SystemMatrix = AA;
XNew = AA\f;
pwNew = XNew(2:2:end);
% % figure(1)
% % plot(pwNew);
% % figure(2)
% % plot(XNew(1:2:end));


for i = 1:nNodes-1
    solve(pwNew(i+1)-pwNew(i), AlphaStab);
    solve(pwNew(i+1)-1, AlphaStab)
end


C = [1:nNodes-1; 2:nNodes]';

for nod = 1:nNodes
    [NeigElements, trash] = find(C == nod);
    Nod(nod).NeigElement = NeigElements;
    Nod(nod).NeigNodes = unique( C(NeigElements,:));
    Nod(nod).h = he*sum(length(NeigElements))/2;
end

for nod = 1:nNodes
    gradNod = sym(zeros( length(Nod(nod).NeigNodes), 1));
    
    for neigElem = Nod(nod).NeigElement'
        Celem = C(neigElem,:);
        gradE = gradN*(he/2) / Nod(nod).h;
        
        ind = [];
        for jj = 1:2
            tt = find(Nod(nod).NeigNodes == Celem(jj));
            ind = [ind, tt];
        end
        gradNod(ind) = gradNod(ind)+gradE;
    end
    Nod(nod).gradNod = gradNod;
end


A1 = sym(zeros(2*nNodes, 2*nNodes));
A2 = A1;

for nod = 1:nNodes
    
   
    index = Nod(nod).NeigNodes;
   
    A1(2*index-1, 2*index-1) = A1(2*index-1, 2*index-1) + Mbulk*Nod(nod).gradNod*Nod(nod).gradNod'*Nod(nod).h;

    for el = [Nod(nod).NeigElement]'
        N = sym(zeros(1,length(Nod(nod).NeigNodes) ));
        for mm = 1:2
            ind = find( Nod(nod).NeigNodes == C(el,mm) );
            N(ind) = 1*he/8;
        end
        ind = find( nod == Nod(nod).NeigNodes);
        N(ind) = 3*he/8;
        
        Q = Nod(nod).gradNod*N;
        A1(2*index-1, 2*index) = A1(2*index-1, 2*index) - Q;
        A1(2*index, 2*index-1) = A1(2*index, 2*index-1) - Q';
    end
    
    
    A2(2*index, 2*index ) = A2( 2*index, 2*index ) - k*Nod(nod).gradNod*Nod(nod).gradNod'*Nod(nod).h;
    
end

for i = 1:nNodes-1
    ii = 2*(i-1);
    A1( ii+[2,4], ii+[2,4] ) = A1( ii+[2,4], ii+[2,4] ) + Mse;    
end



Xn = zeros(2*nNodes,1);
f = Xn; 
f(1) = 1;


AA = A1+dt*A2;
f = A1*Xn+f;
f(2*nNodes-1) = 0;
AA(2,:) = 0; AA(2,2) = 1;
AA(2*nNodes-1,:) = 0; AA(2*nNodes-1,2*nNodes-1) = 1;



XNew = AA\f;
pwNew = XNew(2:2:end);


for i = 1:nNodes-1
    solve(pwNew(i+1)-pwNew(i), AlphaStab);
    solve(pwNew(i+1)-1, AlphaStab)
    vpa(solve(pwNew(i+1)-1, AlphaStab))
    
end

return;


figure(1)
hold on
plot(pwNew);
hold off

figure(2)
hold on
plot(XNew(1:2:end));
hold off;




