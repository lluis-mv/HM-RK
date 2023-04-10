
syms k positive
syms Mbulk positive
syms AlphaStab positive
syms he positive
syms dt positive

syms x real
nNodes = 6;


% AlphaStab = 1;
N = [(he-x)/he, (x)/he]'
Me = int(N*N', x, 0, he)
Mse = AlphaStab*he/12*[1, -1; -1, 1];
gradN = diff(N, x);
He = int( gradN*gradN', x, 0, he);



% Construct elemental solution
M = sym(zeros(nNodes, nNodes));
H = sym(zeros(nNodes, nNodes));

for i = 1:nNodes-1

    index = [i, i+1];
    M(index, index) = M(index,index)+Me/Mbulk;
    M(index, index) = M(index,index)+Mse;
    H(index, index) = H(index,index)+He*k;
end
Helem = H;

A = M + dt*H;
Matrix = A;

pwOld = ones(nNodes,1);

f = M*pwOld; f(1) = 0;
A(1,:) = 0;
A(1,1) = 1;
pwNew = A\f;
pwNew = simplify(pwNew);

disp('ELEMENTAL')
for i = 1:nNodes-1
    simplify(solve(pwNew(i+1)-pwNew(i), AlphaStab))
end



% OK.
% Now nodal...
M = 0*M;
H = 0*H;
for i = 1:nNodes-1
    index = [i, i+1];
    M(index, index) = M(index,index)+Me/Mbulk;
    M(index, index) = M(index,index)+Mse;
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

for nod = 1:nNodes
    index = Nod(nod).NeigNodes;
    H(index, index) =H(index, index) + k*Nod(nod).gradNod*Nod(nod).gradNod'*Nod(nod).h;
end
Hnodal = H;

A = M + dt*H;
MatrixNodal = A;
M2 = M;

pwOld = ones(nNodes,1);

f = M*pwOld; f(1) = 0;
A(1,:) = 0;
A(1,1) = 1;
pwNodal = A\f;
pwNodal = simplify(pwNodal);

disp('NODAL')

for i = 1:nNodes-1
    simplify(solve(pwNodal(i+1)-pwNodal(i), AlphaStab))
end



return;
pp = [pwNew, pwNodal];
pp = subs(pp, Mbulk, 1)
pp = subs(pp, k, 1)
pp = subs(pp, AlphaStab, 0)
pp = subs(pp, he, 1/(nNodes-1));
pp = subs(pp, dt, 1E-12)
pp = eval(pp);
plot(pp);



