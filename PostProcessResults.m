
function [] = PostProcessResults(X, T, NodalData, GPInfo, time, WriteMesh, XFILE)


nNodes = size(X,1);


if (WriteMesh)
% write The Mesh
fid = fopen( ['a', XFILE, '.msh'], 'w' );

if ( size(T,2) == 3)
    fprintf(fid, 'MESH  dimension 2  ElemType Triangle Nnode 3 \n');
elseif ( size(T,2) == 6)
    fprintf(fid, 'MESH  dimension 2  ElemType Triangle Nnode 6 \n');
end

fprintf(fid, 'Coordinates \n ');
for i = 1:nNodes
    fprintf(fid, '%i %e %e %e \n', [i, X(i,:), 0]);
end
fprintf(fid, 'END Coordinates \n');

if ( size(T,2) == 3)
    ElemFormat = ['%i %i %i %i \n'];
elseif ( size(T,2) == 6)
    ElemFormat = ['%i %i %i %i %i %i %i \n'];
end

fprintf(fid, 'elements \n');
for i = 1:nElem
    fprintf( fid, ElemFormat, [i, T(i,:)]);
end
fprintf(fid, 'END elements \n');
fclose(fid);
end



fid = fopen(['a', XFILE, '.res'], 'a' );
fprintf(fid, 'Gid Post Results File 1.0 \n');



fprintf( fid, 'GaussPoints "GP"  ElemType Triangle  \n');
fprintf(fid, ' Number Of Gauss Points: 2 \n');
fprintf(fid, 'Natural Coordinates: Internal \n');
fprintf(fid, 'End GaussPoints \n');

time = num2str(rand());

fprintf( fid, ['Result "DISPLACEMENT" "HM-RK" ', time,' Vector OnNodes \n']);
fprintf( fid, 'ComponentNames "X-DISPLACEMENT", "Y-DISPLACEMENT", "Z-DISPLACEMENT" \n');
fprintf( fid, 'Values \n');
for i = 1:nNodes
    fprintf( fid, '%i %e %e %e \n', [i, U(i,:), 0] );
end
fprintf( fid, 'End Values \n');



fprintf(fid, ['Result "Small_strain_tensor" "HM-RK" ', time,'  Matrix OnGaussPoints "GP" \n ']);
fprintf(fid, ' values \n ');
for i = 1:nElem
    Celem = T(i,:);
    [S1, S2, S3, S4] = ComputeStrains( X(Celem, :), U(Celem,:));
    fprintf( fid, '%i %e %e %e %e %e %e \n', [i, S1] );
    fprintf( fid, '   %e %e %e %e %e %e \n', S2 );
    fprintf( fid, '   %e %e %e %e %e %e \n', S3 );
    fprintf( fid, '   %e %e %e %e %e %e \n', S4 );
end
fprintf(fid, ' END values \n');



fprintf(fid, ['Result "ALMANSI_STRAIN_TENSOR" "HM-RK" ', time,'  Matrix OnGaussPoints "GP" \n ']);
fprintf(fid, ' values \n ');
for i = 1:nElem
    Celem = T(i,:);
    [S1, S2, S3, S4] = ComputeAlmansi( X(Celem, :), U(Celem,:));
    fprintf( fid, '%i %e %e %e %e %e %e \n', [i, S1] );
    fprintf( fid, '   %e %e %e %e %e %e \n', S2 );
    fprintf( fid, '   %e %e %e %e %e %e \n', S3 );
    fprintf( fid, '   %e %e %e %e %e %e \n', S4 );
end
fprintf(fid, ' END values \n ');

fprintf(fid, ['Result "difference_strain_tensor" "HM-RK" ', time,'  Matrix OnGaussPoints "GP" \n ']);
fprintf(fid, ' values \n ');
for i = 1:nElem
    Celem = T(i,:);
    [S1, S2, S3, S4] = ComputeDifference( X(Celem, :), U(Celem,:));
    fprintf( fid, '%i %e %e %e %e %e %e \n', [i, S1] );
    fprintf( fid, '   %e %e %e %e %e %e \n', S2 );
    fprintf( fid, '   %e %e %e %e %e %e \n', S3 );
    fprintf( fid, '   %e %e %e %e %e %e \n', S4 );
end
fprintf(fid, ' END values \n ');


fclose(fid);

%system(['mv a', XFILE, '.res ',  XFILE, '.res ']);
movefile(['a', XFILE, '.msh '],  [XFILE, '.msh ']);
movefile(['a', XFILE, '.res '],  [XFILE, '.res ']);
function [S1, S2, S3, S4] = ComputeDifference( xe, ue)


A = StrainsDifference( xe(1,1), xe(2,1), xe(3,1), xe(4,1), ...
    xe(1,2), xe(2,2), xe(3,2), xe(4,2), ...
    ue(1,1), ue(2,1), ue(3,1), ue(4,1), ...
    ue(1,2), ue(2,2), ue(3,2), ue(4,2));

S1 = A(:,1)';
S2 = A(:,2)';
S3 = A(:,3)';
S4 = A(:,4)';


function [S1, S2, S3, S4] = ComputeAlmansi( xe, ue)


A = StrainsAlmansi( xe(1,1), xe(2,1), xe(3,1), xe(4,1), ...
    xe(1,2), xe(2,2), xe(3,2), xe(4,2), ...
    ue(1,1), ue(2,1), ue(3,1), ue(4,1), ...
    ue(1,2), ue(2,2), ue(3,2), ue(4,2));

S1 = A(:,1)';
S2 = A(:,2)';
S3 = A(:,3)';
S4 = A(:,4)';


function [S1, S2, S3, S4] = ComputeStrains( xe, ue)


A = StrainsEpsilon( xe(1,1), xe(2,1), xe(3,1), xe(4,1), ...
    xe(1,2), xe(2,2), xe(3,2), xe(4,2), ...
    ue(1,1), ue(2,1), ue(3,1), ue(4,1), ...
    ue(1,2), ue(2,2), ue(3,2), ue(4,2));

S1 = A(:,1)';
S2 = A(:,2)';
S3 = A(:,3)';
S4 = A(:,4)';


function [X, T ] = ConvertToFEMMesh(rr, zz, U)


nNodes = size(rr, 1)*size(rr, 2);
X = zeros(nNodes, 2);


for i = 1:size(rr,1)
    for j = 1:size(rr,2)
        iNod = i + (j-1)*size(rr,1);
        X(iNod, 1) = rr(i,j);
        X(iNod, 2) = zz(i,j);
    end
end

T = [];

for i = 1:size(rr,1)-1
    for j = 1:size(rr,2)-1
        nNod1 = i + (j-1)*size(rr,1);
        nNod2 = i +(j)*size(rr,1);
        nNod3 = i +1 +(j)*size(rr,1);
        nNod4 = i +1 +(j-1)*size(rr,1);
        
        nNod = [nNod1, nNod2, nNod3, nNod4, nNod1]';
        
        xNod = X(nNod, :);
        
        XNod =  X(nNod, :) - U(nNod,:);
        
        a = polyarea(xNod(:,1), xNod(:,2));
        A = polyarea(XNod(:,1), XNod(:,2));
        
        ratio = abs( (A-a) / A);
        
        
        if (  ratio > 0.25 )
            % do nothing
        else
            T = [T;
                nNod1, nNod2, nNod3, nNod4];
        end
    end
end


function U  = ConvertToVector( u1, u2)

nNodes = size(u1, 1)*size(u1, 2);
U = zeros(nNodes, 2);

for i = 1:size(u1,1)
    for j = 1:size(u1,2)
        nNod = i + (j-1)*size(u1,1);
        U(nNod, 1) = u1(i,j);
        U(nNod, 2) = u2(i,j);
    end
end
