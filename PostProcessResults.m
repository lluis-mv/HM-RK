
function [] = PostProcessResults(X, T, NodalData, GPInfo, time, WriteMesh, XFILE)


nNodes = size(X,1);
nElem = size(T,1);

if (WriteMesh)
    % write The Mesh
    fid = fopen( [XFILE, '.msh'], 'w' );
    
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
    
    
    fid = fopen([XFILE, '.res'], 'w' );
    fprintf(fid, 'Gid Post Results File 1.0 \n');
    fprintf( fid, 'GaussPoints "GP"  ElemType Triangle  \n');
    if ( size(T,2) == 3)
        fprintf(fid, ' Number Of Gauss Points: 1 \n');
        fprintf(fid, 'Natural Coordinates: Internal \n');
        fprintf(fid, 'End GaussPoints \n');
    elseif ( size(T,2) == 6)
        fprintf(fid, ' Number Of Gauss Points: 3 \n');
        fprintf(fid, 'Natural Coordinates: Given \n');
        fprintf( fid, '%e %e \n', [2/3, 1/6] );
        fprintf( fid, '%e %e \n', [1/6, 1/6] );
        fprintf( fid, '%e %e \n', [1/6, 2/3] );
        fprintf(fid, 'End GaussPoints \n');
    end
    fclose(fid);
end



fid = fopen([XFILE, '.res'], 'a' );






time = num2str(time);

U = zeros(nNodes,2);
WP = zeros(nNodes,1);
for i = 1:nNodes
    U(i,1) = NodalData( 3*(i-1)+1);
    U(i,2) = NodalData( 3*(i-1)+2);
    WP(i) = NodalData(3*(i-1)+3);
end

fprintf( fid, ['Result "DISPLACEMENT" "HM-RK" ', time,' Vector OnNodes \n']);
fprintf( fid, 'ComponentNames "X-DISPLACEMENT", "Y-DISPLACEMENT", "Z-DISPLACEMENT" \n');
fprintf( fid, 'Values \n');
for i = 1:nNodes
    fprintf( fid, '%i %e %e %e \n', [i, U(i,:), 0] );
end
fprintf( fid, 'End Values \n');



fprintf( fid, ['Result "WATER_PRESSURE" "HM-RK" ', time,' Scalar OnNodes \n']);
fprintf( fid, 'Values \n');
for i = 1:nNodes
    fprintf( fid, '%i %e \n', [i, WP(i)] );
end
fprintf( fid, 'End Values \n');


fprintf(fid, ['Result "Cauchy_stress_tensor" "HM-RK" ', time,'  Matrix OnGaussPoints "GP" \n ']);
fprintf(fid, ' values \n ');
for i = 1:nElem
    fprintf( fid, '%i %e %e %e %e %e %e \n', [i; GPInfo(i,1).StressNew] );
    for j = 2:size(GPInfo,2)
        fprintf( fid, '   %e %e %e %e %e %e \n', GPInfo(i,j).StressNew );
    end
end
fprintf(fid, ' END values \n ');



fclose(fid);

return;



fprintf(fid, ['Result "Cauchy_stress_tensor" "HM-RK" ', time,'  Matrix OnGaussPoints "GP" \n ']);
fprintf(fid, ' values \n ');





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
