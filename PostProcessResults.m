
function [] = PostProcessResults(HydroMechanical, X, T, NodalData, GPInfo, time, WriteMesh, XFILE)


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
    if ( size(GPInfo,2) == 1)
        fprintf(fid, ' Number Of Gauss Points: 1 \n');
        fprintf(fid, 'Natural Coordinates: Internal \n');
        fprintf(fid, 'End GaussPoints \n');
    elseif ( size(GPInfo,2) == 3)
        fprintf(fid, ' Number Of Gauss Points: 3 \n');
        fprintf(fid, 'Natural Coordinates: Given \n');
        [al, be] = GetWeights(3);
        for i = 1:3
            fprintf( fid, '%e %e \n', [al(i), be(i)] );
        end
        fprintf(fid, 'End GaussPoints \n');
    elseif ( size(GPInfo, 2) == 6)
        [al, be] = GetWeights(6);
        fprintf(fid, ' Number Of Gauss Points: 6 \n');
        fprintf(fid, 'Natural Coordinates: Given \n');
        for i = 1:6
            fprintf( fid, '%e %e \n', [al(i), be(i)] );
        end
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


if (HydroMechanical)
    fprintf( fid, ['Result "WATER_PRESSURE" "HM-RK" ', time,' Scalar OnNodes \n']);
else
    fprintf( fid, ['Result "PRESSURE" "HM-RK" ', time,' Scalar OnNodes \n']);
end
fprintf( fid, 'Values \n');

for i = 1:nNodes
    fprintf( fid, '%i %e \n', [i, WP(i)] );
end
fprintf( fid, 'End Values \n');

if (HydroMechanical)
    fprintf(fid, ['Result "Cauchy_stress_tensor" "HM-RK" ', time,'  Matrix OnGaussPoints "GP" \n ']);
else
    fprintf(fid, ['Result "total_Cauchy_stress_tensor" "HM-RK" ', time,'  Matrix OnGaussPoints "GP" \n ']);
end
fprintf(fid, ' values \n ');
Idev = eye(6)-1/3*[1,1,1,0,0,0]'*[1,1,1,0,0,0];
m = [1,1,1,0,0,0]';
for i = 1:nElem
    fprintf( fid, '%i ', [i]);
    for j = 1:size(GPInfo,2)
        stress = Idev*GPInfo(i,j).StressNew + m*GPInfo(i,j).N*NodalData(GPInfo(i,j).dofsWP);
        fprintf( fid, '   %e %e %e %e %e %e \n', stress );
    end
end
fprintf(fid, ' END values \n ');

if (HydroMechanical)
    fprintf(fid, ['Result "Cauchy_stress_tensor" "HM-RK" ', time,'  Matrix OnGaussPoints "GP" \n ']);
    fprintf(fid, ' values \n ');
    for i = 1:nElem
        fprintf( fid, '%i ', [i]);
        for j = 1:size(GPInfo,2)
            stress = GPInfo(i,j).StressNew;
            fprintf( fid, '   %e %e %e %e %e %e \n', stress );
        end
    end
    fprintf(fid, ' END values \n ');
end



fprintf(fid, ['Result "Strain_tensor" "HM-RK" ', time,'  Matrix OnGaussPoints "GP" \n ']);
fprintf(fid, ' values \n ');
for i = 1:nElem
    fprintf( fid, '%i %e %e %e %e %e %e \n', [i; GPInfo(i,1).StrainNew] );
    for j = 2:size(GPInfo,2)
        fprintf( fid, '   %e %e %e %e %e %e \n', GPInfo(i,j).StrainNew );
    end
end
fprintf(fid, ' END values \n ');

if ( length(GPInfo(1,1).HistoryNew) > 0)
    fprintf(fid, ['Result "Preconsolidation" "HM-RK" ', time,'  Scalar OnGaussPoints "GP" \n ']);
    fprintf(fid, ' values \n ');
    for i = 1:nElem
        fprintf( fid, '%i %e  \n', [i; GPInfo(i,1).HistoryNew(1)] );
        for j = 2:size(GPInfo,2)
            fprintf( fid, '   %e  \n', GPInfo(i,j).HistoryNew(1) );
        end
    end
    fprintf(fid, ' END values \n ');
end


fprintf(fid, ['Result "VonMises" "HM-RK" ', time,'  Scalar OnGaussPoints "GP" \n ']);
fprintf(fid, ' values \n ');
for i = 1:nElem
    fprintf( fid, '%i ' , i );
    for j = 1:size(GPInfo,2)
        stress = Idev*GPInfo(i,j).StressNew;
        J2 = sqrt(3/2)*sqrt( stress(1)^2+ stress(2)^2 + stress(3)^2 + 2*stress(4)^2);
        fprintf( fid, '   %e  \n', J2);
    end
end
fprintf(fid, ' END values \n ');



fclose(fid);

