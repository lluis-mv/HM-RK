function [X, C] = ReadTheMesh(XFILE)
if (nargin == 0)
    XFILE = 'ThisMesh.msh';
end
fid = fopen(XFILE, 'r');

line = fgetl(fid);
line = fgetl(fid);


while (true)
    line = fgetl(fid);
    if ( line(1) == 'E')
        break;
    end
    line1 = str2num(line);
    X(line1(1),:) = line1(2:3);
end


line = fgetl(fid);
line = fgetl(fid);

while (true)
    line = fgetl(fid);
    if ( line(1) == 'E')
        break;
    end
    line1 = str2num(line);
    C(line1(1),:) = line1(2:9);
end
hola = 1;