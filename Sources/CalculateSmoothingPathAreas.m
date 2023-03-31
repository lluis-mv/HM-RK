
function GPElements = CalculateSmoothingPathAreas( Nodes, Elements, GPElements)

if ( size(Elements,2) == 3)
    for i =  1:size(GPElements,1)
        for j = 1:size(GPElements,2)
            GPElements(i,j).AllWeights = GPElements(i,j).Weight * 1/3*ones(1,3);
        end
    end

elseif ( size(Elements,2) == 4)

    for i =  1:size(GPElements,1)

        X = Nodes( Elements(i,:),:);

        [al, be, w]  = GetWeights(4);

        Areas = zeros(1,4);
        for cas = [1:4]
            [al, be, w]  = GetWeights(4);
            if ( cas == 1)
                al = al/2-0.5;
                be = be/2-0.5;
            elseif (cas == 2)
                al = al/2+0.5;
                be = be/2-0.5;
            elseif (cas == 3)
                al = al/2+0.5;
                be = be/2+0.5;
            elseif (cas == 4)
                al = al/2-0.5;
                be = be/2+0.5;
            end

            for ww = 1:length(w)
                alfa = al(ww);
                beta = be(ww);


                Nsmall_chi = [   beta/4 - 1/4,   alfa/4 - 1/4;
                    1/4 - beta/4, - alfa/4 - 1/4;
                    beta/4 + 1/4,   alfa/4 + 1/4;
                    - beta/4 - 1/4,   1/4 - alfa/4];


                J = Nsmall_chi'*X;

                Areas(cas) = Areas(cas)+det(J)*w(ww);
            end
        end

        for j = 1:size(GPElements,2)
            GPElements(i,j).AllWeights = Areas;
        end
    end
end
